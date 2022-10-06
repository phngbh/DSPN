#!/usr/bin/env Rscript

dat = read.csv("/data/20180212_01_KORA.F4_metabolon_komplett_changedIDs.csv", check.names = F, sep = ";", row.names = 1)
meta_anno = read.csv("20120328_02_KORA.S4_F4_annotation_Mnumber_pathways.csv", sep = ";")
clin <- read_sas("/data/gesamt_k05720g_v2007.sas7bdat", NULL) #%>% filter(!is.na(utmnsi))
info = data.frame(ID = rownames(dat), 
                  mnsi = clin$utmnsi[match(rownames(dat), clin$un_metabMetabolon_f4)]) #%>% na.omit()
#dat = dat[info$ID, ] 
keep = apply(dat, 2, function(x) sum(!is.na(x)) >= nrow(dat)/100*30)
dat = dat[,keep]
dat.woc = mutate(dat, outcome = rep(1, times = nrow(dat)))
met = colnames(dat)
source(UnMetImp.R)
dat.imp.knn = UnMetImp(DataFrame = dat.woc, imp_type = 'knn', 
                       group1 = met, outcome = "outcome", use_covars = FALSE , logScale = F)
dat.imp = dat.imp.knn$mids %>% select(-outcome)
mah.dist = mahalanobis(dat.imp,
                       center = colMeans(dat.imp),
                       cov = pseudoinverse(cov(dat.imp)),
                       inverted = T)
mah.dist.std = abs((mah.dist - mean(mah.dist))/sd(mah.dist))
keep = names(mah.dist.std[mah.dist.std<=4])
dat.imp = dat.imp[keep,] %>% t()
info = filter(info, ID %in% colnames(dat.imp))
saveRDS(dat.imp, file = "metabolomics_processed.rds")