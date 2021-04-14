library(tidyverse)
library(haven)
library(corpcor)
source("UnMetImp.R")

### Load data #####

dat = read.csv("20180212_01_KORA.F4_metabolon_komplett_changedIDs.csv", check.names = F, sep = ";", row.names = 1)
met = colnames(dat)
clin <- read_sas("gesamt_k05720g_v2007.sas7bdat", NULL) #%>% filter(!is.na(utmnsi))
info = data.frame(ID = rownames(dat), 
                  mnsi = clin$utmnsi[match(rownames(dat), clin$un_metabMetabolon_f4)]) #%>% na.omit()
#dat = dat[info$ID, ] 
#keep = apply(dat, 2, function(x) sum(!is.na(x)) >= nrow(dat)/100*10)
dat.woc = mutate(dat, outcome = rep(1, times = nrow(dat)))
dat.imp.knn = UnMetImp(DataFrame = dat.woc, imp_type = 'knn', 
                       group1 = met, outcome = "outcome", use_covars = FALSE , logScale = F)
dat.imp = dat.imp.knn$mids %>% select(-outcome)
mah.dist = mahalanobis(dat.imp,
                       center = colMeans(dat.imp),
                       cov = pseudoinverse(cov(dat.imp)),
                       inverted = T)
mah.dist.std = abs((mah.dist - mean(mah.dist))/sd(mah.dist))
keep = names(mah.dist.std[mah.dist.std<=4])
dat.imp = dat.imp[keep,]
info = filter(info, ID %in% rownames(dat.imp))

### Visualize ####

meansd = data.frame(mean = apply(dat.imp,2, mean), 
                    sd = apply(dat.imp,2,sd), 
                    iqr = apply(dat.imp,2,IQR), 
                    med = apply(dat.imp, 2, median))
ggplot(meansd, aes(x=mean,y=sd))+
  geom_point(col = "black",alpha=0.5)+
  stat_smooth(method = "lm",se=F,col="red")

pca = prcomp(dat.imp)
pca_df = data.frame(pca$x, info)
pca_var = summary(pca)$importance
varp = round(pca_var[2,]*100,2)
ggplot(pca_df,aes(x=PC1,y=PC2)) + 
  geom_point(size = 3, alpha = 0.5)+
  xlab(paste("PC1 (",varp[1],"%)")) +
  ylab(paste("PC2 (",varp[2],"%)")) +
  theme(panel.background = element_rect(fill = "white"))

saveRDS(t(dat.imp), file = "metabolomics_processed.rds")
