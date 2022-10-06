#!/usr/bin/env Rscript

dat = read.csv("/data/KORA_F4_993_samples_Expressionsdaten.csv",header = T,
                    check.names = F,row.names = 1) %>% as.matrix()
gene.an = read.csv("Annotation_HumanHT-12v3_final.csv",check.names = F,sep = "\t")
gene.an = gene.an %>%
  mutate(symbol = unlist(mget(gene.an$Probe_Id,illuminaHumanv3SYMBOL)),
         EntrezID = unlist(mget(gene.an$Probe_Id,illuminaHumanv3ENTREZID)),
         gene = unlist(mget(gene.an$Probe_Id,illuminaHumanv3GENENAME)))
# gene.an$symbol = ifelse(is.na(gene.an$symbol), gene.an$ILMN_Gene, gene.an$symbol)
keep = gene.an$Probe_Id[grep("good", gene.an$`QC comment`)]
sex_tr = gene.an$Probe_Id[gene.an$ILMN_CHR %in% c("X","Y")]
clin <- read_sas("gesamt_k05720g_v2007.sas7bdat", NULL)
tech = read.table("/data/Technical_variables/Transcriptomics/KORA_F4_technical_variables.txt", sep = "", header = T)
info = tech %>% mutate(mnsi = clin$utmnsi[match(.$ZZ.NR,clin$un_expr_f4ogtt)])
dat.fil = dat[ !rownames(dat) %in% sex_tr & rownames(dat) %in% keep
               ,as.character(info$ZZ.NR)]
#symbol = gene.an$symbol[match(rownames(dat.fil),gene.an$Probe_Id)]
# dat.fil = avereps(dat.fil,symbol)
# str(dat.fil)
#design = model.matrix(~ 0 + mnsi, data=info)
dat.nobatch = removeBatchEffect(x=dat.fil, 
                                #design = design,
                                batch = info$p_amplification,
                                covariates = info[,c("RIN","sample_storage_time")])

meansd = data.frame(mean = apply(dat.nobatch,1,mean), sd = apply(dat.nobatch,1,sd), 
                    iqr = apply(dat.nobatch,1,IQR))
ggplot(meansd, aes(x=mean,y=sd))+
  geom_point(col = "black",alpha=0.5)+
  stat_smooth(method = "loess",se=F,col="red")

saveRDS(dat, file = "transcriptomics_processed.rds")