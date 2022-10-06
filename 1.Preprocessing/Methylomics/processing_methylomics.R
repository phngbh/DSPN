#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

library(ggplot2, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(haven, quietly = TRUE)
library(ChAMPdata, quietly = TRUE)
library(ChAMP, quietly = TRUE)
library(limma, quietly = TRUE)
library(ggrepel)
library(cowplot)
library(ebGSEA, lib.loc = "Rlib")
library(fgsea, lib.loc = "Rlib")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19, lib.loc = "Rlib")
# library(qqman, lib.loc = "Rlib")
# library(scales)
# library(scattermore, lib.loc = "Rlib")
# library(ggfastman, lib.loc = "Rlib")
library(tidyverse)
library(caret, lib.loc = "Rlib")
library(glmnet)
library(pROC, lib.loc = "Rlib")
library(globaltest)
library(parallel)
library(reactome.db, lib.loc = "Rlib")
library(org.Hs.eg.db)

#####Load data#######

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
#mvalue = readRDS("mvalue.rds")
beta = loadRData("/data/KF4_beta_qn_bmiq.RData")
detP = loadRData("/data/KF4_dp_id.RData")
# detP = detP[rownames(beta), colnames(beta)]
batch = loadRData("/data//Technical_variables/Methylomics/KF4_PlateChip_1727.RData")
pcs = loadRData("/data/Technical_variables/Methylomics/control_probe_pcs_n1727.RData")
pcs.df = as.data.frame(pcs[,1:20]) %>% mutate(ID = rownames(.))
batch$ZZ_nr = as.character(batch$ZZ_nr)
batch.pcs = inner_join(batch, pcs.df, by = c("ZZ_nr" = "ID"))
cell.type = read.csv("/data/Technical_variables/Methylomics/KF4_QN_BMIQ_estimated_cell_distribution_meanimpute473_lessThanOneFALSE.csv", sep = ";")
cell.type$ZZ_nr = as.character(cell.type$ZZ_nr)
info = left_join(batch.pcs, cell.type) %>% column_to_rownames("ZZ_nr")
info$Chip = as.factor(info$Chip)
info$Batch = as.factor(info$Batch)
info = rownames_to_column(info, "Sample_Name")
beta = beta[,info$Sample_Name]
all(info$Sample_Name == colnames(beta))
detP = detP[rownames(beta), colnames(beta)]

####CpG site and sample filtering#####

beta.fil = champ.filter(beta = beta, pd = info, detP = detP, autoimpute = T, ProbeCutoff = 0.05)

#####Technical effect correction######

champ.SVD(beta = beta.fil$beta, pd = beta.fil$pd)

beta.fil.df = beta.fil$beta
meansd = data.frame(mean = apply(beta.fil.df,1,mean), sd = apply(beta.fil.df,1,sd), 
                    iqr = apply(beta.fil.df,1,IQR))
beta.topvar = beta.fil.df[rownames(meansd[meansd$iqr>quantile(meansd$iqr,0.9),]),]
pca = prcomp(t(beta.topvar), scale. = T)
pca_df = data.frame(pca$x, info)
pca_var = summary(pca)$importance
varp = round(pca_var[2,]*100,2)
ggplot(pca_df,aes(x=PC1,y=PC2,col = Plate)) + 
  geom_point(size = 3, alpha = 0.5)+
  xlab(paste("PC1 (",varp[1],"%)")) +
  ylab(paste("PC2 (",varp[2],"%)")) +
  theme(panel.background = element_rect(fill = "white"))

covariates = info[,c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")]
beta.nobatch = removeBatchEffect(x = log(beta.fil.df),
                                   batch = info[,"Plate"],
                                   batch2 = info[,"Chip"],
                                 covariates = covariates) %>% exp()
champ.SVD(beta = beta.nobatch, pd = beta.fil$pd) #Double-check the technical effects after correction
saveRDS(beta.nobatch, file = "methylomics_processed.rds")