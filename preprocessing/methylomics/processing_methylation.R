#!/usr/bin/env Rscript

library(ggplot2, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(haven, quietly = TRUE)
library(ChAMPdata, quietly = TRUE)
library(ChAMP, quietly = TRUE)
library(limma, quietly = TRUE)

#####Load data#######

# args = commandArgs(trailingOnly=TRUE)
# setwd(args[1])

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
load("./KF4_beta_qn_bmiq.RData")
detP = loadRData("./KF4_dp_id.RData")
detP = detP[rownames(beta), colnames(beta)]
batch = loadRData("./Technical_variables/KF4_PlateChip_1727.RData")
pcs = loadRData("./Technical_variables/control_probe_pcs_n1727.RData")
pcs.df = as.data.frame(pcs[,1:20]) %>% mutate(ID = rownames(.))
batch$ZZ_nr = as.character(batch$ZZ_nr)
batch.pcs = inner_join(batch, pcs.df, by = c("ZZ_nr" = "ID"))
clin <- read_sas("./gesamt_k05720g_v2007.sas7bdat", NULL) %>% filter(!is.na(utmnsi))
index = match(batch.pcs$ZZ_nr, clin$un_meth450k_f4)
info = batch.pcs %>% mutate(mnsi = clin$utmnsi[index],
                            age = clin$utalter[index],
                            sex = clin$lcsex[index]) %>% filter(!is.na(mnsi))
cell.type = read.csv("./Houseman/KF4_QN_BMIQ_estimated_cell_distribution_meanimpute473_lessThanOneFALSE.csv", sep = ";")
cell.type$ZZ_nr = as.character(cell.type$ZZ_nr)
info = left_join(info, cell.type)
info$Chip = as.factor(info$Chip)
info$Batch = as.factor(info$Batch)
info$sex = as.factor(info$sex)
beta = beta[, info$ZZ_nr]
detP = beta[, info$ZZ_nr]

####CpG site and sample filtering#####

beta.imp = champ.impute(beta = beta, pd = info)
beta.fil = champ.filter(beta = beta.imp$beta, pd = beta.imp$pd)

#####Technical effect correction######

champ.SVD(beta = beta.fil$beta, pd = beta.fil$pd)
#tech = dplyr::select(beta.fil$pd, -ZZ_nr)
#beta.cor = champ.runCombat(beta = beta.fil$beta, pd = beta.fil$pd, batchname = c("Plate","Chip.Position","sex"), variablename = "mnsi")

beta.fil.df = beta.fil$beta
meansd = data.frame(mean = apply(beta.fil.df,1,mean), sd = apply(beta.fil.df,1,sd), 
                    iqr = apply(beta.fil.df,1,IQR))
beta.topvar = beta.fil.df[rownames(meansd[meansd$iqr>quantile(meansd$iqr,0.9),]),]
pca = prcomp(t(beta.topvar), scale. = T)
pca_df = data.frame(pca$x, info)
pca_var = summary(pca)$importance
varp = round(pca_var[2,]*100,2)
ggplot(pca_df,aes(x=PC1,y=PC3,col = sex)) + 
  geom_point(size = 3, alpha = 0.5)+
  xlab(paste("PC1 (",varp[1],"%)")) +
  ylab(paste("PC2 (",varp[2],"%)")) +
  theme(panel.background = element_rect(fill = "white"))

modmatrix = model.matrix(~ 0 + mnsi, info)
covar = c(paste0("PC",seq(1:20)), "age", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
beta.nobatch = removeBatchEffect(x = log(beta.fil.df), 
                                 batch = info$Plate, batch2 = info$sex,
                                 covariates = info[,covar],
                                 design = modmatrix) %>%
  exp()

saveRDS(beta.nobatch, file = "beta_processed.rds")
