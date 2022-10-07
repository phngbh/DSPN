#!/usr/bin/env Rscript

setwd("/storage/groups/cbm01/workspace/phong.nguyen/DSPN/")
libs <- c("/storage/groups/cbm01/workspace/phong.nguyen/DSPN/Rlib", .libPaths())
.libPaths(libs)

args = commandArgs(trailingOnly=TRUE)

source("./ml3.R")

###Load data####
cat("Load data\n")
gen = bed_to_df_copyNumber("./genomics_inc3") %>% standardise()
tra = readRDS("./tra_inc3.rds") %>% t() %>% standardise()
tra_thr = readRDS("./tra_inc3_thr.rds") %>% t() %>% standardise()
pro = readRDS("./pro_inc3.rds") %>% t() %>% standardise()
pro_thr = readRDS("./pro_inc3_thr.rds") %>% t() %>% standardise()
meta = readRDS("./meta_inc3.rds") %>% t() %>% standardise()
meta_thr = readRDS("./meta_inc3_thr.rds") %>% t() %>% standardise()
meth = readRDS("./meth_inc3.rds") %>% t() %>% standardise()
oli = readRDS("./olink_processed.rds") %>% as.matrix() %>% standardise()
clin = readRDS("./cli_inc3.rds")  %>% as.matrix() %>% standardise()
data  = list(Genomics = gen, Transcriptomics = tra, Proteomics = pro, Metabolomics = meta, Methylomics = meth, Olink = oli, Clinical = clin)
rm(gen)
rm(meth)
rm(tra)
omics_id = readRDS("./omics_id.rds") 
mnsi = readRDS("./mnsi.rds")
inc = dplyr::select(mnsi, inc3) %>% na.omit

#Make overlapped datasets and CV partitions
data_list_inc3 = make_overlap_set(data, omics_id, inc)
data_list_inc3 = data_list_inc3[c("Transcriptomics","Proteomics","Metabolomics","Clinical")]
data_list_inc3_thr = list(Transcriptomics = tra_thr[as.character(omics_id$Transcriptomics)[match(rownames(data_list_inc3$Clinical), as.character(omics_id$Clinical))],], 
                          Proteomics = pro_thr[as.character(omics_id$Proteomics)[match(rownames(data_list_inc3$Clinical), as.character(omics_id$Clinical))],], 
                          Metabolomics = meta_thr[as.character(omics_id$Metabolomics)[match(rownames(data_list_inc3$Clinical), as.character(omics_id$Clinical))],], 
                          Clinical = clin[rownames(data_list_inc3$Clinical),])
rownames(data_list_inc3_thr$Transcriptomics) = rownames(data_list_inc3$Clinical)
rownames(data_list_inc3_thr$Proteomics) = rownames(data_list_inc3$Clinical)
rownames(data_list_inc3_thr$Metabolomics) = rownames(data_list_inc3$Clinical)
inc = inc[rownames(data_list_inc3$Clinical),]
names(inc) = rownames(data_list_inc3$Clinical)
inc = ifelse(inc == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
folds_inc3 = make_cv_list_ens(outcome = inc, times = 100)

# #Fit concatenated model
# fit_cat = fit_concatenation(data_list = data_list_inc3, y = inc, cv_list = folds_inc3, n_cv =  5, p_metric = "wLogLoss", n = as.numeric(args[[1]]))
# saveRDS(fit_ens, paste0("benchmarking/cat_iter",args[[1]],".rds"))

#Fit ensemble model
fit_ens = fit_ensemble(data_list = data_list_inc3, y = inc, cv_list = folds_inc3, n_cv =  5, p_metric = "wLogLoss", n = as.numeric(args[[1]]))
saveRDS(fit_ens, paste0("benchmarking/ensemble_iter",args[[1]],".rds"))

fit_ens = fit_ensemble(data_list = data_list_inc3_thr, y = inc, cv_list = folds_inc3, n_cv =  5, p_metric = "wLogLoss", n = as.numeric(args[[1]]))
saveRDS(fit_ens, paste0("benchmarking/ensemble_thr_iter",args[[1]],".rds"))

#Forward feature selection
fit_fwSel = fit_forwardSelect(data_list = data_list_inc3, y = inc, cv_list = folds_inc3, p_metric = "wLogLoss", n_datasets = 4, n_cv = 5, n = as.numeric(args[[1]]))
saveRDS(fit_fwSel, paste0("benchmarking/fwSel_iter",args[[1]],".rds"))

fit_fwSel = fit_forwardSelect(data_list = data_list_inc3_thr, y = inc, cv_list = folds_inc3, p_metric = "wLogLoss", n_datasets = 4, n_cv = 5, n = as.numeric(args[[1]]))
saveRDS(fit_fwSel, paste0("benchmarking/fwSel_thr_iter",args[[1]],".rds"))