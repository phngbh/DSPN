#!/usr/bin/env Rscript

setwd("/storage/groups/cbm01/workspace/phong.nguyen/DSPN/")
libs <- c("/storage/groups/cbm01/workspace/phong.nguyen/DSPN/Rlib", .libPaths())
.libPaths(libs)

args = commandArgs(trailingOnly=TRUE)

source("./ml3.R")

cat("Cross-sectional DSPN\n\n")

###Load data####
cat("Load data\n")
gen = bed_to_df_copyNumber("./genomics_mnsi3") %>% standardise()
tra = readRDS("./tra_mnsi3.rds") %>% t() %>% standardise()
pro = readRDS("./pro_mnsi3.rds") %>% t() %>% standardise()
meta = readRDS("./meta_mnsi3.rds") %>% t() %>% standardise()
meth = readRDS("./meth_mnsi3.rds") %>% t() %>% standardise()
oli = readRDS("./olink_processed.rds") %>% as.matrix() %>% standardise()
clin = readRDS("./cli_mnsi3.rds")  %>% as.matrix() %>% standardise()
data  = list(Genomics = gen, Transcriptomics = tra, Proteomics = pro, Metabolomics = meta, Methylomics = meth, Olink = oli, Clinical = clin)
rm(gen)
rm(meth)
rm(tra)
omics_id = readRDS("./omics_id.rds") 
mnsi = readRDS("./mnsi.rds")
crs = dplyr::select(mnsi, utmnsi3) %>% na.omit
inc = dplyr::select(mnsi, inc3) %>% na.omit

#Make overlapped datasets and CV partitions
data_list_mnsi3 = make_overlap_set(data, omics_id, crs)
crs = crs[rownames(data_list_mnsi3$Clinical),] 
names(crs) = rownames(data_list_mnsi3$Clinical)
crs = ifelse(crs == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
folds_mnsi3 = make_cv_list_ens(outcome = crs, times = 100)

# #Fit ensemble models
# fit_mnsi3 = fit_ensemble(data_list = data_list_mnsi3, y = crs, cv_list = folds_mnsi3, n = as.numeric(args[[1]]))
# saveRDS(fit_mnsi3,paste0("ensemble_elnet__elnet_w_roc/mnsi3_fold",args[[1]],".rds"))

#Forward feature selection
perf_mnsi3 = fit_forwardSelectFromClinical(data_list = data_list_mnsi3, y = crs, cv_list = folds_mnsi3, p_metric = "wLogLoss", n_datasets = 7, n_cv = 5, n = as.numeric(args[[1]]))
saveRDS(perf_mnsi3, paste0("~/fwSelection_ll_cli_cv5_t80_new_100_3/mnsi3_fold",args[[1]],".rds"))

perf_mnsi3 = fit_forwardSelect(data_list = data_list_mnsi3, y = crs, cv_list = folds_mnsi3, p_metric = "wLogLoss", n_datasets = 7, n_cv = 5, n = as.numeric(args[[1]]))
saveRDS(perf_mnsi3, paste0("~/fwSelection_ll_all_cv5_t80_new_100_3/mnsi3_fold",args[[1]],".rds"))

cat("Incident DSPN\n\n")

###Load data####
cat("Load data\n")
gen = bed_to_df_ens("./genomics_inc3") %>% standardise()
tra = readRDS("./tra_inc3.rds") %>% t() %>% standardise()
pro = readRDS("./pro_inc3.rds") %>% t() %>% standardise()
meta = readRDS("./meta_inc3.rds") %>% t() %>% standardise()
meth = readRDS("./meth_inc3.rds") %>% t() %>% standardise()
oli = readRDS("./olink_processed.rds") %>% as.matrix() %>% standardise()
clin = readRDS("./cli_inc3.rds")  %>% as.matrix() %>% standardise()
data  = list(Genomics = gen, Transcriptomics = tra, Proteomics = pro, Metabolomics = meta, Methylomics = meth, Olink = oli, Clinical = clin)
rm(gen)
rm(meth)
rm(tra)

#Make overlapped datasets and CV partitions
data_list_inc3 = make_overlap_set(data, omics_id, inc)
inc = inc[rownames(data_list_inc3$Clinical),]
names(inc) = rownames(data_list_inc3$Clinical)
inc = ifelse(inc == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
folds_inc3 = make_cv_list_ens(outcome = inc, times = 100)

# #Fit ensemble models
# fit_inc3 = fit_ensemble(data_list = data_list_inc3, y = inc, cv_list = folds_inc3, n = as.numeric(args[[1]]))
# saveRDS(fit_inc3,paste0("ensemble_elnet__elnet_w_roc/inc3_fold",args[[1]],".rds"))

#Forward feature selection
perf_inc3 = fit_forwardSelectFromClinical(data_list = data_list_inc3, y = inc, cv_list = folds_inc3, p_metric = "wLogLoss", n_datasets = 7, n_cv = 5, n = as.numeric(args[[1]]))
saveRDS(perf_inc3, paste0("~/fwSelection_ll_cli_cv5_t80_new_100_3/inc3_fold",args[[1]],".rds"))

perf_inc3 = fit_forwardSelect(data_list = data_list_inc3, y = inc, cv_list = folds_inc3, p_metric = "wLogLoss", n_datasets = 7, n_cv = 5, n = as.numeric(args[[1]]))
saveRDS(perf_inc3, paste0("~/fwSelection_ll_all_cv5_t80_new_100_3/inc3_fold",args[[1]],".rds"))