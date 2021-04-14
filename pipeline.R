#!/usr/bin/env Rscript

setwd("/storage/groups/cbm01/workspace/phong.nguyen_new/DSPN/")
libs <- c("/storage/groups/cbm01/workspace/phong.nguyen_new/DSPN/Rlib", .libPaths())
.libPaths(libs)

args = commandArgs(trailingOnly=TRUE)

source("./training/algorithm.R")

###Load data####
cat("Load data\n")
#gen = readRDS("./genomics_processed.rds") %>% t()
tra = readRDS("./expression_processed.rds") %>% t()
pro = readRDS("./proteins_processed.rds") %>% t()
meta = readRDS("./metabolomics_processed.rds") %>% t()
meth = readRDS("./beta_processed.rds")
oli = readRDS("./olink_processed.rds")
omics = list(Transcriptomics = tra, Proteomics = pro,
             Metabolomics = meta, Methylomics = meth, Olink = oli)
rm(meth)
rm(tra)

#omics = readRDS("./omics.rds")
lab = readRDS("./phenotype.rds") %>% column_to_rownames("zz_nr") %>% dplyr::select(-utmnsi)
keep_pt = apply(lab,2,function(x) sum(is.na(x)) < 200)
lab = lab[,keep_pt, drop=F]
covar = readRDS("./covar.rds")
omics_id = readRDS("./omics_id.rds")
neuro = readRDS("./neuro.rds")

###Manipulate data for training####
cat("Manipulating data for training\n")
data = make_comb_list(omics_list = c("Genomics","Transcriptomics","Proteomics",
                                     "Metabolomics","Methylomics","Olink"),
                      omics = omics,
                      covar = covar,
                      lab = lab,
                      omics_id = omics_id)

#Make reduced datasets
# for (i in 1:(length(data) -1)){
#   for(o in names(data[[i]][["Samples"]])){
#     data[[i]][["Samples"]][[o]] <- data$GenomicsTranscriptomicsProteomicsMetabolomicsMethylomicsOlink$Samples[[o]]
#   }
# }

data_split = make_train_test(omics_list = data)
train = extract_train_test(data, data_split, train = T)
test = extract_train_test(data, data_split, train = F)
folds = make_cv_list(omics_list = train)


### Fit and cross-validate with random forest
cat("Using elastic net\n")
make_fit(train_list = train, cv_list = folds, neuropathy = neuro, covar = covar, lab = lab, cv=T, n = as.numeric(args[1]), method = "elnet", omics = omics)
make_fit(train_list = train, test_list = test, cv_list = folds, neuropathy = neuro, lab = lab, covar = covar, cv=F, n = as.numeric(args[1]), method = "elnet", omics = omics)

### Fit and cross-validate with random forest
cat("Using random forest\n")
make_fit(train_list = train, cv_list = folds, neuropathy = neuro, covar = covar, lab = lab, cv=T, n = as.numeric(args[1]), method = "rf", omics = omics)
make_fit(train_list = train, test_list = test, cv_list = folds, neuropathy = neuro, lab = lab, covar = covar, cv=F, n = as.numeric(args[1]), method = "rf", omics = omics)

cat("Finished\n")
