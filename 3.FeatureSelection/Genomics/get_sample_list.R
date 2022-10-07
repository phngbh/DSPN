#!/usr/bin/env Rscript

sample_list = readRDS("/2.Partitioning/sample_partition.rds")
sample_list_gen = sample_list$featureSel$Genomics

for (i in 1:100){
	write.table(data.frame(V1 = Genomics[[i]], V2 = Genomics[[i]]), file = paste0("./sample_list/samples_",i,".txt"), colnames = F, rownames = F, sep = "\t")
}

sample_list_gen_all = unlist(sample_list$featureSel$Genomics) %>% unique()
write.table(data.frame(V1 = sample_list_gen_all, V2 = sample_list_gen_all), file = paste0("./sample_list/samples_all.txt"), colnames = F, rownames = F, sep = "\t")

train_samples = sample_list$train
clin <- read_sas("/data/gesamt_k05720g_v2007.sas7bdat", NULL)
clin = filter(clin, !utdm100y15 %in% c(6,7,8)) %>% droplevels()
sampleIDs = clin %>% dplyr::select(lg_dnaAxiom_s4f4, 
                           un_expr_f4ogtt,
                           un_protSoma_f4,
                           un_metabMetabolon_f4,
                           un_meth450k_f4,
                           zz_nr)
sampleIDs$Clinical = clin$zz_nr
colnames(sampleIDs) = c("Genomics", "Transcriptomics", "Proteomics",
                        "Metabolomics", "Methylomics", "Olink", "Clinical")
train_samples_gen = sampleIDs$Genomics[match(train_samples, as.character(sampleIDs$Clinical))]
write.table(data.frame(V1 = train_samples_gen, V2 = train_samples_gen), file = paste0("./sample_list/samples_train.txt"), colnames = F, rownames = F, sep = "\t")