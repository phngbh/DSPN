#!/usr/bin/env Rscript

sample_list = readRDS("/2.Partitioning/sample_partition.rds")
sample_list_gen = sample_list$featureSel$Genomics

for (i in 1:100){
	write.table(data.frame(V1 = Genomics[[i]], V2 = Genomics[[i]]), file = paste0("./sample_list/samples_",i,".txt"), colnames = F, rownames = F, sep = "\t")
}

sample_list_gen_all = unlist(sample_list$featureSel$Genomics) %>% unique()
write.table(data.frame(V1 = sample_list_gen_all, V2 = sample_list_gen_all), file = paste0("./sample_list/samples_all.txt"), colnames = F, rownames = F, sep = "\t")