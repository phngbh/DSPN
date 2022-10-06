#!/usr/bin/env Rscript

library(dplyr)

gene_snps = read.table(as.character(args[1]), header = F, sep = "", check.names = F)
snps = dplyr::select(gene_snps, -V1, -V2) %>% t()
gene_annot_df = data.frame(Gene = rep(gene_snps$V1, length(snps[,1])), SNP = snps[,1])

write.table(gene_annot_df, file = paste0("gene_annot.txt"), col.names = F, row.names = F, quote = F, append = T)