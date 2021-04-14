#!/usr/bin/env Rscript

library(dplyr, lib.loc = "./Rlib")

args = commandArgs(trailingOnly=TRUE)

for (i in 1:10){
	if (as.character(args[2]) %in% c("un14","un15","un17","un18","utmnsi")){
		plink = read.table(paste0("./res_neuro/gwas_",args[1],"_",i,".",args[2],".assoc.linear"), header = T, sep = "", check.names = F)
	} else {
		plink = read.table(paste0("./res_neuro/gwas_",args[1],"_",i,".",args[2],".assoc.logistic"), header = T, sep = "", check.names = F)
	}

	plink = plink[order(plink$P),] %>% head(n=1000)
	write.table(data.frame(SNP = plink$SNP), file = paste0("./sig_snps_neuro/",args[1],"_",args[2],"_",i,".txt"), row.names = F, col.names = F, quote = F)
}

if (as.character(args[2]) %in% c("un14","un15","un17","un18","utmnsi")){
	plink = read.table(paste0("./res_neuro/gwas_",args[1],".",args[2],".assoc.linear"), header = T, sep = "", check.names = F)
} else {
	plink = read.table(paste0("./res_neuro/gwas_",args[1],".",args[2],".assoc.logistic"), header = T, sep = "", check.names = F)
}

plink = plink[order(plink$P),] %>% head(n=1000)
write.table(data.frame(SNP = plink$SNP), file = paste0("./sig_snps_neuro/",args[1],"_",args[2],".txt"), row.names = F, col.names = F, quote = F)




