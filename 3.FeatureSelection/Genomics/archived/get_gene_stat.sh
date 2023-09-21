#!/bin/bash

awk '!/#/ {print}' ./annotates/annotate.genes.annot > ./scratch/annotates_tmp.txt

echo Annotating SNPs to genes
while IFS= read -r line; do
	echo $line > gene_tmp.txt
	Rscript make_gene_annot.R gene_tmp.txt $phenotype
done < ./scratch/annotates_tmp.txt

# echo Getting gene mean statistics
# Rscript get_gene_meanStat.R gene_annot_${phenotype}.txt ./gwas/gwas.${phenotype}.assoc.logistic $phenotype
