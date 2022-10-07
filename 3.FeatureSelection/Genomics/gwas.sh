#!/bin/bash

dirlist=(scratch samples_list gwas)
for dir in  ${dirlist[@]} 
do
    if [[ -d ${dir} ]] 
    then
        echo "Folder ${dir} exists, making new folder..." 
        rm -r ${dir}
        mkdir ${dir}
    else
        echo Making ${dir}
        mkdir ${dir} 
    fi
done 

#Get sample list
Rscript get_sample_list.R

#Make phenotype file
Rscript make_phenotype_file.R

#GWAS
echo Doing GWAS and make BED files of the results
for i in {1..100};
do 

	./plink --bfile /1.Processing/Genomics/genomics_processed --extract snps_2_keep.txt --keep ./samples_list/samples_${i}.txt --make-bed --out ./scratch/tmp
	./plink --bfile ./scratch/tmp \
					--logistic \
					--pheno inc3.txt \
					--all-pheno \
					--allow-no-sex \
					--out ./gwas/gwas_sample_$i
	#Rscript make_bed.R ./gwas/gwas.inc3.assoc.logistic inc3
done

./plink --bfile /1.Processing/Genomics/genomics_processed --extract snps_2_keep.txt --keep ./samples_list/samples_all.txt --make-bed --out ./scratch/tmp
./plink --bfile ./scratch/tmp \
					--logistic \
					--pheno inc3.txt \
					--all-pheno \
					--allow-no-sex \
					--out ./gwas/gwas_all
