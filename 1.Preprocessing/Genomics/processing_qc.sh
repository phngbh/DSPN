#!/bin/bash

#USAGE: processing_qc.sh <working directory> <data directory> 

WDIR=$1 #/mnt/nas/global/ame-ige/DSPN_analysis/genotype
DDIR=$2 #/mnt/nas/global/ame-ige/omicsdata/genotypes/KORAS4F4_affyaxiom_n3788/imputed/HRC/HRCPanel_r1.1/chr1_22/vcf

#Set working directory
echo Working directory is ${WDIR} | sed G
cd $WDIR

#Create directories 
echo Creating working folders...

dirlist=(scratch qc vcf binary)
for dir in  ${dirlist[@]} 
do
	if [[ -d ${WDIR}/${dir} ]] 
	then
		echo "Folder ${dir} exists, clean folder..." | sed G
		rm -rf ${dir}
		mkdir ${dir}
	else
		echo Making ${dir} | sed G
		mkdir ${dir} 
	fi
done 


echo Checking input files...

#Check if input files exist 
filename=($(ls $DDIR))
if [[ -z "$filename" ]] 
then
	echo No genotyping file found
	echo Exiting | sed G
	exit 1
 else
   echo Found input files | sed G
fi 


echo Finish checking | sed G

echo Starting input transformation... | sed G | sed G

##TRANSFORMING INPUT AND MAKING PLINK BINARY FILE######
for chr in {1..22}
do

	echo "Extracting chromosome ${chr}" | sed G 
	gunzip -c $DDIR/KORAS4F4_HRC_N3788_chr${chr}_notmonomorph.vcf.gz > ./vcf/imputed_chr$chr.vcf 

    echo "Making plink binary fileset for chromosome ${chr}" | sed G
	#Make plink binary files
	./plink --vcf ./vcf/imputed_chr${chr}.vcf \
		  --silent \
		  --no-parents \
		  --no-sex \
		  --no-pheno \
		  --make-bed --out ./binary/chr"$chr"
  
  echo "Done proccessing chromosome ${chr}" | sed G | sed G
  
done 

###QUALITY CONTROL##########

echo Merging all chromosomes | sed G
./plink --silent --bfile ./binary/chr1 --merge-list bifilename.txt --make-bed --out ./binary/allchr

echo Start quality control | sed G

#Remove SNPs with missing rate >10%, MAF<0.05 and high deviation from HWE (p<e-10)
echo "Removing SNPs that have typing rate < 99%, MAF < 1% and significant deviation from HWE (p<1e-10)..." | sed G
./plink --bfile ./binary/allchr \
	    --silent \
	    --snps-only \
	    --geno 0.01 \
	    --maf 0.01 \
	    --hwe 1e-10 \
	    --make-bed --out ./qc/allchr_freqfiltered

echo Removing samples that have heterozygosity rate deviating 3sd from mean...
./plink --bfile ./qc/allchr_freqfiltered \
		--silent \
        --indep-pairwise 50 5 0.2 \
        --out ./qc/indepSNP_allchr
./plink --silent --bfile ./qc/allchr_freqfiltered --extract ./qc/indepSNP_allchr.prune.in --het --out ./qc/het_allchr
Rscript --no-save heterozygosity_outliers_list.R 
sed 's/"// g' ./qc/fail-het-qc-allchr.txt | awk '{print$1, $2}'> ./scratch/het_fail_ind.txt
./plink --silent --bfile ./qc/allchr_freqfiltered --remove ./scratch/het_fail_ind.txt --make-bed --out ./qc/allchr_QCed


echo "Finish QC check for chromosome ${chr}. Final QCed files at ./qc/allchr_QCed" | sed G | sed G

###FILTERING#####
echo Filter for samples that have clinical records | sed G

./plink --bfile ./qc/allchr_QCed --keep sampleIDs.txt --make-bed --out ./genomics_processed


