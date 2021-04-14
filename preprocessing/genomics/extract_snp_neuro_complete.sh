#!/bin/bash

#SBATCH --error=/storage/groups/cbm01/workspace/phong.nguyen_new/DSPN/Genomics/sbatch_log/DSPN_%J.err
#SBATCH --output=/storage/groups/cbm01/workspace/phong.nguyen_new/DSPN/Genomics/sbatch_log/DSPN_%J.out
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=1
#SBATCH --exclude=ibis216-010-[001-007,011-012]
#SBATCH -p cpu_p
#SBATCH -t 3-00:00:00
#SBATCH --nice=10000
#SBATCH --job-name=DSPN

omics=$1


for i in {1..10}; do
	#Extract trainging samples
	./plink --bfile ./genomics --keep sample/${omics}_$i.txt --make-bed --out ./scratch/tmp_${omics}_$i

	#Extract LD-pruned SNPs
	./plink --bfile ./scratch/tmp_${omics}_$i --indep-pairwise 50 5 0.5 --out ./scratch/${omics}_$i
	./plink --bfile ./scratch/tmp_${omics}_$i --extract ./scratch/${omics}_$i.prune.in --make-bed --out LDpruned/${omics}_$i

	#GWAS for continuous phenotypes
	./plink --bfile LDpruned/${omics}_$i \
			--linear \
			--pheno neuro_con.txt \
			--all-pheno  \
			--allow-no-sex \
			--pfilter 5e-2 \
			--out ./res_neuro/gwas_${omics}_$i


	#GWAS for categorical phenotypes
	./plink --bfile LDpruned/${omics}_$i \
			--logistic \
			--pheno neuro_cat.txt \
			--all-pheno  \
			--1  \
			--allow-no-sex \
			--pfilter 5e-2 \
			--out ./res_neuro/gwas_${omics}_$i
done

#Extract trainging samples
./plink --bfile ./genomics --keep sample/${omics}.txt --make-bed --out ./scratch/tmp_${omics}

#Extract LD-pruned SNPs
./plink --bfile ./scratch/tmp_${omics} --indep-pairwise 50 5 0.5 --out ./scratch/${omics}
./plink --bfile ./scratch/tmp_${omics} --extract ./scratch/${omics}.prune.in --make-bed --out LDpruned/${omics}

#GWAS for continuous phenotypes
./plink --bfile LDpruned/${omics} \
		--linear \
		--pheno neuro_con.txt \
		--all-pheno  \
		--allow-no-sex \
		--pfilter 5e-2 \
		--out ./res_neuro/gwas_${omics}


#GWAS for categorical phenotypes
./plink --bfile LDpruned/${omics} \
		--logistic \
		--pheno neuro_cat.txt \
		--all-pheno  \
		--1  \
		--allow-no-sex \
		--pfilter 5e-2 \
		--out ./res_neuro/gwas_${omics}

#Extract top significant GWAS SNPs
while IFS= read -r line; do
	Rscript extract_sig_snps_complete.R ${omics} ${line}
done < neuro.txt


#Make small bed files 
for file in ./sig_snps_neuro/${omics}_*.txt; do
    
	name=${file#"./sig_snps_neuro/"}
	name=${name%".txt"}
	./plink --bfile genomics --extract $file --make-bed --out ./small_bed_neuro/${name}

done 
