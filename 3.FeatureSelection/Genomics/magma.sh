#!/bin/bash

dirlist=(subsets annotates magma_gene magma_geneset)
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

echo Annotating SNPs to genes...  | sed G

./magma --annotate window=2,0.5 --snp-loc ./genomics.bim --gene-loc ./NCBI37.3/NCBI37.3.gene.loc --out ./annotates/annotates

#### Cross-sectional DSPN #######

echo Starting process for prevalent DSPN phenotype | sed G

for i in {1..100};
do
	###Removing training samples
	./plink --silent --bfile ./genomics --keep ./samples_mnsi3/samples_${i}.txt --make-bed --out ./subsets/mnsi3_${i}

	echo Running gene-level analysis... | sed G
  # Gene-level analysis for 1 bootstrapped subset
	./magma --bfile ./subsets/mnsi3_${i} --gene-annot ./annotates/annotates.genes.annot --pheno file=mnsi3.txt --out ./magma_gene/mnsi3_${i}
  
	echo Running gene set-level analysis... | sed G
  # Gene-set level analysis for 1 bootstrapped subset
	./magma --gene-results ./magma_gene/mnsi3_${i}.genes.raw --model alpha=0.3 --settings gene-info --set-annot geneset.txt --out ./magma_geneset/mnsi3_${i}

done

# Extract entire feature selection set
./plink --silent --bfile ./genomics --keep samples_mnsi3.txt --make-bed --out fs_mnsi3

# Gene analysis for entire feature selection set
./magma --bfile ./fs_mnsi3 --gene-annot ./annotates/annotates.genes.annot --pheno file=mnsi3.txt --out ./magma_gene/mnsi3

# Final GSEA to extract the leading edge SNPs
Rscript final_gsea.R mnsi3

# Extract final bed file 
./plink --silent -bfile ./genomics --extract leadingEdge_snps_mnsi3.txt --make-bed --out genomics_mnsi3

#### Incident DSPN #######

echo Starting process for incident DSPN phenotype | sed G

for i in {1..100};
do
	###Removing training samples
	./plink --silent --bfile ./genomics --keep ./samples_inc3/samples_${i}.txt --make-bed --out ./subsets/inc3_${i}

	echo Running gene-level analysis... | sed G
  # Gene-level analysis for 1 bootstrapped subset
	./magma --bfile ./subsets/inc3_${i} --gene-annot ./annotates/annotates.genes.annot --pheno file=inc3.txt --out ./magma_gene/inc3_${i}
  
	echo Running gene set-level analysis... | sed G
  # Gene-set level analysis for 1 bootstrapped subset
	./magma --gene-results ./magma_gene/inc3_${i}.genes.raw --model alpha=0.3 --settings gene-info --set-annot geneset.txt --out ./magma_geneset/inc3_${i}

done

# Extract entire feature selection set
./plink --silent --bfile ./genomics --keep samples_inc3.txt --make-bed --out fs_inc3

# Gene analysis for entire feature selection set
./magma --bfile ./fs_inc3 --gene-annot ./annotates/annotates.genes.annot --pheno file=inc3.txt --out ./magma_gene/inc3

# Final GSEA to extract the leading edge SNPs
Rscript final_gsea.R inc3

# Extract final bed file 
./plink --silent -bfile ./genomics --extract leadingEdge_snps_inc3.txt --make-bed --out genomics_inc3