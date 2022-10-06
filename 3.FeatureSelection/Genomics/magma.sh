#!/bin/bash

dirlist=(subsets annotates)
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

./magma --annotate window=2,0.5 --snp-loc /1.Processing/Genomics/genomics_processed.bim --gene-loc ./NCBI37.3/NCBI37.3.gene.loc --out ./annotates/annotate 