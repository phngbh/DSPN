#!/bin/bash

./plink --bfile /1.Processing/Genomics/genomics_processed --extract edge_snps.txt --keep ./samples_list/samples_train.txt --make-bed --out genomics_selected