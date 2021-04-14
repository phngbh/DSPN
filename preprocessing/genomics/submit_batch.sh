#!/bin/bash

# Where are the scripts located
WDIR="/storage/groups/cbm01/workspace/phong.nguyen_new/DSPN/Genomics"
SCRIPT=$1

dirlist=(scratch LDpruned res_neuro sig_snps_neuro small_bed_neuro)
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

echo Starting arrayscript job ...
echo USER: $USER
# some initial info
echo It is now:
date
echo
echo Running on machine
hostname
echo
echo Operating system
uname -r

##############################
### The actual task submission
##############################

# code to execute with passing $TASKID to the R script

while IFS= read -r omics; do
	echo Executing omics combination ${omics}
	sbatch $WDIR/$SCRIPT ${omics}
	sleep 5
done < omics.txt

echo Submitted all jobs !