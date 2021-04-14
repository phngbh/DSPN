#!/bin/bash

#SBATCH --error=/storage/groups/cbm01/workspace/phong.nguyen_new/DSPN/sbatch_log/DSPN_%J.err
#SBATCH --output=/storage/groups/cbm01/workspace/phong.nguyen_new/DSPN/sbatch_log/DSPN_%J.out
#SBATCH --mem-per-cpu=60G
#SBATCH --cpus-per-task=1
#SBATCH --exclude=ibis216-010-[001-007]
#SBATCH -p cpu_p
#SBATCH -t 3-00:00:00
#SBATCH --nice=10000
#SBATCH --job-name=DSPN

echo Starting time is: $(date)

WDIR=/storage/groups/cbm01/workspace/phong.nguyen_new/DSPN
IMDIR=/storage/groups/cbm01/tools/alexander.ohnmacht/r-studio-charliecloud-master-b76b94e8c5040cd58fb0ffd1a463fd7409bb886e/exports
IMAGE=$1
TASK=$2
SCRIPT=$WDIR/training.R


echo Directory in lscratch: /localscratch/$USER/$IMAGE
echo Executing create-docker on $SLURM_JOB_NODELIST
mkdir -p /localscratch/$USER;
cd /localscratch/$USER;
#Check if directory exists
if [[ ! -d "/localscratch/$USER/$IMAGE" ]]
then
	echo "Container $IMAGE does not exists on your filesystem, making container"
	echo Created directory for local environment with docker image: $IMAGE	
	ch-tar2dir $IMDIR/$IMAGE.tar.gz /localscratch/$USER;
	echo Set up the docker image !
else
  echo Container $IMAGE already exists in localscratch
fi

sleep 120 


echo Running container and script 
ch-run -b /localscratch/:/localscratch -b /storage/groups/:/storage/groups /localscratch/$USER/$IMAGE/ -- Rscript $SCRIPT $TASK

# echo Deleting the image... 
# rm -rf /localscratch/$USER/$IMAGE/ 

echo Finishing time is $(date)