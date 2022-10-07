#!/bin/bash

#SBATCH --error=/lustre/groups/cbm01/workspace/phong.nguyen/DSPN/sbatch_log/DSPN_%J.err
#SBATCH --output=/lustre/groups/cbm01/workspace/phong.nguyen/DSPN/sbatch_log/DSPN_%J.out
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=1
#SBATCH --nodelist=cpusrv[21-28]
#SBATCH -p cpu_p
#SBATCH -t 3-00:00:00
#SBATCH --nice=10000
#SBATCH --job-name=DSPN

echo Starting time is: $(date)

WDIR=/lustre/groups/cbm01/workspace/phong.nguyen/DSPN
IMDIR=/lustre/groups/cbm01/tools/work_on_servers/docker/images
IMAGE=$1
TASK=$2
SCRIPT=$3


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
ch-run -b /localscratch/:/localscratch -b /lustre/groups/:/storage/groups -b /home/icb/phong.nguyen:/home/phong.nguyen /localscratch/$USER/$IMAGE/ -- Rscript /storage/groups/cbm01/workspace/phong.nguyen/DSPN/$SCRIPT $TASK

# echo Deleting the image... 
# rm -rf /localscratch/$USER/$IMAGE/ 

echo Finishing time is $(date)