#!/bin/bash

# Where are the scripts located
WDIR="/storage/groups/cbm01/workspace/phong.nguyen_new/DSPN/"
SCRIPT=$1
DOCKER=$2

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

for i in 38 50 56;
do
	echo Executing omics combination $i
	sbatch $WDIR/$SCRIPT $DOCKER $i
	sleep 5
done

echo Submitted all jobs !