#!/bin/bash

# Where are the scripts located
WDIR="/lustre/groups/cbm01/workspace/phong.nguyen/DSPN/"
JOBSCRIPT=$1
DOCKER=$2
EXECSCRIPT=$3

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

for i in {1..100};
do
	echo Executing task $i
	sbatch $WDIR/$JOBSCRIPT $DOCKER $i $EXECSCRIPT
	sleep 5
done

echo Submitted all jobs !