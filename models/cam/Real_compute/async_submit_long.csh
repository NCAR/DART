#!/bin/tcsh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# Shell script to submit script which does repeated assim segments
# use 'qsub -l nodes=# async_submit_long.csh'  to make it run on # nodes

###  Job name
#PBS -N Exp9
###  Declare job non-rerunable
#PBS -r n
###  Output files
#PBS -e Exp9.err
#PBS -o Exp9.log
###  (small(20min), medium(2hr), long(12hr), verylong(72hr))
#PBS -q verylong
#
### This job's working directory; must cd to it, or it will run in /home...
cd $PBS_O_WORKDIR
###  Output to confirm job characteristics
echo Running $PBS_JOBNAME on host `hostname`
echo Time is `date`
echo Directory is `pwd`

echo This job runs on the following processors:
cat "$PBS_NODEFILE"

set node = `head -1 $PBS_NODEFILE` 
echo rsh $node "csh $PBS_O_WORKDIR/async_long_run.csh $PBS_O_WORKDIR $PBS_NODEFILE" 
rsh $node "csh $PBS_O_WORKDIR/async_long_run.csh $PBS_O_WORKDIR $PBS_NODEFILE" &

wait
