#!/bin/bash

# Submits num_cyles of a tiegcm ensemble run. 
#   assimilation each cycle
#
# run5.0.pbs is an array job
#  To check exit code for each arry job 
#  qhist -j JOB_ID -l


num_cycles=25
MODEL_RUNS=$(qsub run5.0.pbs)

echo "Submitted " $MODEL_RUNS

for (( i=1; i<$num_cycles; i++))
do
  FILTER_RUN=$(qsub -W depend=afterok:$MODEL_RUNS submit_filter.pbs)
  MODEL_RUNS=$(qsub -W depend=afterok:$FILTER_RUN run5.0.pbs)
  echo "Submitted cycle" $i
done
