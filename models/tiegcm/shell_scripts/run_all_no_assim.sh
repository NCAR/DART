#!/bin/bash

# Submits num_cyles of a tiegcm ensemble run. 
# Tiegcm stop start only, no assimilation
#
# run5.0.pbs is an array job
#   To check exit code for each arry job 
#   qhist -j JOB_ID -l

num_cycles=25
MODEL_RUNS=$(qsub run5.0.pbs)

echo "Submitted " $MODEL_RUNS

for (( i=1; i<$num_cycles; i++))
do
  MODEL_RUNS=$(qsub -W depend=afterok:$MODEL_RUNS run5.0.pbs)
  echo "Submitted " $MODEL_RUNS
done


