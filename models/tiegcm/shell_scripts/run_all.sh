#!/bin/bash

# To check exit code for each arry job 
# qhist -j 2667681 -l


MODEL_RUNS=$(qsub run2.5.pbs)

echo "Submitted " $MODEL_RUNS

qsub -W depend=afterok:$MODEL_RUNS submit_filter.pbs
