#!/bin/bash

#PBS -A P86850054
#PBS -N compile
#PBS -j oe
#PBS -k eod
#PBS -q regular
#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=36

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

time ./run_all_quickbuilds.sh
