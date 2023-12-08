#!/bin/bash

#PBS -A P86850054
#PBS -N build-everything
#PBS -j oe
#PBS -k eod
#PBS -q main
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=128:mpiprocs=128

export TMPDIR=/glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

time ./run_all_quickbuilds.sh
