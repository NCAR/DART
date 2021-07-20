#!/bin/bash -l

#PBS -A P86850054
#PBS -N debug
#PBS -j oe
#PBS -k eod
#PBS -q regular
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mpiprocs=1

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

module load arm-forge
export MPI_SHEPHERD=true
ddt --connect ./filter
