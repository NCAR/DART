#!/bin/csh

#
#BSUB -n 81
##BSUB -n 16
#BSUB -J convert
#BSUB -o output.driver
#BSUB -e output.driver
#BSUB -q regular
#BSUB -P P64000510
#BSUB -W 15
#BSUB -R "span[ptile=16]"

cd $PWD

rm -rf core*
rm -f obs_seq.00??
mpirun.lsf ./gsi_to_dart  >&! stdout
