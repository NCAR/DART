#!/bin/ksh -x
#
# Script to run matlab
module load idl
bsub -Is -q caldera -W 1:00 -n 1 -P P19010000 idl
exit
