#!/bin/csh
#

#PBS  -N pangu_test
#PBS  -A xxxxxxxx
#PBS  -q main
#PBS  -l select=1:ncpus=4:mem=32G
#PBS  -l walltime=00:10:00
#PBS  -o pangu_test
#PBS  -j oe

setenv          LAUNCHCMD mpiexec

cd $path_to_your_working_directory/

${LAUNCHCMD} ./filter >& a.out
 
