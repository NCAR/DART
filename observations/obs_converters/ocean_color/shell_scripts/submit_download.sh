#!/bin/bash

#PBS -N nasa_chlora_download
#PBS -A p93300012
#PBS -l select=1:ncpus=4:mpiprocs=1
#PBS -l walltime=05:00:00
#PBS -q casper
#PBS -m abe

./get_ocdata.sh > download_stderr.txt 2> download_stdout.txt