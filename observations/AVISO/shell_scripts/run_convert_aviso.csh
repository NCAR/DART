#!/bin/csh -f
#
## DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Convert a series of sea surface topography files to DART observation sequence
# files. The list of satellites and years is defined IN the convert_aviso_2.py
# script. The directory names specifying the input and output directories is
# also specified in the convert_aviso_2.py script.
#
# This script can be submitted as a batch job or run interactively.
#
#BSUB -n 1
#BSUB -o out.%J
#BSUB -e err.%J
#BSUB -W 23:59
#BSUB -q geyser
#BSUB -P P86850054
#BSUB -J convert_aviso
#
#PBS -N convert_aviso      
#PBS -o convert_aviso.out
#PBS -e convert_aviso.err
#PBS -l walltime=01:00:00
#PBS -q economy 
#PBS -l select=1:ncpus=1
#PBS -A P86850054 

# account for the fact that some PBS/SLURM implementations need
# to be told where to run.
if (${?PBS_O_WORKDIR}) then
  cd ${PBS_O_WORKDIR}
endif

if ( -e input.nml ) then
   # just use existing input.nml
else
   # Set the default value for the observation error standard deviation.
   # If convert_aviso_2.py is modified to call convert_aviso with a second
   # arguments specifying specific observation errors, this default value
   # is ignored. 

   set OBS_ERR_STD=0.03

   sed -e "s/<OBS_ERR_STD>/${OBS_ERR_STD}/g" < input.nml.template >! input.nml
endif

python convert_aviso_2.py || exit 1

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
