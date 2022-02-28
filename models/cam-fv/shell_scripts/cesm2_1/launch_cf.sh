#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$

# for command file jobs.
# Sidd Ghosh Feb 22, 2017
# Slurm added by Kevin Raeder July 6, 2019

# On casper's PBS the PMI_RANK variable is not defined,
# but  OMPI_COMM_WORLD_LOCAL_RANK is (different 
# from OMPI_COMM_WORLD_RANK, which is also defined)
# Also, PBS_O_WORKDIR is defined, so it goes into the wrong section.
# The system launch_cf.sh (cheyenne only) tests directly for -z {env_var_name}
# export 

if [ ! -z "$PMI_RANK" ]; then
   line=$(expr $PMI_RANK + 1)
#    echo "launch_cf.sh using PMI_RANK with line = $line"
elif [ ! -z "$OMPI_COMM_WORLD_RANK" ]; then
   line=$(expr $OMPI_COMM_WORLD_RANK + 1)
#    echo "launch_cf.sh using OMPI_COMM_WORLD_RANK with line = $line"
else
   echo "Batch environment is unknown"
   exit 11
fi

INSTANCE=$(sed -n ${line}p $1)

# The following command showed that 563 tasks are launched within .3 seconds.
# echo "launching $INSTANCE at "; date --rfc-3339=ns

eval "$INSTANCE"

# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$
