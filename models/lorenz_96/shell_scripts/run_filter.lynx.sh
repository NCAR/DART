#!/bin/sh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# The general example run script is 'run_filter.csh' in this same
# directory.  This is a specialized version of the script that runs
# under MOAB/Torque on the NCAR Cray test system.
#
#
# first, it appears this submit script must be bash or sh; csh doesn't work
#
# to submit this job:  qsub < run_filter.lynx
# to check the queues: qstat
#
##=============================================================================
# an explanation of the directives follows:
# -N        job name
# -e <arg>  filename for standard error
# -o <arg>  filename for standard out 
# -q <arg>  queue name
# -l mppwidth=xx:walltime=hh:mm:ss    
#           xx is number of mpi tasks, walltime is max run time.
##=============================================================================
#PBS -N filter
#PBS -e filter.err
#PBS -o filter.log
#PBS -q regular
#PBS -l mppwidth=4,walltime=00:30:00

echo 'Starting execution of the run_filter script for the Lorenz 96 model.'
echo " "

# this script appears to start executing in your home dir.
echo starting directory is `pwd`
echo " "

# cd down to where the executables were built.
cd DART/models/lorenz_96/work
echo current working directory now `pwd`
echo " "

# the home directories are NOT mounted on the execution nodes,
# so everything you need must be copied over to /ptmp BEFORE
# starting filter.   (they are visible from the head node that is
# running this script, but once you start aprun, you must have
# your executables and data files under /ptmp.)
cp -fv filter input.nml obs_seq.out filter_ics  /ptmp/${PBS_O_LOGNAME}/work
echo " "

cd /ptmp/${PBS_O_LOGNAME}/work
echo current working directory now `pwd`
echo " "

echo ensuring the 4 files required to run this job are here:
ls -l filter input.nml obs_seq.out filter_ics
echo " "

# and now, run with 4 MPI tasks
aprun -n 4 ./filter

# the log files will be written to whatever dir this script 
# was originally submitted from.

echo " "
echo 'Ending execution of the run_filter script for the Lorenz 96 model.'
echo " "

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

