#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# test_dart.csh can be run from the command line or a batch system.
#               This compiles many of the programs (but not all) and
#               runs a limited number of tests.
#
#==========================================================================
# SLURM directives              sbatch test_dart.csh
#
# sinfo     information about the whole slurm system
# squeue    information about running jobs
# sbatch    submitting a job
# scancel   killing a job
#
#SBATCH --ignore-pbs
#SBATCH --job-name dart_test
#SBATCH -t 2:00:00
#SBATCH -A P86850054
#SBATCH -p dav
#SBATCH -o dart_test.log
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#
# for mpi tests:
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
# for serial tests:
#SB#### --ntasks=1
#SB#### --ntasks-per-node=1
#
#==========================================================================
# PBS directives                qsub test_dart.csh
#
# qstat    information about the running job
# qdel     killing a job
# qsub     submitting a job
# 
#PBS -N dart_test     
#PBS -l walltime=03:00:00
#PBS -A P86850054 
#PBS -j oe
#PBS -m ae
#
# for mpi tests:
#PBS -q regular 
#PBS -l select=1:ncpus=36:mpiprocs=36
# for serial tests:
#P## -q share
#P## -l select=1:ncpus=1

#==========================================================================

set clobber
setenv MPIFLAG '-mpi'

if ( $#argv > 0 ) then
  if ( "$argv[1]" == "-mpi" ) then
    setenv MPIFLAG '-mpi'
  else if ("$argv[1]" == "-nompi") then
    setenv MPIFLAG '-nompi'
  else
    echo "Unrecognized argument to $0: $argv[1]"
    echo "Usage: $0 [ -mpi | -nompi ]"
    echo " default is to run tests using MPI."
    exit -1
  endif
endif

# set any batch system specific items here
if ($?SLURM_JOB_ID) then
  # e.g. casper
  setenv MPICMD "srun"
else if ($?PBS_NODEFILE) then
  # e.g. cheyenne
  setenv MPICMD "mpiexec_mpt"
else
  # other (no queue system, e.g. openmpi on laptop)
  setenv MPICMD "mpirun -n 2"
endif

# if your system supports different options or needs to
# use a different location for these commands, set them here.
# they will be inherited by the other test scripts.
setenv REMOVE 'rm -f'
setenv RMDIR  'rmdir'
setenv COPY   'cp -p'
setenv MOVE   'mv -f'

# require we start running this from the developer_tests dir
if ( ! -d ../models ) then
   echo "../models directory does not exist. $0 must be run from"
   echo "the developer_tests directory."
   exit 2
endif

# cd to the top level DART directory and 
# record where we are running this script
cd ..
set DARTHOME = `pwd`

if ( ! $?host) then
   setenv host `uname -n`
endif
echo "Running $0 on $host"
echo "The top-level DART directory is $DARTHOME"


#----------------------------------------------------------------------
#----------------------------------------------------------------------
# setup complete
#----------------------------------------------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------

echo
echo
echo "=================================================================="
echo "DART tests begin at "`date`
echo "=================================================================="
echo
echo

#----------------------------------------------------------------------

echo
echo
echo "=================================================================="
echo "Building and testing supported models starting at "`date`
echo "=================================================================="
echo
echo

cd ${DARTHOME}/models

./run_tests.csh $MPIFLAG -mpicmd "$MPICMD"

echo
echo
echo "=================================================================="
echo "Supported model tests complete at "`date`
echo "=================================================================="
echo
echo

#----------------------------------------------------------------------

echo
echo
echo "=================================================================="
echo "Building and testing observation converters starting at "`date`
echo "=================================================================="
echo
echo

echo "Not all observation converters are expected to build; you may"
echo "not have all the necessary supporting libraries.  So errors here"
echo "are not fatal."

cd ${DARTHOME}/observations/obs_converters

./run_tests.csh

echo
echo
echo "=================================================================="
echo "Observation converter tests complete at "`date`
echo "=================================================================="
echo
echo

#----------------------------------------------------------------------

echo
echo
echo "=================================================================="
echo "Building and testing support programs starting at "`date`
echo "=================================================================="
echo
echo

cd ${DARTHOME}/assimilation_code/programs

./run_tests.csh $MPIFLAG -mpicmd "$MPICMD"

echo
echo
echo "=================================================================="
echo "Support program tests complete at "`date`
echo "=================================================================="
echo
echo

#----------------------------------------------------------------------

echo
echo
echo "=================================================================="
echo "Building and running developer tests starting at "`date`
echo "=================================================================="
echo
echo

cd ${DARTHOME}/developer_tests

./run_tests.csh $MPIFLAG -mpicmd "$MPICMD"

echo
echo
echo "=================================================================="
echo "Developer tests complete at "`date`
echo "=================================================================="
echo
echo


#----------------------------------------------------------------------


echo
echo
echo "=================================================================="
echo "DART tests complete at "`date`
echo "=================================================================="
echo
echo

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

