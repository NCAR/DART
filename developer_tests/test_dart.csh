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
#>@todo FIXME ... implement some method to run mpi executables. 
#
#==========================================================================
# SLURM directives              sbatch test_batch.csh
#
# sinfo     information about the whole slurm system
# squeue    information about running jobs
# sbatch    submitting a job
# scancel   killing a job
#
#SBATCH --ignore-pbs
#SBATCH --job-name dart_test
#SBATCH --ntasks-per-node=1
#SBATCH -t 2:00:00
#SBATCH -n 1
#SBATCH -A P86850054
#SBATCH -p dav
#SBATCH -o dart_test.log
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#
#==========================================================================
# PBS directives                qsub test_batch.csh
#
#PBS -N dart_test     
#PBS -l walltime=02:00:00
#PBS -q share 
#PBS -l select=1:ncpus=1
#PBS -A P86850054 
#PBS -j oe
#PBS -m abe

set clobber
setenv MPIFLAG '-nompi'

if ( $#argv > 0 ) then
  if ( "$argv[1]" == "-mpi" ) then
    setenv MPIFLAG '-mpi'
  else if ("$argv[1]" == "-nompi") then
    setenv MPIFLAG '-nompi'
  else if ("$argv[1]" == "-default") then
    setenv MPIFLAG '-default'
  else
    echo "Unrecognized argument to $0: $argv[1]"
    echo "Usage: $0 [ -mpi | -nompi | -default ]"
    echo " default is to run tests without using MPI."
    exit -1
  endif
endif

# Since this script does not launch the mpi executables with an mpirun etc.,
# this whole script can only possibly test the serial implementation.
# So - in batch mode, just force it to be nompi.
#>@todo FIXME ... implement some method to run mpi executables. 

if ($?SLURM_JOB_ID) then
    setenv MPIFLAG '-nompi'
else if ($?PBS_NODEFILE) then
    setenv MPIFLAG '-nompi'
endif

# cd to the start of the DART directory
cd ..

if ( ! -d models ) then
   echo "models directory does not exist. $0 must be run from the developer_tests"
   echo "directory -- please try again."
   exit 2
else
   set DARTHOME = `pwd`
endif

echo "The top-level DART directory (DARTHOME) is $DARTHOME"

#----------------------------------------------------------------------
# See if some necessary environment variables are set.
# We'd like to have a short hostname but uname can be configured very
# differently from host to host.
#----------------------------------------------------------------------

if ( ! $?host) then
   setenv host `uname -n`
endif
echo "Running $0 on $host"

#----------------------------------------------------------------------
# Not all unix systems support the same subset of flags; try to figure
# out what system we are running on and adjust accordingly.
#----------------------------------------------------------------------

set OSTYPE = `uname -s`
switch ( ${OSTYPE} )
   case IRIX64:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      breaksw
   case AIX:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      breaksw
   default:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -vp'
      setenv   MOVE 'mv -fv'
      breaksw
endsw

#----------------------------------------------------------------------

echo
echo
echo "=================================================================="
echo "=================================================================="
echo "Building and testing supported models starting at "`date`
echo "=================================================================="
echo "=================================================================="
echo
echo

cd ${DARTHOME}/models
if ( 1 == 1 ) then
  ./buildall.csh $MPIFLAG
endif

echo
echo
echo "=================================================================="
echo "=================================================================="
echo "Supported model testing complete at "`date`
echo "=================================================================="
echo "=================================================================="
echo
echo

#----------------------------------------------------------------------

echo
echo
echo "=================================================================="
echo "=================================================================="
echo "Testing observation converters starting at "`date`
echo "=================================================================="
echo "=================================================================="
echo
echo

echo "Not all observation converters are expected to build; you may"
echo "not have all the necessary supporting libraries.  So errors here"
echo "are not fatal."

cd ${DARTHOME}/observations/obs_converters
if ( 1 == 1 ) then
  ./buildall.csh
endif

echo
echo
echo "=================================================================="
echo "=================================================================="
echo "Observation converter testing complete at "`date`
echo "=================================================================="
echo "=================================================================="
echo
echo

#----------------------------------------------------------------------

echo
echo
echo "=================================================================="
echo "=================================================================="
echo "Building and testing support programs starting at "`date`
echo "=================================================================="
echo "=================================================================="
echo
echo

cd ${DARTHOME}/assimilation_code/programs
if ( 1 == 1 ) then
  ./buildall.csh $MPIFLAG
endif

echo
echo
echo "=================================================================="
echo "=================================================================="
echo "Support program testing complete at "`date`
echo "=================================================================="
echo "=================================================================="
echo
echo

#----------------------------------------------------------------------

echo
echo
echo "=================================================================="
echo "=================================================================="
echo "Running developer tests starting at "`date`
echo "=================================================================="
echo "=================================================================="
echo
echo

cd ${DARTHOME}/developer_tests
if ( 1 == 1 ) then
  ./run_dev_tests.csh
endif

echo
echo
echo "=================================================================="
echo "=================================================================="
echo "Developer tests complete at "`date`
echo "=================================================================="
echo "=================================================================="
echo
echo

#----------------------------------------------------------------------

if ( $?MPI ) then

   #echo
   #echo
   #echo "=================================================================="
   #echo "=================================================================="
   #echo "MPI testing complete  at "`date`
   #echo "=================================================================="
   #echo "=================================================================="
   #echo
   #echo

   echo "No MPI tests yet ... stopping."

   #echo
   #echo
   #echo "=================================================================="
   #echo "=================================================================="
   #echo "MPI testing complete  at "`date`
   #echo "=================================================================="
   #echo "=================================================================="
   #echo
   #echo

else

   echo "MPI tests not requested ... stopping."

endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

