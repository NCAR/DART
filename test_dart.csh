#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$


set clobber
setenv MPIFLAG '-nompi'

if ( $#argv > 0 ) then
  if ( "$argv[1]" == "-mpi" ) then
    setenv MPIFLAG '-mpi'
  else if ("$argv[1]" == "-nompi") then
    setenv MPIFLAG '-nompi'
  else
    echo "Unrecognized argument to $0: $argv[1]"
    echo "Usage: $0 [ -mpi | -nompi ]"
    echo " default is to run tests without using MPI."
    exit -1
  endif
endif

# cd to the start of the DART directory
cd ..

if ( ! -d models ) then
   echo "models does not exist. $0 must be run from the developer_tests"
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

echo "Running DART test on $host"

#----------------------------------------------------------------------
# Compile 'filter' for a wide range of models.
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
echo "Model testing complete at "`date`
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
echo "Testing location modules starting at "`date`
echo "=================================================================="
echo "=================================================================="
echo
echo

cd ${DARTHOME}/developer_tests/location
if ( 1 == 1 ) then
  ./testall.csh
endif

echo
echo
echo "=================================================================="
echo "=================================================================="
echo "Location module testing complete at "`date`
echo "=================================================================="
echo "=================================================================="
echo
echo

echo "SKIPPING Testing single-threaded lorenz_96 (L96) at "`date`
#echo "Testing single-threaded lorenz_96 (L96) at "`date`
echo "=================================================================="
echo

exit 0

echo "=================================================================="

if ! ( $?MPI ) then

   echo "MPI not enabled ... stopping."

else

   echo "No MPI tests yet ... stopping."

   #echo "=================================================================="
   #echo "testing MPI complete  at "`date`
   #echo "=================================================================="
   #echo

endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

