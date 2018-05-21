#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# build and test all the programs given in the list.
# usage: [ -mpi | -nompi | -default ]
#
#----------------------------------------------------------------------

set usingmpi=no

if ( $#argv > 0 ) then
  if ( "$argv[1]" == "-mpi" ) then
    set usingmpi=yes
  else if ( "$argv[1]" == "-default" ) then
    set usingmpi=default
  else if ( "$argv[1]" == "-nompi" ) then
    set usingmpi=no
  else
    echo "Unrecognized argument to $0: $argv[1]"
    echo "Usage: $0 [ -mpi | -nompi | -default ]"
    echo " default is to run tests without MPI"
    exit -1
  endif
endif

# set the environment variable MPI to anything in order to enable the
# MPI builds and tests.  set the argument to the build scripts so it
# knows which ones to build.
if ( "$usingmpi" == "yes" ) then
  echo "Will be building with MPI enabled"
  set QUICKBUILD_ARG='-mpi'
else if ( "$usingmpi" == "default" ) then
  echo "Will be building with the default MPI settings"
  set QUICKBUILD_ARG=''
else if ( "$usingmpi" == "no" ) then
  echo "Will NOT be building with MPI enabled"
  set QUICKBUILD_ARG='-nompi'
else
  echo "Internal error: unrecognized value of usingmpi; should not happen"
  exit -1
endif

#----------------------------------------------------------------------

if ( ! $?REMOVE) then
   setenv REMOVE 'rm -f'
endif

if ( ! $?host) then
   setenv host `uname -n`
endif

echo "Running DART programs test on $host"

#----------------------------------------------------------------------

set programdir = `pwd`

# set the list of programs to include here

set DO_THESE_PROGRAMS = ( \
  compare_states \
  gen_sampling_err_table \
  system_simulation \
)

#----------------------------------------------------------------------
# Compile all executables for each program.
#----------------------------------------------------------------------

@ programnum = 1

foreach PROGRAM ( $DO_THESE_PROGRAMS ) 
    
    echo
    echo
    echo "=================================================================="
    echo "=================================================================="
    echo "Compiling $PROGRAM starting at "`date`
    echo "=================================================================="
    echo "=================================================================="
    echo
    echo

    cd ${programdir}/${PROGRAM}/work 
    set FAILURE = 0

    ./quickbuild.csh ${QUICKBUILD_ARG} || set FAILURE = 1

    @ programnum = $programnum + 1

    echo
    echo
    echo "=================================================================="
    echo "=================================================================="
    if ( $FAILURE ) then
      echo "ERROR - unsuccessful build of $PROGRAM at "`date`
    else
      echo "End of successful build of $PROGRAM at "`date`
    endif
    echo "=================================================================="
    echo "=================================================================="
    echo
    echo

end

echo
echo $programnum programs built.
echo

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
