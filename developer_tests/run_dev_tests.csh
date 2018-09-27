#!/bin/csh 
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# build and test all the models given in the list.
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
if ( ! $?COPY) then
   setenv COPY 'cp -f'
endif

if ( ! $?host) then
   setenv host `uname -n`
endif

echo "Running DART developer tests on $host"

#----------------------------------------------------------------------

set TOPDIR = `pwd`

# collect any directory that has a quickbuild.csh script

set HAS_TESTS = `ls */work/quickbuild.csh`

#----------------------------------------------------------------------
# Compile and run all executables 
#----------------------------------------------------------------------

@ testnum = 1

foreach TESTFILE ( $HAS_TESTS ) 
    
    set TESTDIR = `dirname $TESTFILE`

    echo
    echo
    echo "=================================================================="
    echo "=================================================================="
    echo "Compiling tests in $TESTDIR starting at "`date`
    echo "=================================================================="
    echo "=================================================================="
    echo
    echo

    cd ${TESTDIR}
    set FAILURE = 0

    ./quickbuild.csh ${QUICKBUILD_ARG} || set FAILURE = 1

    @ testnum = $testnum + 1

    echo
    echo
    echo "=================================================================="
    echo "=================================================================="
    if ( $FAILURE ) then
      echo "ERROR - unsuccessful build in $TESTDIR at "`date`
      cd $TOPDIR
      cycle
    else
      echo "End of successful build in $TESTDIR at "`date`
    endif
    echo "=================================================================="
    echo "=================================================================="
    echo
    echo

    echo
    echo
    echo "=================================================================="
    echo "=================================================================="
    echo "Running tests in $TESTDIR starting at "`date`
    echo "=================================================================="
    echo "=================================================================="
    echo
    echo

    foreach TARGET ( mkmf_* )

         set PROG = `echo $TARGET | sed -e 's#mkmf_##'`
         echo Starting $PROG
         ./$PROG 
         echo Finished $PROG

    end

    echo
    echo
    echo "=================================================================="
    echo "=================================================================="
    echo "Done running tests in $TESTDIR at "`date`
    echo "=================================================================="
    echo "=================================================================="
    echo
    echo

    cd $TOPDIR

end


# special for locations
cd $TOPDIR/location
if ( 1 == 1 ) then
  ./testall.csh
endif
cd $TOPDIR


echo
echo $testnum developer tests run.
echo

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
