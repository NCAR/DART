#!/bin/csh 
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# build and test all the models given in the list.
#
# usage: [ -mpi | -nompi ] [ -mpicmd name_of_mpi_launch_command ]
#
#----------------------------------------------------------------------

set usingmpi=no
set MPICMD=""
set LOGDIR=`pwd`/testing_logs

if ( $#argv > 0 ) then
  if ( "$argv[1]" == "-mpi" ) then
    set usingmpi=yes
  else if ( "$argv[1]" == "-nompi" ) then
    set usingmpi=no
  else
    echo "Unrecognized argument to $0: $argv[1]"
    echo "Usage: $0 [ -mpi | -nompi ]  [ -mpicmd name_of_mpi_launch_command ]"
    echo " default is to run tests without MPI"
    exit -1
  endif
  shift
endif

if ( $#argv > 1 ) then
  if ( "$argv[1]" == "-mpicmd" ) then
    set MPICMD = "$argv[2]"
  else
    echo "Unrecognized argument to $0: $argv[1]"
    echo "Usage: $0 [ -mpi | -nompi ]  [ -mpicmd name_of_mpi_launch_command ]"
    echo " default is to run tests without MPI"
    exit -1
  endif
  shift
endif

# set the environment variable MPI to anything in order to enable the
# MPI builds and tests.  set the argument to the build scripts so it
# knows which ones to build.
if ( "$usingmpi" == "yes" ) then
  echo "Building with MPI support."
  set QUICKBUILD_ARG='-mpi'
  if ( ! $?MPICMD) then
    set MPICMD='mpirun -n 2'
  endif
else if ( "$usingmpi" == "no" ) then
  echo "Building WITHOUT MPI support."
  set QUICKBUILD_ARG='-nompi'
  set MPICMD=""
else
  echo "Internal error: unrecognized value of usingmpi; should not happen"
  exit -1
endif

#----------------------------------------------------------------------

if ( ! $?REMOVE) then
   setenv REMOVE 'rm -f'
endif
if ( ! $?REMOVE_DIR) then
   setenv REMOVE_DIR 'rmdir'
endif
if ( ! $?COPY) then
   setenv COPY 'cp -f'
endif
if ( ! $?MOVE) then
   setenv MOVE 'mv -f'
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

mkdir -p $LOGDIR
\rm -f $LOGDIR/*
echo putting build and run logs in $LOGDIR

@ testnum = 0

foreach TESTFILE ( $HAS_TESTS ) 
    
    set TESTDIR = `dirname $TESTFILE`

    echo
    echo
    echo "=================================================================="
    echo "Compiling tests in $TESTDIR starting at "`date`
    echo "=================================================================="
    echo
    echo

    cd ${TESTDIR}
    set FAILURE = 0

    ( ./quickbuild.csh ${QUICKBUILD_ARG} > ${LOGDIR}/buildlog.${TESTDIR}.out ) || set FAILURE = 1

    @ testnum = $testnum + 1

    echo
    echo
    if ( $FAILURE ) then
      echo "=================================================================="
      echo "ERROR - unsuccessful build in $TESTDIR at "`date`
      echo "=================================================================="
      cd $TOPDIR
      continue
    else
      echo "=================================================================="
      echo "End of successful build in $TESTDIR at "`date`
      echo "=================================================================="
      echo
      echo
      echo "=================================================================="
      echo "Running tests in $TESTDIR starting at "`date`
      echo "=================================================================="
      echo
      echo
  
      foreach TARGET ( mkmf_* )
  
           \rm -f *.o *.mod
           \rm -f Makefile input.nml.*_default .cppdefs

           set PROG = `echo $TARGET | sed -e 's#mkmf_##'`
           echo "++++++++++++++++++"
           echo Starting $PROG
           ( ${MPICMD} ./$PROG  > ${LOGDIR}/runlog.${TESTDIR}.out ) || set FAILURE = 1
           echo Finished $PROG
           echo
           if ( $FAILURE ) then
              echo "ERROR - unsuccessful run in $TESTDIR at "`date`
           else
              \rm -f $PROG
           endif
           echo "++++++++++++++++++"
  
      end

      cd $TOPDIR

    endif

    echo
    echo
    echo "=================================================================="
    echo "Done running tests in $TESTDIR at "`date`
    echo "=================================================================="
    echo
    echo

end

#----------------------------------------------------------------------
# Compile and run all location tests
#----------------------------------------------------------------------

# special for locations
cd $TOPDIR/location

echo "=================================================================="
echo "Running location tests starting at "`date`
echo "=================================================================="

set FAILURE = 0

( ./testall.csh > ${LOGDIR}/location_tests.out ) || set FAILURE = 1

echo
echo
echo "=================================================================="
echo "Done running location tests at "`date`
echo "=================================================================="
echo
echo

cd $TOPDIR


echo
echo $testnum developer tests run.
echo

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
