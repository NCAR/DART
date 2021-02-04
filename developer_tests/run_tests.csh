#!/bin/csh 
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# build and test all the models given in the list.
#
# usage: [ -mpi | -nompi ] [ -mpicmd name_of_mpi_launch_command ]
#
#----------------------------------------------------------------------

# prevent shell warning messages about no files found when trying
# to remove files using wildcards.
set nonomatch

set usingmpi=no
set MPICMD=""

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
  echo "MPI programs will be started with: $MPICMD"
else if ( "$usingmpi" == "no" ) then
  echo "Building WITHOUT MPI support."
  set QUICKBUILD_ARG='-nompi'
  set MPICMD=""
else
  echo "Internal error: unrecognized value of usingmpi; should not happen"
  exit -1
endif

# prevent shell warning messages about no files found when trying
# to remove files using wildcards.
set nonomatch

if ( ! $?REMOVE) then
   setenv REMOVE 'rm -f'
endif

if ( ! $?host) then
   setenv host `uname -n`
endif

echo
echo
echo "=================================================================="
echo "Start of DART developer_tests at "`date`
echo "=================================================================="
echo
echo
echo "Running developer_tests on $host"

set LOGDIR=`pwd`/testing_logs
${REMOVE} -r $LOGDIR
mkdir -p $LOGDIR
echo "build and run logs are in: $LOGDIR"

#----------------------------------------------------------------------

set TOPDIR = `pwd`

# collect any directory that has a quickbuild.csh script

set HAS_TESTS = `ls */work/quickbuild.csh`

#----------------------------------------------------------------------
# Compile and run all executables 
#----------------------------------------------------------------------

@ testnum = 0

foreach TESTFILE ( $HAS_TESTS ) 
    
    set TESTDIR = `dirname $TESTFILE`
    set LOGNAME = `echo $TESTDIR | sed -e 's;/[^/]*$;;' -e 's;/;_;g'`

    echo
    echo
    echo "------------------------------------------------------------------"
    echo "Compiling tests in $TESTDIR starting at "`date`
    echo "------------------------------------------------------------------"
    echo
    echo

    cd ${TESTDIR}
    set FAILURE = 0

    ( ./quickbuild.csh ${QUICKBUILD_ARG} > ${LOGDIR}/buildlog.${LOGNAME}.out ) || set FAILURE = 1

    @ testnum = $testnum + 1

    echo
    echo
    if ( $FAILURE ) then
      echo "------------------------------------------------------------------"
      echo "ERROR - unsuccessful build in $TESTDIR at "`date`
      echo "------------------------------------------------------------------"
      cd $TOPDIR
      continue
    else
      echo "------------------------------------------------------------------"
      echo "Running tests in $TESTDIR"
      echo "------------------------------------------------------------------"
      echo
      echo
  
      foreach TARGET ( mkmf_* )
  

           set FAILURE = 0
           set PROG = `echo $TARGET | sed -e 's#mkmf_##'`
           set THISLOG = ${LOGDIR}/runlog.${LOGNAME}.${PROG}.out 
           echo Starting $PROG

           if ( -f using_mpi_for_$PROG && -f ${PROG}.in ) then
              ( ${MPICMD} ./$PROG < ${PROG}.in >> $THISLOG ) || set FAILURE = 1
           else if ( -f using_mpi_for_$PROG ) then
              ( ${MPICMD} ./$PROG              >> $THISLOG ) || set FAILURE = 1
           else if ( -f ${PROG}.in ) then
              (           ./$PROG < ${PROG}.in >> $THISLOG ) || set FAILURE = 1
           else
              (           ./$PROG              >> $THISLOG ) || set FAILURE = 1
           endif
         
           if ( $FAILURE ) then
              echo "ERROR - unsuccessful run of $PROG at "`date`

              switch ( $PROG )
                 case stacktest
                    echo "        stacktest intentionally fails."
                    echo "        stacktest continually reallocates a stack array until it fails."
                 breaksw
                 case test_diag_structure
                    echo "        test_diag_structure intentionally fails."
                    echo "        test_diag_structure tries to allocate 11 domains - max is 10."
                 breaksw
                 case test_state_structure
                    echo "        test_state_structure fails."
                    echo "        test_state_structure not fully tested."
                 breaksw
                 default
                    echo "unexpected error"
                 breaksw
              endsw

           else
              ${REMOVE} $PROG
              echo "Successful run of $PROG"
           endif
  
      end

      ${REMOVE} Makefile input.nml.*_default .cppdefs *.o *.mod

      cd $TOPDIR

    endif

    echo
    echo
    echo "------------------------------------------------------------------"
    echo "Done running tests in $TESTDIR at "`date`
    echo "------------------------------------------------------------------"
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

( ./run_tests.csh > ${LOGDIR}/location_tests.out ) || set FAILURE = 1

echo
echo
echo "=================================================================="
echo "Done running location tests at "`date`
echo "=================================================================="
echo
echo

cd $TOPDIR

echo
echo "$testnum developer_tests run."
echo
echo "=================================================================="
echo "End of DART developer_tests at "`date`
echo "=================================================================="

exit 0
