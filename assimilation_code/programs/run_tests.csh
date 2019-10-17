#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# build and test all the programs given in the list.
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

echo "Running DART programs test on $host"

#----------------------------------------------------------------------

set programdir = `pwd`

# set the list of programs to include here

# FIXME: note that an important set of programs has no testing done
# on them yet.  programs currently in this directory which aren't
# run otherwise in any test scripts include:
#
#  closest_member_tool
#  compare_states
#  compute_error
#  fill_inflation_restart
#  gen_sampling_err_table
#  integrate_model
#  obs_assim_count
#  obs_common_subset
#  obs_diag
#  obs_impact_tool
#  obs_keep_a_few
#  obs_loop
#  obs_selection
#  obs_seq_coverage
#  obs_seq_to_netcdf
#  obs_sequence_tool
#  obs_seq_verify
#  obs_total_error
#  obs_utils
#  perturb_single_instance


# expand these tests.
set DO_THESE_PROGRAMS = ( \
  compare_states \
  system_simulation \
)

#----------------------------------------------------------------------
# Compile all executables for each program.
#----------------------------------------------------------------------

echo
echo
echo "=================================================================="
echo "Starting tests of dart programs at "`date`
echo "=================================================================="
echo
echo

mkdir -p $LOGDIR
\rm -f $LOGDIR/*
echo putting build and run logs in $LOGDIR

@ programnum = 0

foreach PROGRAM ( $DO_THESE_PROGRAMS ) 
    
    echo
    echo
    echo "=================================================================="
    echo "Compiling $PROGRAM starting at "`date`
    echo "=================================================================="
    echo
    echo

    cd ${programdir}/${PROGRAM}/work 
    set FAILURE = 0

    ( ./quickbuild.csh ${QUICKBUILD_ARG} > ${LOGDIR}/buildlog.$PROGRAM.out ) || set FAILURE = 1

    @ programnum = $programnum + 1

    echo
    echo
    if ( $FAILURE ) then
      echo "=================================================================="
      echo "ERROR - unsuccessful build of $PROGRAM at "`date`
      switch ( $PROGRAM )
         case obs_sampling_err
            echo " This build expected to fail if running in reduced precision"
            echo " by defining r8 same as r4. If this is not the case you are"
            echo " testing, you have other problems."
            echo
         breaksw
         default
            echo " unexpected error"
         breaksw
      endsw
      echo "=================================================================="
      echo
      echo
      continue
    else
      echo "=================================================================="
      echo "End of successful build of $PROGRAM at "`date`
      echo "=================================================================="
      echo
      echo

# FIXME add this back when we can successfully run tests on these programs

#      echo
#      echo
#      echo "=================================================================="
#      echo "Running tests for $PROGRAM starting at "`date`
#      echo "=================================================================="
#      echo
#      echo
#  
#      foreach TARGET ( mkmf_* )
#  
#           \rm -f *.o *.mod
#           \rm -f Makefile input.nml.*_default .cppdefs
#
            set FAILURE = 0
#           set PROG = `echo $TARGET | sed -e 's#mkmf_##'`
#           echo "++++++++++++++++++"
#           echo Starting $PROG
#           ( ${MPICMD} ./$PROG  > ${LOGDIR}/runlog.$PROG.out ) || set FAILURE = 1
#           echo Finished $PROG
#           echo "++++++++++++++++++"
#           echo
#
#           \rm -f $PROG
#  
#      end
#
#      echo
#      echo
#      echo "=================================================================="
#      echo "Done with tests of $PROGRAM at "`date`
#      echo "=================================================================="
#      echo
#      echo

    endif


end

echo
echo $programnum programs built.
echo

echo
echo
echo "=================================================================="
echo "Ending tests of dart programs at "`date`
echo "=================================================================="
echo
echo
exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
