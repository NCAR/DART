#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

set SNAME = $0
set clobber

switch ( $#argv )
   case 0:
      # supplying no args - ok
      breaksw
      
   default:
      echo " "
      echo "usage: $SNAME:t"
      echo " "
      echo "This script compiles all available models for this branch of the DART code."
      echo "This must be run from the top-level 'DART' directory."
      exit 1
      breaksw
endsw

if ( ! -d models ) then
   echo "models does not exist. $SNAME:t must be run from the top-level"
   echo "DART directory -- please try again."
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
# set the environment variable MPI to anything in order to enable the
# MPI builds and tests.  set the argument to the build scripts so it
# knows which ones to build.
if ( $?MPI ) then
  echo "Will be building with MPI enabled"
  setenv QUICKBUILD_ARG -mpi
else
  echo "Will NOT be building with MPI enabled"
  setenv QUICKBUILD_ARG -nompi
endif
#----------------------------------------------------------------------

echo "Running DART test on $host"

#----------------------------------------------------------------------
# Compile 'filter' for a wide range of models.
# CAVEATS:
#    The PBL_1d model relies on routines that require special compiler
#    flags to be recognized as F90 code. Since these are compiler-specific,
#    I have not figured out a way yet to automate this. 
#
#----------------------------------------------------------------------

@ modelnum = 10

if ( 1 == 1 ) then
foreach MODEL ( \
  bgrid_solo \
  cam-fv \
  cice \
  cm1 \
  lorenz_63 \
  lorenz_96 \
  mpas_atm \
  POP \
  ROMS \
  wrf )
  # many models intentionally omitted.
    
    echo "=================================================================="
    echo "Compiling $MODEL at "`date`
    echo ""

    cd ${DARTHOME}/models/${MODEL}/work

    ./quickbuild.csh ${QUICKBUILD_ARG} || exit 3

    echo "Trying to run pmo for model $MODEL as a test"
    ./perfect_model_obs

    echo "Removing the newly-built objects ..."
    ${REMOVE} *.o *.mod 
    ${REMOVE} Makefile input.nml.*_default .cppdefs
    foreach TARGET ( mkmf_* )
      set PROG = `echo $TARGET | sed -e 's#mkmf_##'`
      ${REMOVE} $PROG
    end

    @ modelnum = $modelnum + 1
end
endif

echo
echo
echo
echo "=================================================================="
echo "Testing observation converters at "`date`
echo "=================================================================="
echo

echo "Not all observation converters are expected to build; you may"
echo "not have all the necessary supporting libraries.  So errors here"
echo "are not fatal."

cd ${DARTHOME}/observations
if ( 1 == 1 ) then
  ./buildall.csh
endif

echo
echo "=================================================================="
echo "Observation converter testing complete at "`date`
echo "=================================================================="
echo

echo
echo
echo
echo "=================================================================="
echo "Testing location modules at "`date`
echo "=================================================================="
echo

cd ${DARTHOME}/location
if ( 1 == 1 ) then
  ./testall.csh
endif

echo
echo "=================================================================="
echo "Location module testing complete at "`date`
echo "=================================================================="
echo



echo
echo
echo
echo "=================================================================="
echo "Testing single-threaded bgrid_solo at "`date`
echo "=================================================================="
echo

set MODEL = bgrid_solo

cd ${DARTHOME}/models/${MODEL}/work

# Save the 'original' files so we can reinstate them as needed

${COPY} input.nml     input.nml.$$
${COPY} obs_seq.in   obs_seq.in.$$
${COPY} perfect_ics perfect_ics.$$
${COPY} filter_ics   filter_ics.$$

# Begin by compiling all programs; need to stop if an error is detected
./quickbuild.csh -nompi || exit 91

 
# Run the perfect model and the filter
./perfect_model_obs  || exit 92
./filter             || exit 93

echo "Removing the newly-built objects ..."
${REMOVE} filter_restart perfect_restart
${REMOVE} input.nml perfect_ics filter_ics
${REMOVE} obs_seq.in obs_seq.out obs_seq.final
${REMOVE} True_State.nc Prior_Diag.nc Posterior_Diag.nc
${REMOVE} *.o *.mod 
${REMOVE} Makefile input.nml.*_default .cppdefs
foreach TARGET ( mkmf_* )
  set PROG = `echo $TARGET | sed -e 's#mkmf_##'`
  ${REMOVE} $PROG
end

# Reinstate the 'original' files so we can run this again if we need to.

${MOVE}   input.nml.$$   input.nml
${MOVE}  obs_seq.in.$$ obs_seq.in
${MOVE} perfect_ics.$$ perfect_ics
${MOVE}  filter_ics.$$  filter_ics

echo
echo "=================================================================="
echo "Single-threaded testing of bgrid_solo complete at "`date`
echo "=================================================================="
echo

echo
echo "=================================================================="
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

