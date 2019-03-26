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
  echo "Building with MPI support."
  set QUICKBUILD_ARG='-mpi'
else if ( "$usingmpi" == "default" ) then
  echo "Building with the default MPI settings"
  set QUICKBUILD_ARG=''
else if ( "$usingmpi" == "no" ) then
  echo "Building WITHOUT MPI support."
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
if ( ! $?MOVE) then
   setenv MOVE 'mv -f'
endif

if ( ! $?host) then
   setenv host `uname -n`
endif

echo "Running DART model test on $host"

#----------------------------------------------------------------------

set modeldir = `pwd`

# set the list of models to include here

set DO_THESE_MODELS = ( \
  9var \
  POP \
  ROMS \
  bgrid_solo \
  cam-fv \
  cice \
  clm \
  cm1 \
  forced_lorenz_96 \
  ikeda \
  lorenz_04 \
  lorenz_63 \
  lorenz_84 \
  lorenz_96 \
  lorenz_96_2scale \
  mpas_atm \
  null_model \
  simple_advection \
  template \
  wrf \
)

#----------------------------------------------------------------------
# Compile all executables for each model.
#----------------------------------------------------------------------

@ modelnum = 1

foreach MODEL ( $DO_THESE_MODELS ) 
    
    echo
    echo
    echo "=================================================================="
    echo "=================================================================="
    echo "Compiling $MODEL starting at "`date`
    echo "=================================================================="
    echo "=================================================================="
    echo
    echo

    cd ${modeldir}/${MODEL}/work 
    set FAILURE = 0

    ./quickbuild.csh ${QUICKBUILD_ARG} || set FAILURE = 1

    @ modelnum = $modelnum + 1

    echo
    echo
    echo "=================================================================="
    echo "=================================================================="
    if ( $FAILURE ) then
      echo "ERROR - unsuccessful build of $MODEL at "`date`
    else
      echo "End of successful build of $MODEL at "`date`
    endif
    echo "=================================================================="
    echo "=================================================================="
    echo
    echo

end

#----------------------------------------------------------------------
# Run PMO and filter if possible.  Save and restore the original input.nml.
#----------------------------------------------------------------------

@ modelnum = 1

foreach MODEL ( $DO_THESE_MODELS ) 
    
    echo
    echo
    echo "=================================================================="
    echo "=================================================================="
    echo "Testing $MODEL starting at "`date`
    echo "=================================================================="
    echo "=================================================================="
    echo
    echo
 
    cd ${modeldir}/${MODEL}/work
    set FAILURE = 0
    echo "Current directory is " `pwd`

    @ ncdlfiles = `ls *.cdl | wc -l`

    if ( "$MODEL" == "template" ) then
      echo skipping tests of the template directory

    else 
      # save original input.nml & obs seq files here
      set SAVEDIR = saveme.test_dart
      mkdir -p ${SAVEDIR}
      ${COPY} input.nml obs_seq.* ${SAVEDIR}

      if ( -f workshop_setup.csh ) then
  
        echo "Trying to run workshop_setup.csh for model $MODEL as a test"
        ./workshop_setup.csh || set FAILURE = 1
        echo "Re-running workshop_setup.csh to test overwriting files for model $MODEL"
        ./workshop_setup.csh || set FAILURE = 1
  
      else
        echo "Trying to run pmo for model $MODEL as a test"
        echo "Will generate NetCDF files from any .cdl files first."
        # try not to error out if no .cdl files found
        if ( $ncdlfiles > 0 ) then
           foreach i ( *.cdl )
             set base = `basename $i .cdl`
             if ( -f ${base}.nc ) continue
             ncgen -o ${base}.nc $i
           end
        endif
        # assumes the executables from the first pass are still here
        ./perfect_model_obs || set FAILURE = 1
        echo "Rerunning PMO to test for file overwrite"
        ./perfect_model_obs || set FAILURE = 1
        # if possible, try running filter here as well?
      endif
  
      if ( -f model_mod_check ) then
        echo "Trying to run model_mod_check for model $MODEL as a test"
        ./model_mod_check || set FAILURE = 1 
      endif
    endif

    echo "Removing the newly-built objects"
    ${REMOVE} *.o *.mod 
    ${REMOVE} Makefile input.nml.*_default .cppdefs
    foreach TARGET ( mkmf_* )
      set PROG = `echo $TARGET | sed -e 's#mkmf_##'`
      ${REMOVE} $PROG
    end

    echo "Restoring original input.nml and obs_seq files"
    ${MOVE} ${SAVEDIR}/* .
    ${REMOVE} ${SAVEDIR}
  
    if ( $ncdlfiles > 0 ) then
      foreach i ( *.cdl )
        set base = `basename $i .cdl`
        if ( -f ${base}.nc ) rm ${base}.nc
      end
    endif

    @ modelnum = $modelnum + 1

    echo
    echo
    echo "=================================================================="
    echo "=================================================================="
    if ( $FAILURE ) then
      echo "ERROR - unsuccessful test of $MODEL at "`date`
    else
      echo "End of succesful test of $MODEL at "`date`
    endif
    echo "=================================================================="
    echo "=================================================================="
    echo
    echo

end

echo
echo $modelnum models tested.
echo

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
