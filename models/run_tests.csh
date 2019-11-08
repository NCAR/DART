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

# set the quickbuild argument
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

echo "Running DART model test on $host"

#----------------------------------------------------------------------

set modeldir = `pwd`

# set the list of models to include here

set DO_THESE_MODELS = ( \
  9var \
  POP \
  ROMS \
  FESOM \
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
  noah \
  null_model \
  simple_advection \
  template \
  wrf \
)

#----------------------------------------------------------------------
# either run the workshop setup or quickbuild/run then clean
#---------------------------------------------------------------------

echo
echo
echo "=================================================================="
echo "Starting tests of model directory at "`date`
echo "=================================================================="
echo
echo

${REMOVE} ${LOGDIR}/buildlog.*.out ${LOGDIR}/runlog.*.out
mkdir -p ${LOGDIR}
echo "build and run logs are in: $LOGDIR"


@ modelnum = 0

foreach MODEL ( $DO_THESE_MODELS ) 
    
    echo
    echo
    echo "=================================================================="
    echo "Testing $MODEL starting at "`date`
    echo "=================================================================="
    echo
    echo

    cd ${modeldir}/${MODEL}/work
    set FAILURE = 0
    echo "Current directory is " `pwd`

    @ ncdlfiles = `ls *.cdl | wc -l`

    if ( "$MODEL" == "template" ) then
      echo "skipping tests of the template directory"
      continue
    endif

    # save original input.nml & obs seq files here
    set SAVEDIR = saveme.test_dart
    mkdir -p ${SAVEDIR}
    ${COPY} input.nml obs_seq.* ${SAVEDIR}

    # If there is a testing namelist, use it.
    if ( -f input.nml.testing ) then
       ${COPY} input.nml.testing input.nml
    endif

    if ( -f workshop_setup.csh ) then

      echo "Trying to run workshop_setup.csh for model $MODEL as a test"
      ( ./workshop_setup.csh >  ${LOGDIR}/buildlog.${MODEL}.out ) || set FAILURE = 1
      echo
      echo "Re-running workshop_setup.csh to test overwriting files for model $MODEL"
      ( ./workshop_setup.csh >> ${LOGDIR}/buildlog.${MODEL}.out ) || set FAILURE = 1
      echo

    else
      echo "building executables for $MODEL"

      ( ./quickbuild.csh ${QUICKBUILD_ARG} > ${LOGDIR}/buildlog.${MODEL}.out ) || set FAILURE = 1
      echo

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
      # assumes the executables from quickbuild are here
      if ( -f using_mpi_for_perfect_model_obs ) then
         ( $MPICMD ./perfect_model_obs >  ${LOGDIR}/runlog.${MODEL}.out ) || set FAILURE = 1
         echo "Rerunning PMO to test for output file overwrite"
         ( $MPICMD ./perfect_model_obs >> ${LOGDIR}/runlog.${MODEL}.out ) || set FAILURE = 1
      else
         (         ./perfect_model_obs >  ${LOGDIR}/runlog.${MODEL}.out ) || set FAILURE = 1
         echo "Rerunning PMO to test for output file overwrite"
         (         ./perfect_model_obs >> ${LOGDIR}/runlog.${MODEL}.out ) || set FAILURE = 1
      endif
      # FIXME: if possible, try running filter here as well?
    endif

    if ( -f model_mod_check ) then
      echo "Trying to run model_mod_check for model $MODEL as a test"
      if ( -f using_mpi_for_model_mod_check ) then
         ( $MPICMD ./model_mod_check >> ${LOGDIR}/runlog.${MODEL}.out ) || set FAILURE = 1 
      else
         (         ./model_mod_check >> ${LOGDIR}/runlog.${MODEL}.out ) || set FAILURE = 1 
      endif
    endif

    echo "Removing the newly-built objects and executables"
    ${REMOVE} *.o *.mod 
    ${REMOVE} Makefile input.nml.*_default .cppdefs
    foreach TARGET ( mkmf_* )
      set PROG = `echo $TARGET | sed -e 's#mkmf_##'`
      ${REMOVE} $PROG
    end

    echo "Restoring original input.nml and obs_seq files"
    ${MOVE} ${SAVEDIR}/* .
    ${REMOVE_DIR} ${SAVEDIR}
  
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
    if ( $FAILURE ) then
      echo "ERROR - unsuccessful test of $MODEL at "`date`

      switch ( $MODEL )
         case FESOM
            echo "Note that because the FESOM-native code explicitly types reals,"
            echo "the DART mechanism of being able to run in reduced precision by"
            echo "defining real(r8) to be the same as real(r4) via 'types_mod.f90'"
            echo "is not supported. Please check to make sure this is the reason"
            echo "this test is failing."
         breaksw
         case clm
            echo "CLM is expected to fail on this branch."
         breaksw
         default
            echo "unexpected error"
         breaksw
      endsw

    else
      echo "End of succesful test of $MODEL at "`date`
    endif
    echo "=================================================================="
    echo
    echo

end

echo
echo "$modelnum models tested."
echo

echo
echo
echo "=================================================================="
echo "Ending tests of model directory at "`date`
echo "=================================================================="
echo
echo
exit 0

