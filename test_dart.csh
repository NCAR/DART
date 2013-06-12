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
      # supplying no args - default is 'TestDir'
      # The results of each run-time test will be stored in this directory.
      # This will facilitate checks across platforms.

      set BASEOUTPUTDIR = TestDir
      breaksw
   case 1:
      # supplying one argument -- the base output directory.
      # The results of each run-time test will be stored in this directory.
      # This will facilitate checks across platforms.
      set BASEOUTPUTDIR = $1

      breaksw
   default:
      echo " "
      echo "usage: $SNAME:t OutputDirectory"
      echo " "
      echo "This script compiles ?every? program unit for a wide range of models and then"
      echo "does relatively extensive tests of the L96 programs with a variety of options."
      echo "The L96 tests are best run from a 'clean' starting point - i.e. one that"
      echo "is as close to the distribution state as possible. Picking up in the middle"
      echo "is not particularly easy. I always run the script end-to-end. TJH"
      echo " "
      echo "An attempt is made to preserve the state of the input.nml, obs_seq.xx?, and"
      echo "some initial conditions files. If the script completes without errors, your"
      echo "original files are reinstated. If not, your lorenz_96/work directory is a mess."
      echo " "
      echo "This must be run from the top-level 'DART' directory."
      echo " "
      echo "The OutputDirectory will contain the output of the L96 runs and is "
      echo "intended to be useful for comparing results from multiple compilers/platforms."
      echo "If the directory exists, it is guaranteed that the contents will be removed ..."
      echo "If you specify a relative filename for the directory, look in "
      echo "DART/models/lorenz_96/work/ ..."
      echo " "
      echo "Some of the tests should be bit-wise reproducible -- if this fails ... "
      echo "$SNAME:t will abort."
      echo " "
      echo "This is a pretty verbose process, so if you are logging the output,"
      echo "make sure you have plenty of space:"
      echo " "
      echo "./$SNAME:t TestDir |& tee DART_test.log"
      echo " "
      echo "can easily result in a 750 Kb log file"
      exit 0
      breaksw
endsw

if ( ! -d models/lorenz_96 ) then
   echo "models/lorenz_96 does not exist. $SNAME:t must be run from the top-level"
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
# We need to run the editor in batch mode.  If you have 'vim' it needs
# one flag; if you have only the older vanilla 'vi' you need another.
# On several systems 'vi' is a link to 'vim' and uses the newer syntax
# so you cannot distinguish which flag will be needed based only on name.
# First try to run 'vim' by full name and then back off to plain 'vi' 
# if it is not found.  Punt if neither is found.
#---------------------------------------------------------------------- 
set VI_EXE = `which vim` 
if ( -x "${VI_EXE}" ) then
   setenv VI 'vim -e'
else
   set VI_EXE = `which vi` 
   if ( -x "${VI_EXE}" ) then
      setenv VI 'vi -s'
   else
      echo "Neither the vim nor the vi editor were found.  This script"
      echo "cannot continue unless it can use one of them to update"
      echo "the test input namelist files."
      exit 2
   endif
endif

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
  9var \
  am2 \
  bgrid_solo \
  cam \
  #cosmo \
  forced_lorenz_96 \
  #gitm \
  ikeda \
  lorenz_04 \
  lorenz_63 \
  lorenz_84 \
  lorenz_96 \
  lorenz_96_2scale \
  MITgcm_annulus \
  MITgcm_ocean \
  mpas_atm \
  mpas_ocn \
  NAAPS \
  NCOMMAS \
  null_model \
  #PBL_1d \
  pe2lyr \
  POP \
  #rose \
  simple_advection \
  template \
  tiegcm \
  wrf )
  # intentionally omitted:
  #  forced_barot MITgcm_annulus PBL_1d rose
    
    echo "=================================================================="
    echo "Compiling $MODEL at "`date`
    echo ""

    cd ${DARTHOME}/models/${MODEL}/work

    ./quickbuild.csh ${QUICKBUILD_ARG} || exit 3

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

#----------------------------------------------------------------------
# Lots of tests for L96, and one for bgrid_solo (for the 3d loc stuff)
#----------------------------------------------------------------------

cd ${DARTHOME}

echo
echo
echo
echo "=================================================================="
echo "Checking matlab support ... "
echo "=================================================================="
echo

# If matlab and the netcdf functions exist, generate figures on-the-fly.
# If not, archive the data.

which matlab > /dev/null
set MatlabExists=$status

if ( $MatlabExists == 0 ) then

   set fname = dart.matlabcheck.$$
   touch $fname

   echo "fname = '"$fname"';"        >! batchscript.m
   echo "addpath ${DARTHOME}/matlab" >> batchscript.m 
   echo "ChecknetCDFuse(fname);"     >> batchscript.m
   echo "quit"                       >> batchscript.m

   matlab -nosplash -nodesktop -r batchscript
   ${REMOVE} batchscript.m

   set MatlabResult = `tail -1 $fname`

   if ($MatlabResult != 0) then
      set MatlabExists = -1
      echo "Matlab can not be run, see $fname"
   else
      ${REMOVE} $fname
   endif

endif

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
echo "Testing single-threaded lorenz_96 (L96) at "`date`
echo "=================================================================="
echo

set MODEL = lorenz_96

cd ${DARTHOME}/models/${MODEL}/work

# Save the 'original' files so we can reinstate them as needed

${COPY} input.nml     input.nml.$$
${COPY} obs_seq.in   obs_seq.in.$$
${COPY} obs_seq.out obs_seq.out.$$
${COPY} perfect_ics perfect_ics.$$
${COPY} filter_ics   filter_ics.$$

# The results of each run-time test will be stored in this directory.

if ( ! -d ${BASEOUTPUTDIR} ) then
   mkdir -p ${BASEOUTPUTDIR}
endif

foreach TEST ( spun_up 10hour 10hour.mres baseline async2 \
               Fasync0 Fasync0mres \
               prior_osi_1 prior_osi_2 prior_osi_3 )
   echo -n "Cleaning out ${BASEOUTPUTDIR}/${TEST} ... "
   ${REMOVE} ${BASEOUTPUTDIR}/${TEST}
   echo "done."
end

# Begin by compiling all programs; need to stop if an error is detected
./quickbuild.csh -nompi || exit 29

${REMOVE} *.o *.mod 

echo
echo
echo
echo "=================================================================="
echo "NewExperiment: spun_up"
echo "Setup appropriate namelists for a 4x20 1000-step test run."
echo "=================================================================="
echo

set EXP = ${BASEOUTPUTDIR}/spun_up
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous

echo "input.nml initially contains:"
cat input.nml
echo
echo "END of initial input.nml"
echo

${REMOVE} obs_seq.in obs_seq.out obs_seq.final temp_input vi_script

# Create an obs_sequence file
./create_obs_sequence < ../random_obs.input || exit 30
echo 'set_def.out'	       >! temp_input
echo '1'                       >> temp_input
echo '1000'                    >> temp_input
echo '0 0'                     >> temp_input
echo '0 3600'                  >> temp_input
echo 'obs_seq.in'              >> temp_input

echo "create_fixed_network_seq input is "
cat temp_input
echo

./create_fixed_network_seq < temp_input      || exit 31

# Need to modify rest of input.nml for test run
echo ':0'                             >! vi_script
echo '/ens_size'                      >> vi_script
echo ':s/20/80/'                      >> vi_script
echo '/num_groups'                    >> vi_script
echo ':s/1/4/'                        >> vi_script
echo '/save_reg_diagnostics'          >> vi_script
echo ':s/false/true/'                 >> vi_script
echo ':g/silence/s/false/true/'       >> vi_script
echo ':wq'                            >> vi_script

(${VI} input.nml < vi_script)

echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

# Run the perfect model and the filter
./perfect_model_obs  || exit 32
./filter             || exit 33

# Establish baseline for the run.
if ( $MatlabExists == 0 ) then
   ${REMOVE} batchscript.m
   echo "plot_total_err"             >! batchscript.m
   echo "print -dpsc DART_fig1.ps"   >> batchscript.m
   echo "quit"                       >> batchscript.m
   matlab -nosplash -nojvm -r batchscript
   ${REMOVE} batchscript.m
   ${MOVE} DART_fig1.ps ${EXP}
endif

# These are needed for comparison to other configurations.

${COPY} perfect_restart       perfect_ics.spun_up
${COPY}  filter_restart        filter_ics.spun_up

# Now, just save 'everything'
# Even the stuff from the spinup

${MOVE} obs_seq.in                 ${EXP}

${COPY} input.nml                  ${EXP}
${MOVE} obs_seq.out                ${EXP}
${MOVE} True_State.nc              ${EXP}
${MOVE} perfect_restart            ${EXP}

${MOVE} Prior_Diag.nc              ${EXP}
${MOVE} Posterior_Diag.nc          ${EXP}
${MOVE} obs_seq.final              ${EXP}
${MOVE} filter_restart             ${EXP}
${MOVE} reg_diagnostics            ${EXP}
${MOVE} configuration_being_tested ${EXP}
${MOVE} dart_log.out               ${EXP}

# General cleanup

${REMOVE} temp_input vi_script

ls -lrt | tail -30

echo
echo
echo
echo "=================================================================="
echo "NewExperiment: 10hour"
echo "Prepare to set up a sequence of 10 hour runs for testing: '10hour'"
echo "=================================================================="
echo "Just do 10 hours and output filter restarts in both the single and"
echo "multiple file format for later testing.  To reproduce across a"
echo "variety of options, need binary restart files."
echo

# Create an obs_sequence file for 10 hour tests

./create_obs_sequence < ../random_obs.input  || exit 40

echo 'set_def.out'             >! temp_input
echo '1'                       >> temp_input
echo '10'                      >> temp_input
echo '0 0'                     >> temp_input
echo '0 3600'                  >> temp_input
echo 'obs_seq.in'              >> temp_input
echo ':wq'                     >> temp_input

echo "create_fixed_network_seq input is "
cat temp_input
echo " "

./create_fixed_network_seq < temp_input      || exit 41

set EXP = ${BASEOUTPUTDIR}/10hour
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous

echo ':0'                                 >! vi_script
echo '/restart_in_file_name'              >> vi_script
echo ':s/perfect_ics/perfect_ics.spun_up' >> vi_script
echo '/filter_nml'                        >> vi_script
echo '/restart_in_file_name'              >> vi_script
echo ':s/filter_ics/filter_ics.spun_up'   >> vi_script
echo '/write_binary_restart_files'        >> vi_script
echo ':s/.false./.true./'                 >> vi_script
echo ':wq'                                >> vi_script

(${VI} input.nml < vi_script)

echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo
echo "input.nml now contains:"
cat input.nml
echo
echo "END of input.nml"

./perfect_model_obs   || exit 42
./filter              || exit 43

# These are needed for comparison to other configurations.

${COPY} perfect_restart       perfect_ics.10hour
${COPY}  filter_restart        filter_ics.10hour

# Now, just save 'everything'

${COPY} input.nml                  ${EXP}
${COPY} obs_seq.out                ${EXP}
${MOVE} True_State.nc              ${EXP}
${MOVE} perfect_restart            ${EXP}

${MOVE} Prior_Diag.nc              ${EXP}
${MOVE} Posterior_Diag.nc          ${EXP}
${MOVE} obs_seq.final              ${EXP}
${MOVE} filter_restart             ${EXP}
${MOVE} reg_diagnostics            ${EXP}
${MOVE} configuration_being_tested ${EXP}
${MOVE} dart_log.out               ${EXP}

${REMOVE} vi_script temp_input

ls -lrt | tail -30

echo
echo
echo
echo "=================================================================="
echo "NewExperiment: 10hour.mres"
echo "Now run the filter again to produce multiple restart files:"
echo "=================================================================="
echo

set EXP = ${BASEOUTPUTDIR}/10hour.mres
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested
set MULTIRESTART = $EXP

${COPY} input.nml input.nml.previous

echo ':0'                                >! vi_script
echo '/single_restart_file_out'          >> vi_script
echo ':s/true/false/'                    >> vi_script
echo ':wq'                               >> vi_script

(${VI} input.nml < vi_script)

echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter  || exit 53

${COPY} input.nml                  ${EXP} 
${MOVE} Prior_Diag.nc              ${EXP}
${MOVE} Posterior_Diag.nc          ${EXP}
${MOVE} obs_seq.final              ${EXP}
${MOVE} filter_restart.*           ${EXP}
${MOVE} reg_diagnostics            ${EXP}
${MOVE} configuration_being_tested ${EXP}
${MOVE} dart_log.out               ${EXP}

${REMOVE} vi_script

ls -lrt | tail -30

echo
echo
echo
echo "=================================================================="
echo "NewExperiment: baseline"
echo "Now do a second 10 hour run from the end of the first 10 hour"
echo "Also change the filter back to only produce single restart file"
echo "This is the configuration that sets the standard."
echo "=================================================================="
echo

set EXP = ${BASEOUTPUTDIR}/baseline
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
set BASELINE = ${EXP}
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous

echo ':0'                                         >! vi_script
echo '/restart_in_file_name'                      >> vi_script
echo ':s/perfect_ics.spun_up/perfect_ics.10hour'  >> vi_script
echo '/filter_nml'                                >> vi_script
echo '/restart_in_file_name'                      >> vi_script
echo ':s/filter_ics.spun_up/filter_ics.10hour'    >> vi_script
echo '/single_restart_file_out'                   >> vi_script
echo ':s/false/true/'                             >> vi_script
echo ':wq'                                        >> vi_script

(${VI} input.nml < vi_script)

echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./perfect_model_obs  || exit 62
./filter             || exit 63

${COPY} input.nml                  ${EXP}
${MOVE} obs_seq.out                ${EXP}
${MOVE} True_State.nc              ${EXP}
${MOVE} perfect_restart            ${EXP}

${MOVE} Prior_Diag.nc              ${EXP}
${MOVE} Posterior_Diag.nc          ${EXP}
${MOVE} obs_seq.final              ${EXP}
${MOVE} filter_restart             ${EXP}
${MOVE} reg_diagnostics            ${EXP}
${MOVE} configuration_being_tested ${EXP}
${MOVE} dart_log.out               ${EXP}

${REMOVE} vi_script

ls -lrt | tail -30

echo
echo
echo
echo "=================================================================="
echo "NewExperiment: async2"
echo "Test two async options"
echo "async == 0 == model advance by subroutine (already done)"
echo "async == 2 == model advance with F90 calls to shell script"
echo "=================================================================="
echo

${COPY} ../shell_scripts/advance_model.csh .

set EXP = ${BASEOUTPUTDIR}/async2
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous

echo ':0'                            >! vi_script
echo '/perfect_model_obs_nml'        >> vi_script
echo '/async'                        >> vi_script
echo ':s/0/2/'                       >> vi_script
echo ':wq'                           >> vi_script

(${VI} input.nml < vi_script)

echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./perfect_model_obs  || exit 82

diff perfect_restart  ${BASELINE}/perfect_restart || exit 84
diff     obs_seq.out  ${BASELINE}/obs_seq.out     || exit 85
#diff  True_State.nc  ${BASELINE}/True_State.nc   || exit 86

${COPY} input.nml                  ${EXP}
${COPY} advance_model.csh          ${EXP}
${MOVE} obs_seq.out                ${EXP}
${MOVE} True_State.nc              ${EXP}
${MOVE} perfect_restart            ${EXP}
${MOVE} configuration_being_tested ${EXP}
${MOVE} dart_log.out               ${EXP}
${MOVE} assim_model_state_ud?????  ${EXP}

${REMOVE} vi_script

ls -lrt | tail -30

echo
echo
echo
echo "=================================================================="
echo "NewExperiment: Fasync0"
echo "Next, start checking filter options for this case"
echo "Begin by checking single versus multiple restarts"
echo "=================================================================="
echo

set EXP = ${BASEOUTPUTDIR}/Fasync0
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} ${BASELINE}/obs_seq.out    .

./filter  || exit 73

diff      obs_seq.final  ${BASELINE}/obs_seq.final     || exit 10
diff     filter_restart  ${BASELINE}/filter_restart    || exit 11
diff    reg_diagnostics  ${BASELINE}/reg_diagnostics   || exit 12

# Determine if the run was correct.
# can fire off a matlab batch job if matlab exists.
if ( $MatlabExists == 0 ) then

   ${COPY} ${BASELINE}/True_State.nc  .

   echo "plot_total_err"             >! batchscript.m
   echo "print -dpsc DART_fig2.ps"   >> batchscript.m
   echo "quit"                       >> batchscript.m
   matlab -nosplash -nojvm -r batchscript
   ${REMOVE} batchscript.m
   ${REMOVE} True_State.nc
   ${MOVE} DART_fig2.ps ${EXP}
endif

${COPY} input.nml                  ${EXP}
${MOVE} Prior_Diag.nc              ${EXP}
${MOVE} Posterior_Diag.nc          ${EXP}
${MOVE} obs_seq.final              ${EXP}
${MOVE} filter_restart             ${EXP}
${MOVE} reg_diagnostics            ${EXP}
${MOVE} configuration_being_tested ${EXP}
${MOVE} dart_log.out               ${EXP}

ls -lrt | tail -30

echo
echo
echo
echo "=================================================================="
echo "NewExperiment: Fasync0mres"
echo "Now use multiple input files (named filter_restart.00xx)"
echo "=================================================================="
echo

set EXP = ${BASEOUTPUTDIR}/Fasync0mres
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous

${COPY} ${MULTIRESTART}/filter_restart.*  .

echo ':0'                                    >! vi_script
echo '/filter_nml'                           >> vi_script
echo '/restart_in_file_name'                 >> vi_script
echo ':s/filter_ics.10hour/filter_restart/'  >> vi_script
echo '/single_restart_file_in'               >> vi_script
echo ':s/true/false/'                        >> vi_script
echo ':wq'                                   >> vi_script

(${VI} input.nml < vi_script)

echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter  || exit 83

diff      obs_seq.final  ${BASELINE}/obs_seq.final     || exit 10
diff     filter_restart  ${BASELINE}/filter_restart    || exit 11
diff    reg_diagnostics  ${BASELINE}/reg_diagnostics   || exit 12
#diff     Prior_Diag.nc  ${BASELINE}/Prior_Diag.nc
#diff Posterior_Diag.nc  ${BASELINE}/Posterior_Diag.nc

${COPY} input.nml                  ${EXP}
${MOVE} Prior_Diag.nc              ${EXP}
${MOVE} Posterior_Diag.nc          ${EXP}
${MOVE} obs_seq.final              ${EXP}
${MOVE} filter_restart             ${EXP}
${MOVE} reg_diagnostics            ${EXP}
${MOVE} configuration_being_tested ${EXP}
${MOVE} dart_log.out               ${EXP}

${REMOVE} vi_script

ls -lrt | tail -30

echo
echo
echo
echo "=================================================================="
echo "NewExperiment: prior_osi_1"
echo "Test the observation space inflation option"
echo "=================================================================="
echo

set EXP = ${BASEOUTPUTDIR}/prior_osi_1
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous

echo ':0'                                     >! vi_script
echo '/ens_size'                              >> vi_script
echo ':s/80/20/'                              >> vi_script
echo '/num_groups'                            >> vi_script
echo ':s/4/1/'                                >> vi_script
echo '/inf_flavor'                            >> vi_script
echo ':s/2/1/'                                >> vi_script
echo '/inf_initial'                           >> vi_script
echo ':s/1\.0/1\.05/'                         >> vi_script
echo ':wq'                                    >> vi_script

(${VI} input.nml < vi_script)

echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter  || exit 123

${COPY} input.nml                    ${EXP}
${MOVE} Prior_Diag.nc                ${EXP}
${MOVE} Posterior_Diag.nc            ${EXP}
${MOVE} obs_seq.final                ${EXP}
${MOVE} filter_restart               ${EXP}
${MOVE} configuration_being_tested   ${EXP}
${MOVE} prior_inflate_restart        ${EXP}
${MOVE} prior_inflate_diag           ${EXP}
${MOVE} dart_log.out                 ${EXP}

${REMOVE} vi_script

ls -lrt | tail -30

echo
echo
echo
echo "=================================================================="
echo "NewExperiment: prior_osi_2"
echo "Test the spatially varying state inflation option"
echo "=================================================================="
echo

set EXP = ${BASEOUTPUTDIR}/prior_osi_2
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous

echo ':0'                                     >! vi_script
echo '/inf_flavor'                            >> vi_script
echo ':s/1/2/'                                >> vi_script
echo ':wq'                                    >> vi_script

(${VI} input.nml < vi_script)

echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter  || exit 123

${COPY} input.nml                    ${EXP}
${MOVE} Prior_Diag.nc                ${EXP}
${MOVE} Posterior_Diag.nc            ${EXP}
${MOVE} obs_seq.final                ${EXP}
${MOVE} filter_restart               ${EXP}
${MOVE} configuration_being_tested   ${EXP}
${MOVE} prior_inflate_restart        ${EXP}
${MOVE} prior_inflate_diag           ${EXP}
${MOVE} dart_log.out                 ${EXP}

${REMOVE} vi_script

ls -lrt | tail -30

echo
echo
echo
echo "=================================================================="
echo "NewExperiment: prior_osi_3"
echo "Test the fixed state space inflation option"
echo "=================================================================="
echo

set EXP = ${BASEOUTPUTDIR}/prior_osi_3
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous

echo ':0'                                     >! vi_script
echo '/inf_flavor'                            >> vi_script
echo ':s/2/3/'                                >> vi_script
echo ':wq'                                    >> vi_script

(${VI} input.nml < vi_script)

echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter  || exit 123

${COPY} input.nml                    ${EXP}
${MOVE} Prior_Diag.nc                ${EXP}
${MOVE} Posterior_Diag.nc            ${EXP}
${MOVE} obs_seq.final                ${EXP}
${MOVE} filter_restart               ${EXP}
${MOVE} configuration_being_tested   ${EXP}
${MOVE} prior_inflate_restart        ${EXP}
${MOVE} prior_inflate_diag           ${EXP}
${MOVE} dart_log.out                 ${EXP}

${REMOVE} vi_script

ls -lrt | tail -30

echo
echo
echo
echo "=================================================================="
echo "Single-threaded testing of lorenz_96 complete at "`date`
echo "=================================================================="
echo

echo "=================================================================="
echo "cleaning up ... and restoring original input.nml, ics ... "
${REMOVE} filter_restart.*
${REMOVE} input.nml perfect_ics perfect_ics.spun_up perfect_ics.10hour
${REMOVE} obs_seq.in filter_ics  filter_ics.spun_up  filter_ics.10hour
${REMOVE} obs_seq.out set_def.out  True_State.nc input.nml.previous
${REMOVE} *.o *.mod 
${REMOVE} Makefile input.nml.*_default .cppdefs
foreach TARGET ( mkmf_* )
  set PROG = `echo $TARGET | sed -e 's#mkmf_##'`
  ${REMOVE} $PROG
end

# Reinstate the 'original' files so we can run this again if we need to.

${MOVE}   input.nml.$$   input.nml
${MOVE}  obs_seq.in.$$ obs_seq.in
${MOVE} obs_seq.out.$$ obs_seq.out
${MOVE} perfect_ics.$$ perfect_ics
${MOVE}  filter_ics.$$  filter_ics

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

