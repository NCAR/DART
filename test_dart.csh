#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$

set SNAME = $0
set clobber

switch ( $#argv )
   case 1:
      # supplying one argument -- the base output directory.
      # The results of each run-time test will be stored in this directory.
      # This will facilitate checks across platforms.

      breaksw
   default:
      echo " "
      echo "usage: $SNAME:t OutputDirectory"
      echo " "
      echo "This script compiles several program units for a wide range of models and then"
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
      setenv REMOVE 'rm -rvf'
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
      echo "the test input files."
      exit 2
   endif
endif

echo "Running DART test on $host"

#----------------------------------------------------------------------
# Compile 'filter' for a wide range of models.
#----------------------------------------------------------------------

@ modelnum = 31

if ( 1 == 1 ) then
foreach MODEL ( 9var lorenz_63 lorenz_84 lorenz_96 lorenz_96_2scale \
    lorenz_04 forced_lorenz_96 bgrid_solo cam wrf pe2lyr ) 
    # null_model PBL_1d rose MITgcm_annulus sccm )

    echo "=================================================================="
    echo "Compiling $MODEL at "`date`
    echo ""

    cd ${DARTHOME}/models/${MODEL}/work

    ${REMOVE} ../../../obs_def/obs_def_mod.f90 
    ${REMOVE} ../../../obs_kind/obs_kind_mod.f90
    ${REMOVE} preprocess
    ${REMOVE} *.o *.mod

    @ makenum  = 1

    csh mkmf_preprocess
    make || exit 1
    ./preprocess

    foreach PROG ( create_obs_sequence create_fixed_network_seq \
                   perfect_model_obs filter )

       ${REMOVE} ${PROG} Makefile input.nml.${PROG}_default .cppdefs

       csh mkmf_${PROG}  || exit $modelnum
       make              || exit $makenum

       @ makenum  = $makenum  + 1

       ${REMOVE} ${PROG} Makefile input.nml.${PROG}_default .cppdefs

    end

    ${REMOVE} *.o *.mod

   @ modelnum = $modelnum + 1
end
endif

#----------------------------------------------------------------------
# Lots of tests for L96
#----------------------------------------------------------------------

cd ${DARTHOME}

# If matlab and the netcdf functions exist, generate figures on-the-fly.

which matlab > /dev/null
set MatlabExists=$status

if ( $MatlabExists == 0 ) then

   set fname = dart.matlabcheck.$$

   echo "fname = '"$fname"';"        >! batchscript.m
   echo "addpath ${DARTHOME}/matlab" >> batchscript.m 
   echo "ChecknetCDFuse(fname);"     >> batchscript.m
   echo "quit"                       >> batchscript.m
   matlab -nosplash -nojvm -r batchscript
   ${REMOVE} batchscript.m

   set MatlabResult = `tail -1 $fname`

   if ($MatlabResult != 0) then
      set MatlabExists = -1
      echo "Matlab can not be run, see $fname"
   endif

endif

echo ""
echo "=================================================================="
echo "Testing lorenz_96 (L96) at "`date`
echo "=================================================================="
echo ""

cd ${DARTHOME}/models/lorenz_96/work

# Save the 'original' files so we can reinstate them later

${COPY} input.nml     input.nml.$$
${COPY} obs_seq.in   obs_seq.in.$$
${COPY} obs_seq.out obs_seq.out.$$
${COPY} perfect_ics perfect_ics.$$
${COPY} filter_ics   filter_ics.$$

# The results of each run-time test will be stored in this directory.

set BASEOUTPUTDIR = $1
if ( ! -d ${BASEOUTPUTDIR} ) then
   mkdir -p ${BASEOUTPUTDIR}
endif

foreach TEST ( spun_up 10hour 10hour.mres baseline out_of_core async2 \
               async3 Fasync0 async0mres3dom async0mres5dom \
               async2mres3dom async3mres3dom async0mres0prll \
               osi svsi sssi )
   echo -n "Cleaning out ${BASEOUTPUTDIR}/${TEST} ... "
   ${REMOVE} ${BASEOUTPUTDIR}/${TEST}
   echo "done."
end

# Make sure that all .o, .mod and executables are gone
${REMOVE} *.o *.mod assim_region create_fixed_network_seq create_obs_seq filter
${REMOVE} integrate_model perfect_model_obs ../../../obs_kind/obs_kind_mod.f90
${REMOVE} merge_obs_seq obs_diag smoother ../../../obs_def/obs_def_mod.f90

# Begin by compiling all programs; need to stop if an error is detected
csh mkmf_preprocess               || exit   1
make                              || exit   2
./preprocess                      || exit   3
# Preprocess must be done before rest of makes
csh mkmf_assim_region             || exit   4
make                              || exit   5
csh mkmf_create_fixed_network_seq || exit   6
make                              || exit   7
csh mkmf_create_obs_sequence      || exit   8
make                              || exit   9
csh mkmf_filter                   || exit  10
make                              || exit  11
csh mkmf_integrate_model          || exit  12
make                              || exit  13
csh mkmf_perfect_model_obs        || exit  14
make                              || exit  15
csh mkmf_merge_obs_seq            || exit  16
make                              || exit  17
csh mkmf_obs_diag                 || exit  18
make                              || exit  19
csh mkmf_smoother                 || exit  20
make                              || exit  21

echo ""
echo "=================================================================="
echo "Setup appropriate namelists for a 4x20 1000-step test run."
echo "=================================================================="
echo ""

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
echo " "
echo "END of initial input.nml"
echo " "

# Create an obs_sequence file
${REMOVE} obs_seq.in obs_seq.out obs_seq.final temp_input vi_script
./create_obs_sequence < ../random_obs.input || exit 30
echo 'set_def.out'	       >! temp_input
echo '1'                       >> temp_input
echo '1000'                    >> temp_input
echo '0 0'                     >> temp_input
echo '0 3600'                  >> temp_input
echo 'obs_seq.in'              >> temp_input

echo "create_fixed_network_seq input is "
cat temp_input
echo " "

./create_fixed_network_seq < temp_input      || exit 31

# Need to modify rest of input.nml for test run
echo ':0'                             >! vi_script
echo '/ens_size'                      >> vi_script
echo ':s/20/80/'                      >> vi_script
echo '/num_groups'                    >> vi_script
echo ':s/1/4/'                        >> vi_script
echo ':wq'                            >> vi_script
(${VI} input.nml < vi_script ) || exit 98
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

# Even the stuff from the spinup

${MOVE} obs_seq.in  ${EXP}

# General cleanup

${REMOVE} temp_input vi_script go_end_filter

ls -lrt | tail -30    # TJH debug

echo ""
echo "------------------------------------------------------------------"
echo "Prepare to set up a sequence of 10 hour runs for testing"
echo "------------------------------------------------------------------"
echo "In first case, just do 10 hours and output filter restarts in both"
echo "the single file and multiple file format for later testing"
echo "To reproduce across a variety of options, need num_domains >= 2"
echo "and binary restart files to be written."
echo "=================================================================="
echo ""

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

echo " "
echo "=================================================================="
echo "To reproduce across a variety of options, need num_domains 2 or greater"
echo "and binary restart files. Start from previous restart and create "
echo "single  previous restart and create new restart."
echo "=================================================================="
echo " "

set EXP = ${BASEOUTPUTDIR}/10hour
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous
${REMOVE} vi_script

echo ':0'                                 >! vi_script
echo '/num_domains'                       >> vi_script
echo ':s/1/3/'                            >> vi_script
echo ':0'                                 >> vi_script
echo '/restart_in_file_name'              >> vi_script
echo ':s/perfect_ics/perfect_ics.spun_up' >> vi_script
echo '/filter_nml'                        >> vi_script
echo '/restart_in_file_name'              >> vi_script
echo ':s/filter_ics/filter_ics.spun_up'   >> vi_script
echo '/write_binary_restart_files'        >> vi_script
echo ':s/.false./.true./'                 >> vi_script
echo ':wq'                                >> vi_script
(${VI} input.nml < vi_script ) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

echo ""
echo "=================================================================="
echo "Run the perfect model and the filter to produce single restart files:"
echo "perfect_ics.10hour and filter_ics.10hour"
echo "=================================================================="
echo ""

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

${REMOVE} vi_script temp_input go_end_filter

ls -lrt | tail -30    # TJH debug

echo ""
echo "=================================================================="
echo "Now run the filter again to produce multiple restart files:"
echo "filter_restart.01"
echo "=================================================================="
echo ""
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
(${VI} input.nml < vi_script ) || exit 98
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

${REMOVE} vi_script go_end_filter

ls -lrt | tail -30    # TJH debug

echo " "
echo "=================================================================="
echo "Now do a second 10 hour run from the end of the first 10 hour"
echo "Also change the filter back to only produce single restart file"
echo "This is the configuration that sets the standard."
echo "=================================================================="
echo " "
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
(${VI} input.nml < vi_script ) || exit 98
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

${REMOVE} vi_script go_end_filter

ls -lrt | tail -30    # TJH debug

echo ""
echo "=================================================================="
echo "Test storing the ensemble on disk instead of in core"
echo "=================================================================="
echo ""
set EXP = ${BASEOUTPUTDIR}/out_of_core
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous

echo ':0'                             >! vi_script
echo '/ensemble_manager'              >> vi_script
echo '/in_core'                       >> vi_script
echo ':s/true/false/'                 >> vi_script
echo ':wq'                            >> vi_script
(${VI} input.nml < vi_script ) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./perfect_model_obs  || exit 72

 diff perfect_restart             ${BASELINE}/perfect_restart  || exit 74
 diff     obs_seq.out             ${BASELINE}/obs_seq.out      || exit 75
#diff  True_State.nc.out_of_core  ${BASELINE}/True_State.nc    || exit 76

${COPY} input.nml                      ${EXP}
${MOVE} obs_seq.out                    ${EXP}
${MOVE} True_State.nc                  ${EXP}
${MOVE} perfect_restart                ${EXP}
${MOVE} ens_manager_ens_file.0001.0001 ${EXP}
${MOVE} configuration_being_tested     ${EXP}
${MOVE} dart_log.out                   ${EXP}

${REMOVE} vi_script go_end_filter

ls -lrt | tail -30    # TJH debug

#-----------------------------------------------------------------------
# Need to get the scripts 
#-----------------------------------------------------------------------
${COPY} ../shell_scripts/*.csh .
${COPY} ${DARTHOME}/shell_scripts/advance_ens.csh      advance_ens.csh
${COPY} ${DARTHOME}/shell_scripts/assim_filter.csh    assim_filter.csh
${COPY} ${DARTHOME}/shell_scripts/filter_server.csh  filter_server.csh

echo ""
echo "=================================================================="
echo "Test two async options"
echo "async == 0 == model advance by subroutine"
echo "async == 2 == model advance with F90 calls to shell script"
echo "async == 3 == model advance when signalled by semaphore file"
echo "Change to async=2 and in_core=true for perfect model."
echo "=================================================================="
echo ""
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
echo '/in_core'                      >> vi_script
echo ':s/false/true/'                >> vi_script
echo ':wq'                           >> vi_script
(${VI} input.nml < vi_script ) || exit 98
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
${MOVE} obs_seq.out                ${EXP}
${MOVE} True_State.nc              ${EXP}
${MOVE} perfect_restart            ${EXP}
${MOVE} assim_model_state_*        ${EXP}
${MOVE} filter_control             ${EXP}
${MOVE} integrate_model_out_temp1  ${EXP}
${MOVE} configuration_being_tested ${EXP}
${MOVE} member_1                   ${EXP}
${MOVE} dart_log.out               ${EXP}

${REMOVE} vi_script go_end_filter

ls -lrt | tail -30    # TJH debug

echo ""
echo "=================================================================="
echo "Now try async == 3 (filter_server.csh driving the model advance)"
echo "=================================================================="
echo ""
set EXP = ${BASEOUTPUTDIR}/async3
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous
${REMOVE} vi_script

echo ':0'                            >! vi_script
echo '/perfect_model_obs_nml'        >> vi_script
echo '/async'                        >> vi_script
echo ':s/2/3/'                       >> vi_script
echo ':wq'                           >> vi_script
(${VI} input.nml < vi_script ) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter_server.csh &
./perfect_model_obs || exit 92

echo "If you're still here after 10 seconds or so ..."
echo "filter_server is not terminating when it should."
wait

diff perfect_restart  ${BASELINE}/perfect_restart || exit 94
diff     obs_seq.out  ${BASELINE}/obs_seq.out     || exit 95
#diff  True_State.nc  ${BASELINE}/True_State.nc   || exit 96

${COPY} input.nml                  ${EXP}
${MOVE} obs_seq.out                ${EXP}
${MOVE} True_State.nc              ${EXP}
${MOVE} perfect_restart            ${EXP}
${MOVE} assim_model_state_*        ${EXP}
${MOVE} filter_control             ${EXP}
${MOVE} integrate_model_out_temp1  ${EXP}
${MOVE} configuration_being_tested ${EXP}
${MOVE} run_job.log                ${EXP}
${MOVE} dart_log.out               ${EXP}

${REMOVE} vi_script go_end_filter

ls -lrt | tail -30    # TJH debug

echo ""
echo "=================================================================="
echo "Next, start checking filter options for this case"
echo "Begin by checking single versus multiple restarts"
echo "=================================================================="
echo ""

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

${REMOVE} vi_script go_end_filter

ls -lrt | tail -30    # TJH debug

echo ""
echo "=================================================================="
echo "Now use multiple input files (named filter_restart.00xx)"
echo "=================================================================="
echo ""
set EXP = ${BASEOUTPUTDIR}/async0mres3dom
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
(${VI} input.nml < vi_script ) || exit 98
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

${REMOVE} vi_script go_end_filter

ls -lrt | tail -30    # TJH debug

echo ""
echo "=================================================================="
echo "Next switch the number of domains from 3 to 5"
echo "=================================================================="
echo ""
set EXP = ${BASEOUTPUTDIR}/async0mres5dom
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous

echo ':0'                                     >! vi_script
echo '/num_domains'                           >> vi_script
echo ':s/3/5/'                                >> vi_script
echo ':wq'                                    >> vi_script
(${VI} input.nml < vi_script ) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter  || exit 93

diff      obs_seq.final  ${BASELINE}/obs_seq.final     || exit 10
diff     filter_restart  ${BASELINE}/filter_restart    || exit 11

${COPY} input.nml                  ${EXP}
${MOVE} Prior_Diag.nc              ${EXP}
${MOVE} Posterior_Diag.nc          ${EXP}
${MOVE} obs_seq.final              ${EXP}
${MOVE} filter_restart             ${EXP}
${MOVE} reg_diagnostics            ${EXP}
${MOVE} configuration_being_tested ${EXP}
${MOVE} dart_log.out               ${EXP}

${REMOVE} vi_script go_end_filter

ls -lrt | tail -30    # TJH debug

echo ""
echo "=================================================================="
echo "Next switch the number of domains back to 3; try parallel option 2"
echo "=================================================================="
echo ""
set EXP = ${BASEOUTPUTDIR}/async2mres3dom
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous

echo ':0'                                     >! vi_script
echo '/do_parallel'                           >> vi_script
echo ':s/0/2/'                                >> vi_script
echo '/num_domains'                           >> vi_script
echo ':s/5/3/'                                >> vi_script
echo ':wq'                                    >> vi_script
(${VI} input.nml < vi_script ) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter || exit 103

diff      obs_seq.final  ${BASELINE}/obs_seq.final     || exit 10
diff     filter_restart  ${BASELINE}/filter_restart    || exit 11

${COPY} input.nml                    ${EXP}
${MOVE} Prior_Diag.nc                ${EXP}
${MOVE} Posterior_Diag.nc            ${EXP}
${MOVE} obs_seq.final                ${EXP}
${MOVE} filter_restart               ${EXP}
${MOVE} configuration_being_tested   ${EXP}
${MOVE} filter_assim_region__in[123] ${EXP}
${MOVE} filter_assim_region_out[123] ${EXP}
${MOVE}   assim_region_out_temp[123] ${EXP}
${MOVE} filter_assim_obs_seq         ${EXP}
${MOVE} assim_region_control         ${EXP}
${MOVE} region_[123]                 ${EXP}
${MOVE} dart_log.out                 ${EXP}

${REMOVE} vi_script go_end_filter

ls -lrt | tail -30    # TJH debug

echo ""
echo "=================================================================="
echo "Try parallel option 3"
echo "=================================================================="
echo ""
set EXP = ${BASEOUTPUTDIR}/async3mres3dom
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous

echo ':0'                                     >! vi_script
echo '/do_parallel'                           >> vi_script
echo ':s/2/3/'                                >> vi_script
echo ':wq'                                    >> vi_script
(${VI} input.nml < vi_script ) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter_server.csh &
./filter  || exit 113

echo "If you're still here after 10 seconds or so ..."
echo "filter_server is not terminating when it should."
wait

diff      obs_seq.final  ${BASELINE}/obs_seq.final     || exit 10
diff     filter_restart  ${BASELINE}/filter_restart    || exit 11

${COPY} input.nml                    ${EXP}
${MOVE} Prior_Diag.nc                ${EXP}
${MOVE} Posterior_Diag.nc            ${EXP}
${MOVE} obs_seq.final                ${EXP}
${MOVE} filter_restart               ${EXP}
${MOVE} configuration_being_tested   ${EXP}
${MOVE} filter_assim_region__in[123] ${EXP}
${MOVE} filter_assim_region_out[123] ${EXP}
${MOVE}   assim_region_out_temp[123] ${EXP}
${MOVE} filter_assim_obs_seq         ${EXP}
${MOVE} assim_region_control         ${EXP}
${MOVE} run_job.log                  ${EXP}
${MOVE} dart_log.out                 ${EXP}

${REMOVE} vi_script go_end_filter

ls -lrt | tail -30    # TJH debug

echo ""
echo "=================================================================="
echo "Go back to parallel option 0 and proceed"
echo "=================================================================="
echo ""
set EXP = ${BASEOUTPUTDIR}/async0mres0prll
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous

echo ':0'                                     >! vi_script
echo '/do_parallel'                           >> vi_script
echo ':s/3/0/'                                >> vi_script
echo ':wq'                                    >> vi_script
(${VI} input.nml < vi_script ) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter  || exit 123

diff      obs_seq.final  ${BASELINE}/obs_seq.final     || exit 10
diff     filter_restart  ${BASELINE}/filter_restart    || exit 11
diff    reg_diagnostics  ${BASELINE}/reg_diagnostics   || exit 12

${COPY} input.nml                    ${EXP}
${MOVE} Prior_Diag.nc                ${EXP}
${MOVE} Posterior_Diag.nc            ${EXP}
${MOVE} obs_seq.final                ${EXP}
${MOVE} filter_restart               ${EXP}
${MOVE} configuration_being_tested   ${EXP}
${MOVE} reg_diagnostics              ${EXP}
${MOVE} dart_log.out                 ${EXP}

${REMOVE} vi_script go_end_filter

ls -lrt | tail -30    # TJH debug

echo ""
echo "=================================================================="
echo "Test the observation space inflation option"
echo "=================================================================="
echo ""
set EXP = ${BASEOUTPUTDIR}/osi
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
echo '/cutoff'                                >> vi_script
echo ':s/1000000\.0/0\.2/'                    >> vi_script
echo '/num_domains'                           >> vi_script
echo ':s/3/1/'                                >> vi_script
echo '/do_obs_inflate'                        >> vi_script
echo ':s/false/true/'                         >> vi_script
echo '/obs_inf_sd_initial'                    >> vi_script
echo ':s/0\.0/0\.1/'                          >> vi_script
echo ':wq'                                    >> vi_script
(${VI} input.nml < vi_script ) || exit 98
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
${MOVE} inflate_restart              ${EXP}
${MOVE} inflate_diag                 ${EXP}
${MOVE} dart_log.out                 ${EXP}

${REMOVE} vi_script go_end_filter

ls -lrt | tail -30    # TJH debug

echo ""
echo "=================================================================="
echo "Test the spatially varying state inflation option"
echo "=================================================================="
echo ""
set EXP = ${BASEOUTPUTDIR}/svsi
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous

echo ':0'                                     >! vi_script
echo '/do_obs_inflate'                        >> vi_script
echo ':s/true/false/'                         >> vi_script
echo '/do_varying_ss_inflate'                 >> vi_script
echo ':s/false/true/'                         >> vi_script
echo '/ss_inf_sd_initial'                     >> vi_script
echo ':s/0\.0/0\.1/'                          >> vi_script
echo ':wq'                                    >> vi_script
(${VI} input.nml < vi_script ) || exit 98
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
${MOVE} inflate_restart              ${EXP}
${MOVE} dart_log.out                 ${EXP}

${REMOVE} vi_script go_end_filter

ls -lrt | tail -30    # TJH debug

echo ""
echo "=================================================================="
echo "Test the single state space inflation option"
echo "=================================================================="
echo ""
set EXP = ${BASEOUTPUTDIR}/sssi
if (   -d    ${EXP} ) then
   ${REMOVE} ${EXP}/*
else
   mkdir     ${EXP}
endif
echo $EXP >! configuration_being_tested

${COPY} input.nml input.nml.previous

echo ':0'                                     >! vi_script
echo '/do_varying_ss_inflate'                 >> vi_script
echo ':s/true/false/'                         >> vi_script
echo '/do_single_ss_inflate'                  >> vi_script
echo ':s/false/true/'                         >> vi_script
echo ':wq'                                    >> vi_script
(${VI} input.nml < vi_script ) || exit 98
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
${MOVE} inflate_restart              ${EXP}
${MOVE} inflate_diag                 ${EXP}
${MOVE} dart_log.out                 ${EXP}

${REMOVE} vi_script go_end_filter

ls -lrt | tail -30    # TJH debug

echo "=================================================================="
echo "cleaning up ... and restoring original input.nml, ics ... "

${REMOVE} filter_restart.*
${REMOVE} input.nml perfect_ics perfect_ics.spun_up perfect_ics.10hour
${REMOVE} obs_seq.in filter_ics  filter_ics.spun_up  filter_ics.10hour
${REMOVE} obs_seq.out set_def.out  True_State.nc input.nml.previous

# Reinstate the 'original' files so we can run this again if we need to.

${MOVE}   input.nml.$$   input.nml
${MOVE}  obs_seq.in.$$ obs_seq.in
${MOVE} obs_seq.out.$$ obs_seq.out
${MOVE} perfect_ics.$$ perfect_ics
${MOVE}  filter_ics.$$  filter_ics

echo "=================================================================="
echo ""
echo "Testing complete  at "`date`
echo "=================================================================="
echo ""
