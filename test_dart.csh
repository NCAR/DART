#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section, 
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$

set SNAME = $0
set clobber

switch ( $#argv )
   case 0:
      # supplying no arguments -- echo usage not
      breaksw
   default:
      echo " "
      echo "usage: $SNAME:t"
      echo " "
      echo "This script compiles 'filter' for a wide range of models and then does"
      echo "relatively extensive tests of the L96 programs with a variety of options."
      echo "The L96 tests are best run from a 'clean' starting point - i.e. one that"
      echo "is as close to the distribution state as possible. Picking up in the middle"
      echo "is not particularly easy. I always run the script end-to-end. TJH"
      echo " "
      echo "This must be run from the top-level 'DART' directory."
      echo " "
      echo "This is a pretty verbose process, so if you are logging the output,"
      echo "make sure you have plenty of space:"
      echo " "
      echo "./$SNAME:t |& tee DART_test.log"
      echo " "
      echo "can easily result in a 750 Kb log file"
      exit 1
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
   set host = `uname -n`
endif

setenv REMOVE 'rm -rfv'
setenv COPY   'cp -pv'
setenv MOVE   'mv -v'

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

       rm ${REMOVE} ${PROG} Makefile input.nml.${PROG}_default .cppdefs

    end

    rm -f *.o *.mod

   @ modelnum = $modelnum + 1
end
endif

#----------------------------------------------------------------------
# Lots of tests for L96
#----------------------------------------------------------------------

cd ${DARTHOME}

# If matlab and the netcdf functions exist, generate figures on-the-fly.
# If not, simply archive output to test directories for later.

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

if ( $MatlabExists != 0 ) then
   mkdir -p ${DARTHOME}/Test1  # directory for output
   mkdir -p ${DARTHOME}/Test2  # directory for output
endif

echo ""
echo "=================================================================="
echo "Testing lorenz_96 (L96) at "`date`
echo "=================================================================="
echo ""

cd ${DARTHOME}/models/lorenz_96/work

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
#echo '/assim_tools_nml'               >> vi_script
#echo '/cov_inflate'                   >> vi_script
#echo ':s/-1.0/1.05/'                  >> vi_script
#echo '/cov_inflate_sd'                >> vi_script
#echo ':s/-0.05/0.05'                  >> vi_script
#echo '/sd_lower_bound'                >> vi_script
#echo ':s/-0.05/0.05'                  >> vi_script
echo ':wq'                            >> vi_script
(vi -s vi_script -e input.nml > /dev/null) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

# Run the perfect model and the filter
./perfect_model_obs  || exit 32
./filter             || exit 33

ls -lrt | tail -30    # TJH debug

# Establish baseline for the run.
if ( $MatlabExists == 0 ) then
   ${REMOVE} batchscript.m
   echo "plot_total_err"              > batchscript.m
   echo "print -dpsc DART_fig1.ps"   >> batchscript.m
   echo "quit"                       >> batchscript.m
   matlab -nosplash -nojvm -r batchscript
   ${REMOVE} batchscript.m
else
   echo "Matlab does not exist in your PATH ... archiving output"
   ${COPY} input.nml           ${DARTHOME}/Test1
   ${COPY} True_State.nc       ${DARTHOME}/Test1
   ${COPY} obs_seq.out         ${DARTHOME}/Test1
   ${COPY} perfect_restart     ${DARTHOME}/Test1

   ${COPY} Prior_Diag.nc       ${DARTHOME}/Test1
   ${COPY} Posterior_Diag.nc   ${DARTHOME}/Test1
   ${COPY} obs_seq.final       ${DARTHOME}/Test1
   ${COPY} inflate_restart     ${DARTHOME}/Test1
   ${COPY} filter_restart      ${DARTHOME}/Test1
   ${COPY} reg_diagnostics     ${DARTHOME}/Test1
endif

${MOVE} perfect_restart       perfect_ics.spun_up
${MOVE} True_State.nc          True_State.spun_up.nc

${MOVE}     Prior_Diag.nc      Prior_Diag.spun_up.nc
${MOVE} Posterior_Diag.nc  Posterior_Diag.spun_up.nc
${MOVE}   obs_seq.final     obs_seq.final.spun_up
${MOVE} inflate_restart   inflate_restart.spun_up    # zero size
${MOVE}  filter_restart        filter_ics.spun_up
${MOVE} reg_diagnostics   reg_diagnostics.spun_up

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
${REMOVE} obs_seq.in obs_seq.out obs_seq.final temp_input
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
(vi -s vi_script -e input.nml > /dev/null) || exit 98
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

ls -lrt | tail -30    # TJH debug

${COPY}     obs_seq.out       obs_seq.out.10hour
${MOVE} perfect_restart       perfect_ics.10hour
${MOVE}  filter_restart        filter_ics.10hour
${MOVE} inflate_restart       inflate_ics.10hour    # zero size
${MOVE} reg_diagnostics   reg_diagnostics.10hour

${MOVE}     True_State.nc      True_State.10hour.nc
${MOVE}     Prior_Diag.nc      Prior_Diag.10hour.nc
${MOVE} Posterior_Diag.nc  Posterior_Diag.10hour.nc

echo ""
echo "=================================================================="
echo "Now run the filter again to produce multiple restart files:"
echo "filter_restart.01"
echo "=================================================================="
echo ""

${COPY} input.nml input.nml.previous
${REMOVE} vi_script

echo ':0'                                >! vi_script
echo '/single_restart_file_out'          >> vi_script
echo ':s/true/false/'                    >> vi_script
echo ':wq'                               >> vi_script
(vi -s vi_script -e input.nml > /dev/null) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter  || exit 53

ls -lrt | tail -30    # TJH debug

echo " "
echo "=================================================================="
echo "Now do a second 10 hour run from the end of the first 10 hour"
echo "Also change the filter back to only produce single restart file"
echo "=================================================================="
echo " "

${COPY} input.nml input.nml.previous
${REMOVE} vi_script

echo ':0'                                         >! vi_script
echo '/restart_in_file_name'                      >> vi_script
echo ':s/perfect_ics.spun_up/perfect_ics.10hour'  >> vi_script
echo '/filter_nml'                                >> vi_script
echo '/restart_in_file_name'                      >> vi_script
echo ':s/filter_ics.spun_up/filter_ics.10hour'    >> vi_script
echo '/single_restart_file_out'                   >> vi_script
echo ':s/false/true/'                             >> vi_script
echo ':wq'                                        >> vi_script
(vi -s vi_script -e input.nml > /dev/null) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./perfect_model_obs  || exit 62
./filter             || exit 63

ls -lrt | tail -30    # TJH debug

# Do some tests on perfect_model_obs options first
# Set up baseline output file

${MOVE}     True_State.nc      True_State.baseline.nc
${MOVE}     obs_seq.out       obs_seq.out.baseline
${MOVE} perfect_restart   perfect_restart.baseline
${MOVE}     Prior_Diag.nc      Prior_Diag.baseline.nc
${MOVE} Posterior_Diag.nc  Posterior_Diag.baseline.nc
${MOVE}   obs_seq.final     obs_seq.final.baseline
${MOVE}  filter_restart    filter_restart.baseline
${MOVE} inflate_restart   inflate_restart.baseline    # zero size
${MOVE} reg_diagnostics   reg_diagnostics.baseline


echo ""
echo "=================================================================="
echo "Test storing the ensemble on disk instead of in core"
echo "=================================================================="
echo ""

${COPY} input.nml input.nml.previous
${REMOVE} vi_script

echo ':0'                             >! vi_script
echo '/ensemble_manager'              >> vi_script
echo '/in_core'                       >> vi_script
echo ':s/true/false/'                 >> vi_script
echo ':wq'                            >> vi_script
(vi -s vi_script -e input.nml > /dev/null) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./perfect_model_obs  || exit 72

ls -lrt | tail -30    # TJH debug

${MOVE} True_State.nc      True_State.nc.out_of_core
${MOVE} perfect_restart  perfect_restart.out_of_core
${MOVE} obs_seq.out          obs_seq.out.out_of_core
${MOVE} ens_manager_ens_file.0001.0001 ens_manager_ens_file.0001.0001.out_of_core

diff perfect_restart.out_of_core  perfect_restart.baseline || exit 74
diff     obs_seq.out.out_of_core      obs_seq.out.baseline || exit 75
#diff  True_State.nc.out_of_core    True_State.nc.baseline || exit 76

echo ""
echo "=================================================================="
echo "Test two async options"
echo "async == 0 == model advance by subroutine"
echo "async == 2 == model advance with F90 calls to shell script"
echo "async == 3 == model advance when signalled by semaphore file"
echo "Change to async=2 and in_core=true for perfect model."
echo "=================================================================="
echo ""

#-----------------------------------------------------------------------
# Need to get the scripts 
#-----------------------------------------------------------------------
${COPY} ../shell_scripts/*.csh .
${COPY} ${DARTHOME}/shell_scripts/advance_ens.csh      advance_ens.csh
${COPY} ${DARTHOME}/shell_scripts/assim_filter.csh    assim_filter.csh
${COPY} ${DARTHOME}/shell_scripts/filter_server.csh  filter_server.csh

${COPY} input.nml input.nml.previous
${REMOVE} vi_script

echo ':0'                            >! vi_script
echo '/perfect_model_obs_nml'        >> vi_script
echo '/async'                        >> vi_script
echo ':s/0/2/'                       >> vi_script
echo '/in_core'                      >> vi_script
echo ':s/false/true/'                >> vi_script
echo ':wq'                           >> vi_script
(vi -s vi_script -e input.nml > /dev/null) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./perfect_model_obs  || exit 82

ls -lrt | tail -30    # TJH debug

# This run result in:
# True_State.nc, obs_seq.out, perfect_restart
# assim_model_state_ic1
# filter_control
# assim_model_state_ud1
# integrate_model_out_temp1

${MOVE} True_State.nc       True_State.nc.2
${MOVE} obs_seq.out           obs_seq.out.2
${MOVE} perfect_restart   perfect_restart.2

diff perfect_restart.2  perfect_restart.baseline || exit 84
diff     obs_seq.out.2      obs_seq.out.baseline || exit 85
#diff  True_State.nc.2    True_State.nc.baseline || exit 86

echo ""
echo "=================================================================="
echo "Now try async == 3 (filter_server.csh driving the model advance)"
echo "=================================================================="
echo ""

${COPY} input.nml input.nml.previous
${REMOVE} vi_script

echo ':0'                            >! vi_script
echo '/perfect_model_obs_nml'        >> vi_script
echo '/async'                        >> vi_script
echo ':s/2/3/'                       >> vi_script
echo ':wq'                           >> vi_script
(vi -s vi_script -e input.nml > /dev/null) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

# TJH debug as CAM implements it ... filter_server.csh is terminated by the
# existence of go_end_filter_server ... created by filter.csh ... which is
# not part of the standard scripts ...

./filter_server.csh &
./perfect_model_obs || exit 92

ls -lrt | tail -30    # TJH debug
echo "If you're still here after 10 seconds or so ..."
echo "filter_server is not terminating when it should."
wait


# This run result in:
# True_State.nc, obs_seq.out, perfect_restart
# assim_model_state_ic1
# filter_control
# assim_model_state_ud1
# integrate_model_out_temp1

${MOVE} perfect_restart   perfect_restart.3
${MOVE} obs_seq.out           obs_seq.out.3
${COPY} True_State.nc       True_State.nc.3

diff perfect_restart.3  perfect_restart.baseline || exit 94
diff     obs_seq.out.3      obs_seq.out.baseline || exit 95
#diff  True_State.nc.3    True_State.nc.baseline || exit 96

echo ""
echo "=================================================================="
echo "Next, start checking filter options for this case"
echo "Begin by checking single versus multiple restarts"
echo "=================================================================="
echo ""

${COPY} obs_seq.out.baseline obs_seq.out

./filter  || exit 73

ls -lrt | tail -30    # TJH debug

${MOVE} obs_seq.final              obs_seq.single.final
${MOVE} filter_restart            filter_restart.single
${MOVE} inflate_restart          inflate_restart.single   # zero length
${MOVE} reg_diagnostics          reg_diagnostics.single
${MOVE} Posterior_Diag.nc      Posterior_Diag.single.nc
${COPY}     Prior_Diag.nc          Prior_Diag.single.nc

# Determine if the run was correct.
# can fire off a matlab batch job if matlab exists.
if ( $MatlabExists == 0 ) then
   echo "plot_total_err"             >! batchscript.m
   echo "print -dpsc DART_fig2.ps"   >> batchscript.m
   echo "quit"                       >> batchscript.m
   matlab -nosplash -nojvm -r batchscript
   ${REMOVE} batchscript.m
else
   echo "Matlab does not exist in your PATH ... archiving output"
   ${COPY} input.nml  ${DARTHOME}/Test2
   ${COPY} *.single.* ${DARTHOME}/Test2
endif

echo ""
echo "=================================================================="
echo "Now use multiple input files (named filter_restart.00xx)"
echo "=================================================================="
echo ""

${COPY} input.nml input.nml.previous
${REMOVE} vi_script

echo ':0'                                    >! vi_script
echo '/filter_nml'                           >> vi_script
echo '/restart_in_file_name'                 >> vi_script
echo ':s/filter_ics.10hour/filter_restart/'  >> vi_script
echo '/single_restart_file_in'               >> vi_script
echo ':s/true/false/'                        >> vi_script
echo ':wq'                                   >> vi_script
(vi -s vi_script -e input.nml > /dev/null) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter  || exit 83

ls -lrt | tail -30    # TJH debug

${MOVE}     Prior_Diag.nc        Prior_Diag.in_files.nc
${MOVE} Posterior_Diag.nc    Posterior_Diag.in_files.nc
${MOVE} obs_seq.final               obs_seq.in_files.final
${MOVE}  filter_restart      filter_restart.in_files
${MOVE} reg_diagnostics     reg_diagnostics.in_files

diff       obs_seq.in_files.final        obs_seq.single.final || exit 87
diff      filter_restart.in_files       filter_restart.single || exit 88
#diff    reg_diagnostics.in_files      reg_diagnostics.single || exit 89
#diff      Prior_Diag.nc.in_files        Prior_Diag.single.nc
#diff  Posterior_Diag.nc.in_files    Posterior_Diag.single.nc

echo ""
echo "=================================================================="
echo "Next switch the number of domains from 3 to 5"
echo "=================================================================="
echo ""

${COPY} input.nml input.nml.previous
${REMOVE} vi_script

echo ':0'                                     >! vi_script
echo '/num_domains'                           >> vi_script
echo ':s/3/5/'                                >> vi_script
echo ':wq'                                    >> vi_script
(vi -s vi_script -e input.nml > /dev/null) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter  || exit 93

ls -lrt | tail -30    # TJH debug

diff obs_seq.final                  obs_seq.single.final || exit 97
diff filter_restart          filter_restart.single       || exit 98
#diff      Prior_Diag.nc         Prior_Diag.single.nc
#diff  Posterior_Diag.nc     Posterior_Diag.single.nc

echo ""
echo "=================================================================="
echo "Next switch the number of domains back to 3; try parallel option 2"
echo "=================================================================="
echo ""

${COPY} input.nml input.nml.previous
${REMOVE} vi_script

echo ':0'                                     >! vi_script
echo '/do_parallel'                           >> vi_script
echo ':s/0/2/'                                >> vi_script
echo '/num_domains'                           >> vi_script
echo ':s/5/3/'                                >> vi_script
echo ':wq'                                    >> vi_script
(vi -s vi_script -e input.nml > /dev/null) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter || exit 103

ls -lrt | tail -30    # TJH debug

# This run results in the following:
# filter_assim_obs_seq
# filter_assim_region__in1
# filter_assim_region__in2
# filter_assim_region__in3
# filter_assim_region_out1
# filter_assim_region_out2
# filter_assim_region_out3
# assim_region_control
# assim_region_out_temp1
# assim_region_out_temp2
# assim_region_out_temp3

diff obs_seq.final                   obs_seq.single.final || exit 107
diff filter_restart           filter_restart.single       || exit 108
#diff     Prior_Diag.nc           Prior_Diag.single.nc    || exit 384
#diff Posterior_Diag.nc       Posterior_Diag.single.nc    || exit 385

echo ""
echo "=================================================================="
echo "Try parallel option 3"
echo "=================================================================="
echo ""

${COPY} input.nml input.nml.previous
${REMOVE} vi_script

echo ':0'                                     >! vi_script
echo '/do_parallel'                           >> vi_script
echo ':s/2/3/'                                >> vi_script
echo ':wq'                                    >> vi_script
(vi -s vi_script -e input.nml > /dev/null) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter_server.csh &
./filter  || exit 113

ls -lrt | tail -30    # TJH debug
echo "If you're still here after 10 seconds or so ..."
echo "filter_server is not terminating when it should."
wait

# This run results in the following output files:
# obs_seq.final, filter_restart, Prior_Diag.nc, Posterior_Diag.nc
# filter_assim_obs_seq
# assim_region_control
# filter_server.log
# filter_assim_region__in?   [1,3]
# filter_assim_region_out?   [1,3]
# assim_region_out_temp?     [1,3]

diff       obs_seq.final              obs_seq.single.final || exit 117
diff      filter_restart       filter_restart.single       || exit 118
#diff assim_tools_restart assim_tools_restart.single       || exit 119
#diff      Prior_Diag.nc           Prior_Diag.single.nc
#diff  Posterior_Diag.nc       Posterior_Diag.single.nc

echo ""
echo "=================================================================="
echo "Go back to parallel option 0 and proceed"
echo "=================================================================="
echo ""

${COPY} input.nml input.nml.previous
${REMOVE} vi_script

echo ':0'                                     >! vi_script
echo '/do_parallel'                           >> vi_script
echo ':s/3/0/'                                >> vi_script
echo ':wq'                                    >> vi_script
(vi -s vi_script -e input.nml > /dev/null) || exit 98
echo "input.nml differences from previous run:"
diff input.nml input.nml.previous
echo "input.nml now contains:"
cat input.nml
echo " "
echo "END of input.nml"

./filter  || exit 123

ls -lrt | tail -30    # TJH debug

# This run results in the following output files:
# Prior_Diag.nc, Posterior_Diag.nc, obs_seq.final, filter_restart,
# reg_diagnostics

diff obs_seq.final                    obs_seq.single.final || exit 127
diff filter_restart            filter_restart.single       || exit 128
#diff assim_tools_restart assim_tools_restart.single       || exit 129
#diff     Prior_Diag.nc            Prior_Diag.single.nc
#diff Posterior_Diag.nc        Posterior_Diag.single.nc

echo ""
echo "Testing complete  at "`date`
echo "=================================================================="
echo ""
