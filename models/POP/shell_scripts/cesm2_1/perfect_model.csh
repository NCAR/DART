#!/bin/csh
#
# Copyright 2020 University Corporation for Atmospheric Research
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License.
# Please view the License at http://www.apache.org/licenses/LICENSE-2.0
#
#=========================================================================
# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.
#=========================================================================

echo "`date` -- BEGIN GENERATE POP TRUE STATE"

# As of CESM2.0, the perfect_model.csh is called by CESM - and has
# two arguments: the CASEROOT and the DATA_ASSIMILATION_CYCLE

setenv CASEROOT $1
setenv ASSIMILATION_CYCLE $2

source ${CASEROOT}/DART_params.csh || exit 1

# Python uses C indexing on loops; cycle = [0,....,$DATA_ASSIMILATION_CYCLES - 1]
# "Fix" that here, so the rest of the script isn't confusing.
@ cycle = $ASSIMILATION_CYCLE + 1

# xmlquery must be executed in $CASEROOT.
cd ${CASEROOT}
setenv CASE                     `./xmlquery CASE                     --value`
setenv ENSEMBLE_SIZE            `./xmlquery NINST_LND                --value`
setenv EXEROOT                  `./xmlquery EXEROOT                  --value`
setenv RUNDIR                   `./xmlquery RUNDIR                   --value`
setenv ARCHIVE                  `./xmlquery DOUT_S_ROOT              --value`
setenv TOTALPES                 `./xmlquery TOTALPES                 --value`
setenv STOP_N                   `./xmlquery STOP_N                   --value`
setenv DATA_ASSIMILATION_CYCLES `./xmlquery DATA_ASSIMILATION_CYCLES --value`
setenv TASKS_PER_NODE           `./xmlquery MAX_TASKS_PER_NODE       --value`

#=========================================================================
# Block 1: Determine time of model state ... from file name
# of the form "./${CASE}.pop.r.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#=========================================================================

cd ${RUNDIR}
set FILE = `head -n 1 rpointer.ocn.restart`
# :t is a csh modifier that returns the "tail" of the filename
# essentially stripping away the filepath that prepends the file
# :r is a csh modifier that returns the "root" of the filename 
# essentially stripping away the file suffix
set FILE = $FILE:t:r
# :e is a csh modifier that returns the "extension" of the filename
# since we stripped away the .nc in the previous command, the current
# extension is the datetime in the form YYYY-MM-DD-SSSSS
set OCN_DATE_EXT = $FILE:e
# Piping the datetime to this sed command replaces the hyphens with spaces
set OCN_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
# Piping the indexed substrings to bc strips away the leading zeros
set OCN_YEAR     = `echo $OCN_DATE[1] | bc`
set OCN_MONTH    = `echo $OCN_DATE[2] | bc`
set OCN_DAY      = `echo $OCN_DATE[3] | bc`
set OCN_SECONDS  = `echo $OCN_DATE[4] | bc`
set OCN_HOUR     = `echo $OCN_DATE[4] / 3600 | bc`

echo "valid time of model is $OCN_YEAR $OCN_MONTH $OCN_DAY $OCN_SECONDS (seconds)"
echo "valid time of model is $OCN_YEAR $OCN_MONTH $OCN_DAY $OCN_HOUR (hours)"

#=========================================================================
# Block 2: Populate a run-time directory with the obs_seq needed to run DART.
#=========================================================================

# Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of POP.

set YYYYMM   = `printf %04d%02d                ${OCN_YEAR} ${OCN_MONTH}`
set OBSFNAME = `printf obs_seq.0Z.%04d%02d%02d ${OCN_YEAR} ${OCN_MONTH} ${OCN_DAY}`
set OBS_FILE = ${BASEOBSDIR}/${YYYYMM}/${OBSFNAME}

if (  -e   ${OBS_FILE} ) then
   ${LINK} ${OBS_FILE} obs_seq.in
else
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... no observation file $OBS_FILE"
   exit -2
endif

#=========================================================================
# Block 3: Populate a run-time directory with the input.nml needed to run DART.
#=========================================================================

echo "`date` -- BEGIN COPY BLOCK"

if ( -e ${CASEROOT}/input.nml ) then
   ${COPY} ${CASEROOT}/input.nml .
else
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   exit -3
endif

echo "`date` -- END COPY BLOCK"

set OCN_RESTART_FILENAME = ${CASE}.pop.r.${OCN_DATE_EXT}.nc
${LINK} ./$OCN_RESTART_FILENAME pop.r.nc

#=========================================================================
# Block 4: Advance the model and harvest the synthetic observations.
# 
# Please consult the perfect_model_obs_nml in input.nml in order to 
# interpret this section.
# 
# if write_output_state_to_file = .true. then DART outputs a state file.
# This file is named using the value of output_state_files.
# The default setting is output_state_files = 'perfect_restart.nc' 
# 
# The synthetic observations are also output. This file is named using 
# the value of obs_seq_out_file_name. 
# The default setting is obs_seq_out_filename = 'obs_seq.perfect'
# 
# Additionally, dart_log.out is created. It contains run-time output
# of all DART routines.
#
# If input.nml:perfect_model_obs_nml:single_file_out = .true. , another
# output file is created named true_state.nc. It contains the DART state.
# This is identical to the original POP state but only contains the 
# variables that are part of the DART model state and some additional
# metadata to interpret locations, etc. The DART state-space diagnostics
# reference this file.
#=========================================================================

echo "`date` -- BEGIN POP PERFECT_MODEL_OBS"

${LAUNCHCMD} ${EXEROOT}/perfect_model_obs

if ($status != 0) then
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   exit -4
endif

${MOVE} obs_seq.perfect  ../pop_obs_seq.${OCN_DATE_EXT}.perfect
${MOVE} dart_log.out     ../pop_dart_log.${OCN_DATE_EXT}.out
${MOVE} true_state.nc    ../pop_true_state.${OCN_DATE_EXT}.nc

echo "`date` -- END   POP PERFECT_MODEL_OBS"

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

${REMOVE} dart_log.nml

echo "`date` -- END   GENERATE POP TRUE STATE"

exit 0
