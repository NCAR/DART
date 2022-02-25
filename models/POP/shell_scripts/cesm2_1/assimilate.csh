#!/bin/csh
#
# Copyright 2020 University Corporation for Atmospheric Research
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License.
# Please view the License at http://www.apache.org/licenses/LICENSE-2.0
#
# This script performs an assimilation by directly reading and writing to
# the POP restart file. There is no post-processing step 'dart_to_pop'.

#=========================================================================
# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.
#=========================================================================

echo "`date` -- BEGIN POP_ASSIMILATE"

# As of CESM2.0, the assimilate.csh is called by CESM - and has
# two arguments: the CASEROOT and the DATA_ASSIMILATION_CYCLE

setenv CASEROOT $1
setenv ASSIMILATION_CYCLE $2

source ${CASEROOT}/DART_params.csh || exit 1

# Python uses C indexing on loops; cycle = [0,....,$DATA_ASSIMILATION_CYCLES - 1]
# "Fix" that here, so the rest of the script isn't confusing.
@ cycle = $ASSIMILATION_CYCLE + 1

#=========================================================================
# Get the environmental variables needed for the script
#=========================================================================
cd ${CASEROOT}
setenv CASE               `./xmlquery CASE             --value`
setenv ENSEMBLE_SIZE      `./xmlquery NINST_OCN        --value`
setenv EXEROOT            `./xmlquery EXEROOT          --value`
setenv RUNDIR             `./xmlquery RUNDIR           --value`
setenv ARCHIVE            `./xmlquery DOUT_S_ROOT      --value`
setenv TOTALPES           `./xmlquery TOTALPES         --value`
setenv STOP_N             `./xmlquery STOP_N           --value`

#=========================================================================
# Block 1: Determine time of model state ... from file name of first member
# of the form "./${CASE}.pop_${ensemble_member}.r.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#=========================================================================
cd ${RUNDIR}
set FILE = `head -n 1 rpointer.ocn_0001.restart`
set FILE = $FILE:r
set OCN_DATE_EXT = `echo $FILE:e`
set OCN_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
set OCN_YEAR     = `echo $OCN_DATE[1] | bc`
set OCN_MONTH    = `echo $OCN_DATE[2] | bc`
set OCN_DAY      = `echo $OCN_DATE[3] | bc`
set OCN_SECONDS  = `echo $OCN_DATE[4] | bc`
set OCN_HOUR     = `echo $OCN_DATE[4] / 3600 | bc`

echo "valid time of model is $OCN_YEAR $OCN_MONTH $OCN_DAY $OCN_SECONDS (seconds)"
echo "valid time of model is $OCN_YEAR $OCN_MONTH $OCN_DAY $OCN_HOUR (hours)"

#=========================================================================
# Block 2: Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of POP.
#=========================================================================

set YYYYMMDD = `printf %04d%02d%02d ${OCN_YEAR} ${OCN_MONTH} ${OCN_DAY}`
set YYYYMM   = `printf %04d%02d     ${OCN_YEAR} ${OCN_MONTH}`
set OBSFNAME = obs_seq.0Z.${YYYYMMDD}
set OBS_FILE = ${BASEOBSDIR}/${YYYYMM}/${OBSFNAME}

if ( -e      obs_seq.out ) then
   ${REMOVE} obs_seq.out
endif

if (  -e   ${OBS_FILE} ) then
   ${LINK} ${OBS_FILE} obs_seq.out
else
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... no observation file $OBS_FILE"
   exit -2
endif

#=========================================================================
# Block 3: Populate a run-time directory with the input needed to run DART.
#=========================================================================

echo "`date` -- BEGIN COPY BLOCK"

if (  -e   ${CASEROOT}/input.nml ) then
   ${COPY} ${CASEROOT}/input.nml .
else
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   exit -3
endif

echo "`date` -- END COPY BLOCK"


#=========================================================================
# Block 4: REQUIRED DART namelist settings
#=========================================================================
${REMOVE} restarts_in.txt restarts_out.txt

foreach FILE ( rpointer.ocn_????.restart )
   head -n 1 ${FILE} >> restarts_in.txt
end

# WARNING: this is the part where the files just get overwritten
${COPY} restarts_in.txt restarts_out.txt

# Filter needs a pop_in and a pop.r.nc to start.
# Lots of ways to get the filename

set OCN_RESTART_FILENAME = `head -n 1 rpointer.ocn_0001.restart`

if ( -e pop.r.nc  ) then
    $REMOVE pop.r.nc
endif
if ( -e pop_in  ) then
    $REMOVE pop_in
endif

${LINK} $OCN_RESTART_FILENAME pop.r.nc
${LINK} pop_in_0001           pop_in

#=========================================================================
# Block 5: DART INFLATION
# This stages the files that contain the inflation values.
# The inflation values change through time and should be archived.
#
# This file is only relevant if 'inflation' is turned on -
# i.e. if inf_flavor(1) /= 0 AND inf_initial_from_restart = .TRUE.
#
# filter_nml
# inf_flavor                  = 2,                       0,
# inf_initial_from_restart    = .true.,                  .false.,
#
# The inflation file is essentially a duplicate of the DART model state ...
# For the purpose of this script, they are the output of a previous assimilation,
# so they should be named something like prior_inflate_restart.YYYY-MM-DD-SSSSS
#
# NOTICE: inf_initial_from_restart and inf_sd_initial_from_restart are somewhat
# problematic. During the bulk of an experiment, these should be TRUE, since
# we want to read existing inflation files. However, the first assimilation
# might need these to be FALSE and then subsequently be set to TRUE.
# There are two ways to handle this:
#
# 1) Create the initial files offline (perhaps with 'fill_inflation_restart')
#    and stage them with the appropriate names in the RUNDIR.
#    You must manually remove the pop_inflation_cookie file
#    from the RUNDIR in this case.
#    - OR -
# 2) create a cookie file called RUNDIR/pop_inflation_cookie
#    The existence of this file will cause this script to set the
#    namelist appropriately. This script will 'eat' the cookie file
#    to prevent this from happening for subsequent executions. If the
#    inflation file does not exist for them, and it needs to, this script
#    should die. The CESM_DART_config.csh script automatically creates a
#    cookie file to support this option.
#
# The strategy is to use the LATEST inflation file from the CESM 'rundir'.
# After an assimilation, the new inflation values/files will be moved to
# the CESM rundir to be used for subsequent assimilations. If the short-term
# archiver has worked correctly, only the LATEST files will available. Of
# course, it is not required to have short-term archiving turned on, so ...
#=========================================================================

set  MYSTRING = `grep inf_flavor input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  PRIOR_INF = $MYSTRING[2]
set  POSTE_INF = $MYSTRING[3]

set  MYSTRING = `grep inf_initial_from_restart input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  PRIOR_TF = `echo $MYSTRING[2] | tr '[:upper:]' '[:lower:]'`
set  POSTE_TF = `echo $MYSTRING[3] | tr '[:upper:]' '[:lower:]'`

# IFF we want PRIOR inflation:

if ( $PRIOR_INF > 0 ) then

   if ($PRIOR_TF == false) then
      # we are not using an existing inflation file.
      echo "inf_flavor(1) = $PRIOR_INF, using namelist values."

   else if ( -e pop_inflation_cookie ) then
      # We want to use an existing inflation file, but this is
      # the first assimilation so there is no existing inflation
      # file. This is the signal we need to to coerce the namelist
      # to have different values for this execution ONLY.
      # Since the local namelist comes from CASEROOT each time, we're golden.

      set PRIOR_TF = FALSE

      if (-x ${EXEROOT}/fill_inflation_restart) then
         ${EXEROOT}/fill_inflation_restart
      else
         echo "ERROR: Requested PRIOR inflation restart for the first cycle."
         echo "       There are no existing inflation files available "
         echo "       and ${EXEROOT}/fill_inflation_restart is missing."
         echo "EXITING"
         exit -3
      endif

ex input.nml <<ex_end
g;inf_initial_from_restart ;s;= .*;= .${PRIOR_TF}., .${POSTE_TF}.,;
g;inf_sd_initial_from_restart ;s;= .*;= .${PRIOR_TF}., .${POSTE_TF}.,;
wq
ex_end

   else
      # Look for the output from the previous assimilation

      # Checking for a prior inflation file to use

      (ls -rt1 ${CASE}.pop.output_priorinf_mean.* | tail -n 1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      # If one exists, use it as input for this assimilation
      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest input_priorinf_mean.nc
      else
         echo "ERROR: Requested PRIOR inflation but specified no incoming inflation file."
         echo "ERROR: expected something like ${CASE}.pop.output_priorinf_mean.YYYY-MM-DD-SSSSS.nc"
         exit 127
      endif

      # Checking for a prior inflation sd file to use

      (ls -rt1 ${CASE}.pop.output_priorinf_sd.* | tail -n 1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest input_priorinf_sd.nc
      else
         echo "ERROR: Requested PRIOR inflation but specified no incoming inflation SD file."
         echo "ERROR: expected something like ${CASE}.pop.output_priorinf_sd.YYYY-MM-DD-SSSSS.nc"
         exit 2
      endif

   endif
else
   echo "Prior Inflation           not requested for this assimilation."
endif

# POSTERIOR: We look for the 'newest' and use it - IFF we need it.

if ( $POSTE_INF > 0 ) then

   if ($POSTE_TF == false) then
      # we are not using an existing inflation file.
      echo "inf_flavor(2) = $POSTE_INF, using namelist values."

   else if ( -e pop_inflation_cookie ) then
      # We want to use an existing inflation file, but this is
      # the first assimilation so there is no existing inflation
      # file. This is the signal we need to to coerce the namelist
      # to have different values for this execution ONLY.
      # Since the local namelist comes from CASEROOT each time, we're golden.

      set POSTE_TF = FALSE

      if (-x ${EXEROOT}/fill_inflation_restart) then
         ${EXEROOT}/fill_inflation_restart
         ${MOVE} input_priorinf_mean.nc input_postinf_mean.nc || exit 125
         ${MOVE} input_priorinf_sd.nc   input_postinf_sd.nc   || exit 126

      else
         echo "ERROR: Requested POSTERIOR inflation restart for the first cycle."
         echo "       There are no existing inflation files available "
         echo "       and ${EXEROOT}/fill_inflation_restart is missing."
         echo "EXITING"
         exit 127
      endif

ex input.nml <<ex_end
g;inf_initial_from_restart ;s;= .*;= .${PRIOR_TF}., .${POSTE_TF}.,;
g;inf_sd_initial_from_restart ;s;= .*;= .${PRIOR_TF}., .${POSTE_TF}.,;
wq
ex_end

   else
      # Look for the output from the previous assimilation
      # Checking for a posterior inflation file to use

      (ls -rt1 ${CASE}.pop.output_postinf_mean.* | tail -n 1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest input_postinf_mean.nc
      else
         echo "ERROR: Requested POSTERIOR inflation but specified no incoming inflation file."
         echo "ERROR: expected something like ${CASE}.pop.output_postinf_mean.YYYY-MM-DD-SSSSS.nc"
         exit -5
      endif

      # Checking for a posterior inflation sd file to use

      (ls -rt1 ${CASE}.pop.output_postinf_sd.* | tail -n 1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest input_postinf_sd.nc
      else
         echo "ERROR: Requested POSTERIOR inflation but specified no incoming inflation SD file."
         echo "ERROR: expected something like ${CASE}.pop.output_postinf_sd.YYYY-MM-DD-SSSSS.nc"
         exit 2
      endif

   endif
else
   echo "Posterior Inflation       not requested for this assimilation."
endif

${REMOVE} pop_inflation_cookie



#=========================================================================
# Block 6: Actually run the assimilation. 
# WARNING: this version just overwrites the input - no ability to recover
#
# DART namelist settings required:
# &filter_nml:           async                    = 0,
# &filter_nml:           adv_ens_command          = "no_CESM_advance_script",
# &filter_nml:           input_state_file_list    = "restarts_in.txt"
# &filter_nml:           output_state_file_list   = "restarts_out.txt"
# &filter_nml:           stages_to_write          = 'preassim', 'output'
# &filter_nml:           output_restarts          = .true.
# &filter_nml:           output_mean              = .true.
# &filter_nml:           output_sd                = .true.
# &filter_nml:           write_all_stages_at_end  = .true.
#
# &filter_nml:           obs_sequence_in_name    = 'obs_seq.out'
# &filter_nml:           obs_sequence_out_name   = 'obs_seq.final'
# &filter_nml:           init_time_days          = -1,
# &filter_nml:           init_time_seconds       = -1,
# &filter_nml:           first_obs_days          = -1,
# &filter_nml:           first_obs_seconds       = -1,
# &filter_nml:           last_obs_days           = -1,
# &filter_nml:           last_obs_seconds        = -1,
# &ensemble_manager_nml: single_restart_file_in  = .false.
# &ensemble_manager_nml: single_restart_file_out = .false.
#
#=========================================================================

echo "`date` -- BEGIN FILTER"
${LAUNCHCMD} ${EXEROOT}/filter || exit -6
echo "`date` -- END FILTER"

#=========================================================================
# Block 7: Tag the output with the valid time of the model state
#=========================================================================

foreach FILE ( preassim_mean.nc           preassim_sd.nc \
               preassim_priorinf_mean.nc  preassim_priorinf_sd.nc \
               preassim_postinf_mean.nc   preassim_postinf_sd.nc \
               postassim_mean.nc          postassim_sd.nc \
               postassim_priorinf_mean.nc postassim_priorinf_sd.nc \
               postassim_postinf_mean.nc  postassim_postinf_sd.nc \
               output_mean.nc             output_sd.nc \
               output_priorinf_mean.nc    output_priorinf_sd.nc \
               output_postinf_mean.nc     output_postinf_sd.nc )

   if ( -e $FILE ) then
      set FBASE = $FILE:r
      set FEXT = $FILE:e
      ${MOVE} ${FILE} ${CASE}.pop.${FBASE}.${OCN_DATE_EXT}.${FEXT}
   endif
end

# Tag the observation file and run-time output

${MOVE} obs_seq.final   ${CASE}.pop.obs_seq.final.${OCN_DATE_EXT}
${MOVE} dart_log.out    ${CASE}.pop.dart_log.${OCN_DATE_EXT}.out

#=========================================================================
# Cleanup
#=========================================================================

echo "`date` -- END POP_ASSIMILATE"

exit 0
