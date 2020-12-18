#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.
# Search below for TIMECHECK to see what times this script will
# run.

echo "`date` -- BEGIN GENERATE CAM TRUE STATE"

set nonomatch       # suppress "rm" warnings if wildcard does not match anything

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following
switch ("`hostname`")
   case be*:
      # NCAR "bluefire"
      set   MOVE = '/usr/local/bin/mv -fv'
      set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
      set   LINK = '/usr/local/bin/ln -fvs'
      set REMOVE = '/usr/local/bin/rm -fr'

      set BASEOBSDIR = /glade/proj3/image/Observations/ACARS
   breaksw

   case ys*:
      # NCAR "yellowstone"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set BASEOBSDIR = /glade/p/image/Observations/ACARS
   breaksw

   default:
      # NERSC "hopper"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set BASEOBSDIR = /scratch/scratchdirs/nscollin/ACARS
   breaksw
endsw

#-------------------------------------------------------------------------
# Determine time of model state ... from file name
# of the form "./${CASE}.cam.i.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#-------------------------------------------------------------------------

set FILE = `head -n 1 rpointer.atm`
set FILE = $FILE:r
set ATM_DATE_EXT = `echo $FILE:e`
set ATM_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
set ATM_YEAR     = `echo $ATM_DATE[1] | bc`
set ATM_MONTH    = `echo $ATM_DATE[2] | bc`
set ATM_DAY      = `echo $ATM_DATE[3] | bc`
set ATM_SECONDS  = `echo $ATM_DATE[4] | bc`
set ATM_HOUR     = `echo $ATM_DATE[4] / 3600 | bc`

echo "valid time of model is $ATM_YEAR $ATM_MONTH $ATM_DAY $ATM_SECONDS (seconds)"
echo "valid time of model is $ATM_YEAR $ATM_MONTH $ATM_DAY $ATM_HOUR (hours)"

#-------------------------------------------------------------------------
# Determine if current time is a perfect model time.
# If not, return before doing anything.
#-------------------------------------------------------------------------

## TIMECHECK:
if ( $ATM_HOUR == 0 || $ATM_HOUR == 6 || $ATM_HOUR == 12 || $ATM_HOUR == 18) then
   echo "Hour is $ATM_HOUR so we are generating perfect obs for the atmosphere"
else
   echo "Hour is not 0,6,12 or 18Z so we are skipping generating perfect obs for the atmosphere"
   echo "`date` -- END   GENERATE CAM TRUE STATE"
   exit 0
endif

#-------------------------------------------------------------------------
# Create temporary working directory for the perfect model and go there
#-------------------------------------------------------------------------

set temp_dir = pmo_cam
echo "temp_dir is $temp_dir"

if ( -d $temp_dir ) then
   ${REMOVE} $temp_dir/*
else
   mkdir -p $temp_dir
endif
cd $temp_dir

#-----------------------------------------------------------------------------
# Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of CAM.
#-----------------------------------------------------------------------------

set YYYYMM   = `printf %04d%02d                ${ATM_YEAR} ${ATM_MONTH}`
set OBSFNAME = `printf obs_seq%04d%02d%02d%02d ${ATM_YEAR} ${ATM_MONTH} ${ATM_DAY} ${ATM_HOUR}`
set OBS_FILE = ${BASEOBSDIR}/${YYYYMM}_6H/${OBSFNAME}

if (  -e   ${OBS_FILE} ) then
   ${LINK} ${OBS_FILE} obs_seq.in
else
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... no observation file $OBS_FILE"
   exit -1
endif

#=========================================================================
# Block 1: Populate a run-time directory with the input needed to run DART.
#
# DART namelist settings required:
# &perfect_model_obs_nml:  restart_in_file_name    = 'perfect_ics'
# &perfect_model_obs_nml:  obs_sequence_in_name    = 'obs_seq.in'
# &perfect_model_obs_nml:  obs_sequence_out_name   = 'obs_seq.perfect'
# &perfect_model_obs_nml:  init_time_days          = -1,
# &perfect_model_obs_nml:  init_time_seconds       = -1,
# &perfect_model_obs_nml:  first_obs_days          = -1,
# &perfect_model_obs_nml:  first_obs_seconds       = -1,
# &perfect_model_obs_nml:  last_obs_days           = -1,
# &perfect_model_obs_nml:  last_obs_seconds        = -1,
# &cam_to_dart_nml:        cam_to_dart_output_file = 'dart_ics'
#=========================================================================

if ( ! -e ${CASEROOT}/cam_input.nml ) then
   echo "ERROR ... DART required file ${CASEROOT}/cam_input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/cam_input.nml not found ... ERROR"
   exit -2
endif

sed -e "s#dart_ics#perfect_ics#" < ${CASEROOT}/cam_input.nml >! input.nml

# Turns out the .h0. files are timestamped with the START of the
# run, which is *not* ATM_DATE_EXT ...  I just link to a whatever
# is convenient (since the info is static).

set ATM_INITIAL_FILENAME = "../${CASE}.cam.i.${ATM_DATE_EXT}.nc"
set ATM_HISTORY_FILENAME = `ls -1t ../${CASE}.cam*.h0.* | head -n 1`

${LINK} $ATM_INITIAL_FILENAME caminput.nc
${LINK} $ATM_HISTORY_FILENAME cam_phis.nc

#=========================================================================
# Block 2: Convert 1 CAM restart file to a DART initial conditions file.
# At the end of the block, we have a DART initial condition file  perfect_ics
# that came from the contents of the pointer file ../rpointer.atm
#=========================================================================

echo "`date` -- BEGIN CAM-TO-DART"

${EXEROOT}/cam_to_dart

if ($status != 0) then
   echo "ERROR ... DART died in 'cam_to_dart' ... ERROR"
   echo "ERROR ... DART died in 'cam_to_dart' ... ERROR"
   exit -3
endif

echo "`date` -- END CAM-TO-DART"

#=========================================================================
# Block 3: Advance the model and harvest the synthetic observations.
# output files are:
# True_state.nc   ...... the DART state
# obs_seq.perfect ...... the synthetic observations
# dart_log.out    ...... run-time output of all DART routines
# perfect_restart ...... which we don't need
#=========================================================================

echo "`date` -- BEGIN CAM PERFECT_MODEL_OBS"

${EXEROOT}/perfect_model_obs_cam

if ($status != 0) then
   echo "ERROR ... DART died in 'perfect_model_obs_cam' ... ERROR"
   echo "ERROR ... DART died in 'perfect_model_obs_cam' ... ERROR"
   exit -4
endif

${MOVE} True_State.nc    ../cam_True_State.${ATM_DATE_EXT}.nc
${MOVE} obs_seq.perfect  ../cam_obs_seq.${ATM_DATE_EXT}.perfect
${MOVE} dart_log.out     ../cam_dart_log.${ATM_DATE_EXT}.out

echo "`date` -- END   CAM PERFECT_MODEL_OBS"

#=========================================================================
# Block 4: Update the cam restart file
#=========================================================================

# not needed ... perfect_model_obs does not update the model state.

#=========================================================================
# Block 5: Link the next cam file to the static name needed here.
#=========================================================================

cd ${RUNDIR}

set ATM_INITIAL_FILENAME = ${CASE}.cam.i.${ATM_DATE_EXT}.nc

${LINK} ${ATM_INITIAL_FILENAME} cam_initial.nc || exit -5

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

# Eat the cookie regardless
${REMOVE} ../cam_inflation_cookie
${REMOVE} perfect_ics dart_log.nml

echo "`date` -- END   GENERATE CAM TRUE STATE"

exit 0


