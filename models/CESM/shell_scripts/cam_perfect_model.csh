#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

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
# Determine time of model state ... from file name of first member
# of the form "./${CASE}.cam_${ensemble_member}.i.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#-------------------------------------------------------------------------

set FILE = `head -n 1 rpointer.atm_0001`
set FILE = $FILE:t
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
   echo "`date` -- END CAM_ASSIMILATE"
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

set YYYYMM   = `printf %04d%02d ${ATM_YEAR} ${ATM_MONTH}`
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
#=========================================================================

echo "`date` -- BEGIN COPY BLOCK"

if (  -e   ${CASEROOT}/cam_input.nml ) then
   ${COPY} ${CASEROOT}/cam_input.nml input.nml
else
   echo "ERROR ... DART required file ${CASEROOT}/cam_input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/cam_input.nml not found ... ERROR"
   exit -2
endif

echo "`date` -- END COPY BLOCK"

#=========================================================================
# Block 2: Convert 1 CAM restart file to DART initial conditions file.
# At the end of the block, we have DART initial condition file  perfect_ics
# that came from pointer file ../rpointer.atm_0001
#
# REQUIRED DART namelist settings:
# &perfect_model_obs_nml:  restart_in_file_name    = 'perfect_ics'
# &cam_to_dart_nml:        cam_to_dart_output_file = 'dart_ics'
#=========================================================================

echo "`date` -- BEGIN CAM-TO-DART"

set member = 1

   set MYTEMPDIR = member_${member}
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   # Turns out the .h0. files are timestamped with the START of the
   # run, which is *not* ATM_DATE_EXT ...  I just link to a whatever
   # is convenient (since the info is static).
   # make sure there are no old output logs hanging around
   $REMOVE output.${member}.cam_to_dart

   set ATM_INITIAL_FILENAME = "../../${CASE}.cam.i.${ATM_DATE_EXT}.nc"
   set ATM_HISTORY_FILENAME = `ls -1t ../../${CASE}.cam*.h0.* | head -n 1`
   set     DART_IC_FILENAME = `printf filter_ics.%04d     ${member}`
   set    DART_RESTART_FILE = `printf filter_restart.%04d ${member}`

   sed -e "s#dart_ics#../${DART_IC_FILENAME}#" \
       -e "s#dart_restart#../${DART_RESTART_FILE}#" < ../input.nml >! input.nml

   ${LINK} $ATM_INITIAL_FILENAME caminput.nc
   ${LINK} $ATM_HISTORY_FILENAME cam_phis.nc

   ${EXEROOT}/cam_to_dart >! output.${member}.cam_to_dart 

   if ($status != 0) then
   echo "ERROR ... DART died in 'cam_to_dart' ... ERROR"
   echo "ERROR ... DART died in 'cam_to_dart' ... ERROR"
   exit -6
endif

echo "`date` -- END CAM-TO-DART"

#=========================================================================
# Block 3: Advance the model and harvest the synthetic observations.
# Will result in a single file : 'perfect_restart' which we don't need
# for a perfect model experiment with CESM.
#
# DART namelist settings required:
# &perfect_model_obs_nml:           async                  = 0,
# &perfect_model_obs_nml:           adv_ens_command        = "./no_model_advance.csh",
# &perfect_model_obs_nml:           restart_in_file_name   = 'perfect_ics'
# &perfect_model_obs_nml:           restart_out_file_name  = 'perfect_restart'
# &perfect_model_obs_nml:           obs_sequence_in_name   = 'obs_seq.in'
# &perfect_model_obs_nml:           obs_sequence_out_name  = 'obs_seq.out'
# &perfect_model_obs_nml:           init_time_days         = -1,
# &perfect_model_obs_nml:           init_time_seconds      = -1,
# &perfect_model_obs_nml:           first_obs_days         = -1,
# &perfect_model_obs_nml:           first_obs_seconds      = -1,
# &perfect_model_obs_nml:           last_obs_days          = -1,
# &perfect_model_obs_nml:           last_obs_seconds       = -1,
#
#=========================================================================

# CAM:static_init_model() always needs a caminput.nc and a cam_phis.nc
# for geometry information, etc.

set ATM_INITIAL_FILENAME = ../${CASE}.cam.i.${ATM_DATE_EXT}.nc
set ATM_HISTORY_FILENAME = `ls -1t ../${CASE}.cam*.h0.* | head -n 1`

${LINK} $ATM_INITIAL_FILENAME caminput.nc
${LINK} $ATM_HISTORY_FILENAME cam_phis.nc

echo "`date` -- BEGIN CAM PERFECT_MODEL_OBS"
${LAUNCHCMD} ${EXEROOT}/perfect_model_obs_cam || exit -7
echo "`date` -- END   CAM PERFECT_MODEL_OBS"

${MOVE} Prior_Diag.nc      ../cam_Prior_Diag.${ATM_DATE_EXT}.nc
${MOVE} Posterior_Diag.nc  ../cam_Posterior_Diag.${ATM_DATE_EXT}.nc
${MOVE} obs_seq.final      ../cam_obs_seq.${ATM_DATE_EXT}.final
${MOVE} dart_log.out       ../cam_dart_log.${ATM_DATE_EXT}.out

#=========================================================================
# Block 4: Update the cam restart file
#=========================================================================

# not needed ... perfect_model_obs does not update the model state.

#=========================================================================
# Block 5: Link the next cam file to the static name needed here.
#=========================================================================

cd ${RUNDIR}

   set ATM_INITIAL_FILENAME = ${CASE}.cam.i.${ATM_DATE_EXT}.nc

   ${LINK} ${ATM_INITIAL_FILENAME} cam_initial.nc || exit -9

end

#-------------------------------------------------------------------------
# Cleanup
# we (DART) do not need these files, and CESM does not need them either
# to continue a run.  if we remove them here they do not get moved to
# the short-term archiver.
#-------------------------------------------------------------------------

# ${REMOVE} ${RUNDIR}/*.rh0.*
# ${REMOVE} ${RUNDIR}/*.rs1.*

echo "`date` -- END   GENERATE CAM TRUE STATE"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

