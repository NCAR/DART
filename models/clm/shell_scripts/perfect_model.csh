#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

echo "`date` -- BEGIN GENERATE CLM TRUE STATE"

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

      set BASEOBSDIR = /glade/proj3/image/Observations/FluxTower
   breaksw

   case ys*:
      # NCAR "yellowstone"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set BASEOBSDIR = /glade/p/image/Observations/land
   breaksw

   case lone*:
      # UT lonestar
      set   MOVE = '/bin/mv -fv'
      set   COPY = '/bin/cp -fv --preserve=timestamps'
      set   LINK = '/bin/ln -fvs'
      set REMOVE = '/bin/rm -fr'

      set BASEOBSDIR = ${WORK}/DART/observations/snow/work/obs_seqs
   breaksw

   case la*:
      # LBNL "lawrencium"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set BASEOBSDIR = /your/observation/directory/here
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
# of the form "./${CASE}.clm2.r.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#-------------------------------------------------------------------------

set FILE = `head -n 1 rpointer.lnd`
set FILE = $FILE:r
set LND_DATE_EXT = `echo $FILE:e`
set LND_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
set LND_YEAR     = `echo $LND_DATE[1] | bc`
set LND_MONTH    = `echo $LND_DATE[2] | bc`
set LND_DAY      = `echo $LND_DATE[3] | bc`
set LND_SECONDS  = `echo $LND_DATE[4] | bc`
set LND_HOUR     = `echo $LND_DATE[4] / 3600 | bc`

echo "valid time of model is $LND_YEAR $LND_MONTH $LND_DAY $LND_SECONDS (seconds)"
echo "valid time of model is $LND_YEAR $LND_MONTH $LND_DAY $LND_HOUR (hours)"

#-------------------------------------------------------------------------
# Create temporary working directory for the perfect model and go there
#-------------------------------------------------------------------------

set temp_dir = pmo_clm
echo "temp_dir is $temp_dir"

if ( -d $temp_dir ) then
   ${REMOVE} $temp_dir/*
else
   mkdir -p $temp_dir
endif
cd $temp_dir

#-----------------------------------------------------------------------------
# Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of CLM.
#
# The CLM observations are stowed in two sets of directories.
# If you are stopping every 24 hours or more, the obs are in directories like YYYYMM.
# In all other situations the observations come from directories like YYYYMM_6H.
# The only ugly part here is if the first advance and subsequent advances are
# not the same length. The observations _may_ come from different directories.
#
# The contents of the file must match the history file contents if one is using
# the obs_def_tower_mod or could be the 'traditional' +/- 12Z ... or both.
# Since the history file contains the previous days' history ... so must the obs file.
#-----------------------------------------------------------------------------

if ($STOP_N >= 24) then
   set OBSDIR = `printf %04d%02d    ${LND_YEAR} ${LND_MONTH}`
else
   set OBSDIR = `printf %04d%02d_6H ${LND_YEAR} ${LND_MONTH}`
endif

set OBS_FILE = ${BASEOBSDIR}/${OBSDIR}/obs_seq.${LND_DATE_EXT}

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
# &clm_to_dart_nml:        clm_to_dart_output_file = 'dart_ics'
# &model_nml:              clm_restart_filename        = 'clm_restart.nc'
# &model_nml:              clm_history_filename        = 'clm_history.nc'
# &model_nml:              clm_vector_history_filename = 'clm_vector_history.nc'
#=========================================================================

if ( ! -e ${CASEROOT}/input.nml ) then
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   exit -2
endif

sed -e "s#dart_ics#perfect_ics#" < ${CASEROOT}/input.nml >! input.nml

# DART/CLM routines all need a clm_restart.nc, clm_history.nc, etc.
# The flux tower forward operator looks for a CLM history file with
# an instance number in the filename.

set      LND_RESTART_FILENAME = ${CASE}.clm2.r.${LND_DATE_EXT}.nc
set      LND_HISTORY_FILENAME = ${CASE}.clm2.h0.${LND_DATE_EXT}.nc
set  LND_VEC_HISTORY_FILENAME = ${CASE}.clm2.h2.${LND_DATE_EXT}.nc
set     OBS1_HISTORY_FILENAME = ${CASE}.clm2.h1.${LND_DATE_EXT}.nc
set     OBS2_HISTORY_FILENAME = ${CASE}.clm2_0001.h1.${LND_DATE_EXT}.nc

${LINK} ../$LND_RESTART_FILENAME clm_restart.nc
${LINK} ../$LND_HISTORY_FILENAME clm_history.nc

if (  -e   ../$OBS1_HISTORY_FILENAME) then
   ${LINK} ../$OBS1_HISTORY_FILENAME $OBS2_HISTORY_FILENAME
endif
if (  -e   ../$LND_VEC_HISTORY_FILENAME) then
   ${LINK} ../$LND_VEC_HISTORY_FILENAME clm_vector_history.nc
endif

#=========================================================================
# Block 2: Convert 1 CLM restart file to a DART initial conditions file.
# At the end of the block, we have a DART initial condition file  perfect_ics
# that came from the contents of the pointer file ../rpointer.lnd
#=========================================================================

echo "`date` -- BEGIN CLM-TO-DART"

${EXEROOT}/clm_to_dart

if ($status != 0) then
   echo "ERROR ... DART died in 'clm_to_dart' ... ERROR"
   echo "ERROR ... DART died in 'clm_to_dart' ... ERROR"
   exit -3
endif

echo "`date` -- END CLM-TO-DART"

#=========================================================================
# Block 3: Advance the model and harvest the synthetic observations.
# output files are:
# True_state.nc   ...... the DART state
# obs_seq.perfect ...... the synthetic observations
# dart_log.out    ...... run-time output of all DART routines
# perfect_restart ...... which we don't need
#=========================================================================

echo "`date` -- BEGIN CLM PERFECT_MODEL_OBS"

${EXEROOT}/perfect_model_obs

if ($status != 0) then
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   exit -4
endif

${MOVE} True_State.nc    ../clm_True_State.${LND_DATE_EXT}.nc
${MOVE} obs_seq.perfect  ../clm_obs_seq.${LND_DATE_EXT}.perfect
${MOVE} dart_log.out     ../clm_dart_log.${LND_DATE_EXT}.out

echo "`date` -- END   CLM PERFECT_MODEL_OBS"

#=========================================================================
# Block 4: Update the clm restart file
#=========================================================================

# not needed ... perfect_model_obs does not update the model state.

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

${REMOVE} perfect_ics dart_log.nml

echo "`date` -- END   GENERATE CLM TRUE STATE"

exit 0


