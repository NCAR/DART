#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.
# Search below for TIMECHECK to see what times this script will
# assimilate.

echo "`date` -- BEGIN CAM_ASSIMILATE"

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
      set  LAUNCHCMD = mpirun.lsf
   breaksw

   case ys*:
      # NCAR "yellowstone"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'
      set TASKS_PER_NODE = `echo $LSB_SUB_RES_REQ | sed -ne '/ptile/s#.*\[ptile=\([0-9][0-9]*\)]#\1#p'`
      setenv MP_DEBUG_NOTIMEOUT yes

      set BASEOBSDIR = /glade/p/image/Observations/ACARS
      set  LAUNCHCMD = mpirun.lsf
   breaksw

   default:
      # NERSC "hopper"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set BASEOBSDIR = /scratch/scratchdirs/nscollin/ACARS
      set  LAUNCHCMD = "aprun -n $NTASKS"
   breaksw
endsw

set ensemble_size = ${NINST_ATM}

#-------------------------------------------------------------------------
# Determine time of model state ... from file name of first member
# of the form "./${CASE}.cam_${ensemble_member}.i.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#-------------------------------------------------------------------------

set FILE = `head -n 1 rpointer.atm_0001`
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
# Determine if current time is an assimilation time.
# If not, return before assimilating.
#-------------------------------------------------------------------------

## TIMECHECK:
if ( $ATM_HOUR == 0 || $ATM_HOUR == 6 || $ATM_HOUR == 12 || $ATM_HOUR == 18) then
   echo "Hour is $ATM_HOUR so we are assimilating the atmosphere"
else
   echo "Hour is not 0,6,12 or 18Z so we are skipping the atmosphere assimilation"
   echo "`date` -- END CAM_ASSIMILATE"
   exit 0
endif

#-------------------------------------------------------------------------
# Create temporary working directory for the assimilation and go there
#-------------------------------------------------------------------------

set temp_dir = assimilate_cam
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
   ${LINK} ${OBS_FILE} obs_seq.out
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

# If possible, use the round-robin approach to deal out the tasks.
# Since the ensemble manager is not used by cam_to_dart or dart_to_cam,
# it is OK to set it here and have it used by all routines.

if ($?TASKS_PER_NODE) then
   if ($#TASKS_PER_NODE > 0) then
      ${COPY} input.nml input.nml.$$
      sed -e "s#layout.*#layout = 2#" \
          -e "s#tasks_per_node.*#tasks_per_node = $TASKS_PER_NODE#" input.nml.$$ >! input.nml
      ${REMOVE} input.nml.$$
   endif
endif

#=========================================================================
# Block 2: Stage the files needed for SAMPLING ERROR CORRECTION
#
# The sampling error correction is a lookup table.
# The tables were originally in the DART distribution, but should
# have been staged to $CASEROOT at setup time.
# Each ensemble size has its own (static) file.
# It is only needed if
# input.nml:&assim_tools_nml:sampling_error_correction = .true.,
#=========================================================================

set  MYSTRING = `grep sampling_error_correction input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set SECSTRING = `echo $MYSTRING[2] | tr '[:upper:]' '[:lower:]'`

if ( $SECSTRING == true ) then
   set SAMP_ERR_FILE = ${CASEROOT}/final_full.${ensemble_size}
   if (  -e   ${SAMP_ERR_FILE} ) then
      ${COPY} ${SAMP_ERR_FILE} .
   else
      echo "ERROR: no sampling error correction file for this ensemble size."
      echo "ERROR: looking for ${SAMP_ERR_FILE}"
      exit -3
   endif
else
   echo "Sampling Error Correction not requested for this assimilation."
endif

#=========================================================================
# Block 3: DART INFLATION
# This stages the files that contain the inflation values.
# The inflation values change through time and should be archived.
#
# This file is only relevant if 'inflation' is turned on -
# i.e. if inf_flavor(1) /= 0 AND inf_initial_from_restart = .TRUE.
#
# filter_nml
# inf_flavor                  = 2,                       0,
# inf_initial_from_restart    = .true.,                  .false.,
# inf_in_file_name            = 'prior_inflate_ics',     'post_inflate_ics',
# inf_out_file_name           = 'prior_inflate_restart', 'post_inflate_restart',
# inf_diag_file_name          = 'prior_inflate_diag',    'post_inflate_diag',
#
# NOTICE: the archiving scripts more or less require the names of these
# files to be as listed above. When being archived, the filenames get a
# unique extension (describing the assimilation time) appended to them.
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
#    You must manually remove the cam_inflation_cookie file
#    from the RUNDIR in this case.
#    - OR -
# 2) create a cookie file called RUNDIR/cam_inflation_cookie
#    The existence of this file will cause this script to set the
#    namelist appropriately. This script will 'eat' the cookie file
#    to prevent this from happening for subsequent executions. If the
#    inflation file does not exist for them, and it needs to, this script
#    should die. The CESM_DART_config script automatically creates a cookie
#    file to support this option.
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

# its a little tricky to remove both styles of quotes from the string.

set  MYSTRING = `grep inf_in_file_name input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set  PRIOR_INF_IFNAME = $MYSTRING[2]
set  POSTE_INF_IFNAME = $MYSTRING[3]

set  MYSTRING = `grep inf_out_file_name input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set  PRIOR_INF_OFNAME = $MYSTRING[2]
set  POSTE_INF_OFNAME = $MYSTRING[3]

set  MYSTRING = `grep inf_diag_file_name input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set  PRIOR_INF_DIAG = $MYSTRING[2]
set  POSTE_INF_DIAG = $MYSTRING[3]

# IFF we want PRIOR inflation:

if ( $PRIOR_INF > 0 ) then

   if ($PRIOR_TF == false) then
      # we are not using an existing inflation file.
      echo "inf_flavor(1) = $PRIOR_INF, using namelist values."

   else if ( -e ../cam_inflation_cookie ) then
      # We want to use an existing inflation file, but this is
      # the first assimilation so there is no existing inflation
      # file. This is the signal we need to to coerce the namelist
      # to have different values for this execution ONLY.
      # Since the local namelist comes from CASEROOT each time, we're golden.

      set PRIOR_TF = FALSE

ex input.nml <<ex_end
g;inf_initial_from_restart ;s;= .*;= .${PRIOR_TF}., .${POSTE_TF}.,;
g;inf_sd_initial_from_restart ;s;= .*;= .${PRIOR_TF}., .${POSTE_TF}.,;
wq
ex_end

   else
      # Look for the output from the previous assimilation
      (ls -rt1 ../cam_${PRIOR_INF_OFNAME}.* | tail -n 1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      # If one exists, use it as input for this assimilation
      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest ${PRIOR_INF_IFNAME}
      else
         echo "ERROR: Requested PRIOR inflation but specified no incoming inflation file."
         echo "ERROR: expected something like ../cam_${PRIOR_INF_OFNAME}.YYYY-MM-DD-SSSSS"
         exit -4
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

   else if ( -e ../cam_inflation_cookie ) then
      # We want to use an existing inflation file, but this is
      # the first assimilation so there is no existing inflation
      # file. This is the signal we need to to coerce the namelist
      # to have different values for this execution ONLY.
      # Since the local namelist comes from CASEROOT each time, we're golden.

      set POSTE_TF = FALSE

ex input.nml <<ex_end
g;inf_initial_from_restart ;s;= .*;= .${PRIOR_TF}., .${POSTE_TF}.,;
g;inf_sd_initial_from_restart ;s;= .*;= .${PRIOR_TF}., .${POSTE_TF}.,;
wq
ex_end

   else
      # Look for the output from the previous assimilation
      (ls -rt1 ../cam_${POSTE_INF_OFNAME}.* | tail -n 1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      # If one exists, use it as input for this assimilation
      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest ${POSTE_INF_IFNAME}
      else
         echo "ERROR: Requested POSTERIOR inflation but specified no incoming inflation file."
         echo "ERROR: expected something like ../cam_${POSTE_INF_OFNAME}.YYYY-MM-DD-SSSSS"
         exit -5
      endif
   endif
else
   echo "Posterior Inflation       not requested for this assimilation."
endif

# Eat the cookie regardless
${REMOVE} ../cam_inflation_cookie

#=========================================================================
# Block 4: Convert N CAM initial files to DART initial condition files.
# cam_to_dart is serial code, we can do all of these at the same time
# as long as we can have unique namelists for each of them.
#
# At the end of the block, we have DART initial condition files  filter_ics.[1-N]
# that came from pointer files ../rpointer.atm_[1-N]
#
# REQUIRED DART namelist settings:
# &filter_nml:           restart_in_file_name    = 'filter_ics'
#                        restart_out_file_name   = 'filter_restart'
# &ensemble_manager_nml: single_restart_file_in  = '.false.'
# &cam_to_dart_nml:      cam_to_dart_output_file = 'dart_ics',
# &dart_to_cam_nml:      dart_to_cam_input_file  = 'dart_restart',
#                        advance_time_present    = .false.
#=========================================================================

echo "`date` -- BEGIN CAM-TO-DART"

set member = 1
while ( ${member} <= ${ensemble_size} )

   # Each member will do its job in its own directory.
   # That way, we can do N of them simultaneously -

   set MYTEMPDIR = member_${member}
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   # make sure there are no old output logs hanging around
   $REMOVE output.${member}.cam_to_dart

   set ATM_INITIAL_FILENAME = `printf ${CASE}.cam_%04d.i.${ATM_DATE_EXT}.nc  ${member}`
   set ATM_HISTORY_FILENAME = `printf ${CASE}.cam_%04d.h0.${ATM_DATE_EXT}.nc ${member}`
   set     DART_IC_FILENAME = `printf filter_ics.%04d     ${member}`
   set    DART_RESTART_FILE = `printf filter_restart.%04d ${member}`

   sed -e "s#dart_ics#../${DART_IC_FILENAME}#" \
       -e "s#dart_restart#../${DART_RESTART_FILE}#" < ../input.nml >! input.nml

   ${LINK} ../../$ATM_INITIAL_FILENAME caminput.nc
   ${LINK} ../../$ATM_HISTORY_FILENAME cam_phis.nc

   echo "starting cam_to_dart for member ${member} at "`date`
   ${EXEROOT}/cam_to_dart >! output.${member}.cam_to_dart &

   cd ..

   @ member++
end

wait

set nsuccess = `fgrep 'Finished ... at YYYY' member*/output.[0-9]*.cam_to_dart | wc -l`
if (${nsuccess} != ${ensemble_size}) then
   echo "ERROR ... DART died in 'cam_to_dart' ... ERROR"
   echo "ERROR ... DART died in 'cam_to_dart' ... ERROR"
   exit -6
endif

echo "`date` -- END CAM-TO-DART for all ${ensemble_size} members."

#=========================================================================
# Block 5: Actually run the assimilation.
# Will result in a set of files : 'filter_restart.xxxx'
#
# DART namelist settings required:
# &filter_nml:           async                   = 0,
# &filter_nml:           adv_ens_command         = "no_advance_script",
# &filter_nml:           restart_in_file_name    = 'filter_ics'
# &filter_nml:           restart_out_file_name   = 'filter_restart'
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

# CAM:static_init_model() always needs a caminput.nc and a cam_phis.nc
# for geometry information, etc.

set ATM_INITIAL_FILENAME = ${CASE}.cam_0001.i.${ATM_DATE_EXT}.nc
set ATM_HISTORY_FILENAME = ${CASE}.cam_0001.h0.${ATM_DATE_EXT}.nc

${LINK} ../$ATM_INITIAL_FILENAME caminput.nc
${LINK} ../$ATM_HISTORY_FILENAME cam_phis.nc

echo "`date` -- BEGIN FILTER"
${LAUNCHCMD} ${EXEROOT}/filter_cam || exit -7
echo "`date` -- END FILTER"

${MOVE} preassim.nc      ../cam_preassim.${ATM_DATE_EXT}.nc
${MOVE} analysis.nc      ../cam_analysis.${ATM_DATE_EXT}.nc
${MOVE} obs_seq.final    ../cam_obs_seq.${ATM_DATE_EXT}.final
${MOVE} dart_log.out     ../cam_dart_log.${ATM_DATE_EXT}.out

# Accomodate any possible inflation files
# 1) rename file to reflect current date
# 2) move to RUNDIR so the DART INFLATION BLOCK works next time and
#    that they can get archived.

foreach FILE ( ${PRIOR_INF_OFNAME} ${POSTE_INF_OFNAME} ${PRIOR_INF_DIAG} ${POSTE_INF_DIAG} )
   if ( -e ${FILE} ) then
      ${MOVE} ${FILE} ../cam_${FILE}.${ATM_DATE_EXT}
   else
      echo "No ${FILE} for ${ATM_DATE_EXT}"
   endif
end

#=========================================================================
# Block 6: Update the cam restart files ... simultaneously ...
#
# Each member will do its job in its own directory, which already exists
# and has the required input files remaining from 'Block 4'
#=========================================================================

echo "`date` -- BEGIN DART-TO-CAM"
set member = 1
while ( $member <= $ensemble_size )

   cd member_${member}

   ${REMOVE} output.${member}.dart_to_cam

   echo "starting dart_to_cam for member ${member} at "`date`
   ${EXEROOT}/dart_to_cam >! output.${member}.dart_to_cam &

   cd ..

   @ member++
end

wait

set nsuccess = `fgrep 'Finished ... at YYYY' member*/output.[0-9]*.dart_to_cam | wc -l`
if (${nsuccess} != ${ensemble_size}) then
   echo "ERROR ... DART died in 'dart_to_cam' ... ERROR"
   echo "ERROR ... DART died in 'dart_to_cam' ... ERROR"
   exit -8
endif

echo "`date` -- END DART-TO-CAM for all ${ensemble_size} members."

#=========================================================================
# Block 7: The cam files have now been updated, move them into position.
#
# As implemented, the input filenames are static in the CESM namelists.
# Since the short-term archiver creates unique directories for these,
# we must link the uniquely-named files to static names. When the short-term
# archiver 'restores' the CESM files, the links will still be valid.
#
# IMPORTANT: the DART/models/cam/shell_scripts/st_archive.sh MUST be used
# instead of the CESM st_archive.sh script.
#=========================================================================

cd ${RUNDIR}

set member = 1
while ( ${member} <= ${ensemble_size} )

   set inst_string = `printf _%04d $member`

   set ATM_INITIAL_FILENAME = ${CASE}.cam${inst_string}.i.${ATM_DATE_EXT}.nc

   ${LINK} ${ATM_INITIAL_FILENAME} cam_initial${inst_string}.nc || exit -9

   @ member++

end

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

echo "`date` -- END CAM_ASSIMILATE"

exit 0


