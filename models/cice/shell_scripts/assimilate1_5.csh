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

echo "`date` -- BEGIN ICE_ASSIMILATE"

set nonomatch       # suppress "rm" warnings if wildcard does not match anything

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following
switch ("`hostname`")
   case ys*:
         # NCAR "yellowstone"
         set      MOVE = '/bin/mv -v'
         set      COPY = '/bin/cp -v --preserve=timestamps'
         set      LINK = '/bin/ln -vs'
         set    REMOVE = '/bin/rm -rf'
         set LAUNCHCMD = mpirun.lsf
         set TASKS_PER_NODE = `echo $LSB_SUB_RES_REQ | sed -ne '/ptile/s#.*\[ptile=\([0-9][0-9]*\)]#\1#p'`
         setenv MP_DEBUG_NOTIMEOUT yes
      breaksw
   case linux_system_with_utils_in_other_dirs*:
         # example of pointing this script at a different set of basic commands
         set      MOVE = '/usr/local/bin/mv -v'
         set      COPY = '/usr/local/bin/cp -v --preserve=timestamps'
         set      LINK = '/usr/local/bin/ln -vs'
         set    REMOVE = '/usr/local/bin/rm -fr'
         set LAUNCHCMD = mpirun.lsf
      breaksw
   default:
         # NERSC "hopper"
         set      MOVE = 'mv -v'
         set      COPY = 'cp -v --preserve=timestamps'
         set      LINK = 'ln -vs'
         set    REMOVE = 'rm -fr'
         set LAUNCHCMD = "aprun -n $NTASKS"
      breaksw
endsw

# The bogus strings get replaced when CESM_DART_config is run
setenv    CASEROOT BOGUSCASEROOT
setenv BASEOBSROOT BOGUSBASEOBSDIR

#-------------------------------------------------------------------------
# Get the case-specific variables
#-------------------------------------------------------------------------

cd ${CASEROOT} || exit 1
setenv CASE           `./xmlquery CASE        -value`
setenv ensemble_size  `./xmlquery NINST_ICE   -value`
setenv ICE_COMPONENT  `./xmlquery COMP_ICE    -value`
setenv EXEROOT        `./xmlquery EXEROOT     -value`
setenv RUNDIR         `./xmlquery RUNDIR      -value`
setenv archive        `./xmlquery DOUT_S_ROOT -value`

cd ${RUNDIR}

#-------------------------------------------------------------------------
# Determine time of model state ... from file name of first member
# of the form "./${CASE}.cice_${ensemble_member}.r.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#-------------------------------------------------------------------------

set FILE = `head -n 1 rpointer.ice_0001`
set FILE = $FILE:r
set ICE_DATE_EXT = `echo $FILE:e`
set ICE_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
set ICE_YEAR     = `echo $ICE_DATE[1] | bc`
set ICE_MONTH    = `echo $ICE_DATE[2] | bc`
set ICE_DAY      = `echo $ICE_DATE[3] | bc`
set ICE_SECONDS  = `echo $ICE_DATE[4] | bc`
set ICE_HOUR     = `echo $ICE_DATE[4] / 3600 | bc`

echo "valid time of model is $ICE_YEAR $ICE_MONTH $ICE_DAY $ICE_SECONDS (seconds)"
echo "valid time of model is $ICE_YEAR $ICE_MONTH $ICE_DAY $ICE_HOUR (hours)"

#-------------------------------------------------------------------------
# Create temporary working directory for the assimilation and go there
#-------------------------------------------------------------------------

set temp_dir = assimilate_ice
echo "temp_dir is $temp_dir"

if ( -d $temp_dir ) then
   ${REMOVE} $temp_dir/*
else
   mkdir -p $temp_dir
endif
cd $temp_dir

#-----------------------------------------------------------------------------
# Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of ICE.
#-----------------------------------------------------------------------------
# Make sure the file name structure matches the obs you will be using.
# PERFECT model obs output appends .perfect to the filenames

set YYYYMM   = `printf %04d%02d                ${ICE_YEAR} ${ICE_MONTH}`
if (! -d ${BASEOBSDIR}/${YYYYMM}_CESM) then
   echo "CESM+DART requires 6 hourly obs_seq files in directories of the form YYYYMM_CESM"
   echo "The directory ${BASEOBSDIR}/${YYYYMM}_CESM is not found.  Exiting"
   exit -10
endif

set OBSFNAME = `printf obs_seq.%04d-%02d-%02d-%05d ${ICE_YEAR} ${ICE_MONTH} ${ICE_DAY} ${ICE_SECONDS}`

set OBS_FILE = ${BASEOBSDIR}/${YYYYMM}_CESM/${OBSFNAME}
echo "OBS_FILE = $OBS_FILE"

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

if (  -e   ${CASEROOT}/input.nml ) then
   ${COPY} ${CASEROOT}/input.nml .
   # TODO FIXME ... ripping out comments? not needed.
#  sed -e "/#/d;/^\!/d;/^[ ]*\!/d" ${CASEROOT}/input.nml >! input.nml  || exit 39

else
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   exit -2
endif

echo "`date` -- END COPY BLOCK"

# If possible, use the round-robin approach to deal out the tasks.
# Since the ensemble manager is not used by dart_to_ice,
# it is OK to set it here and have it used by all routines.

if ($?TASKS_PER_NODE) then
   if ($#TASKS_PER_NODE > 0) then
      ${COPY} input.nml input.nml.$$
      sed -e "s#layout.*#layout = 2#" \
          -e "s#tasks_per_node.*#tasks_per_node = $TASKS_PER_NODE#" \
          input.nml.$$ >! input.nml || exit 40
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

set  MYSTRING = `grep 'sampling_error_correction' input.nml`
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
# i.e. if inf_flavor(:) /= 0 AND inf_initial_from_restart = .TRUE.
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
#    You must manually remove the ice_inflation_cookie file
#    from the RUNDIR in this case.
#    - OR -
# 2) create a cookie file called RUNDIR/ice_inflation_cookie
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

set  MYSTRING = `grep 'inf_flavor' input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  PRIOR_INF = $MYSTRING[2]
set  POSTE_INF = $MYSTRING[3]

set  MYSTRING = `grep 'inf_initial_from_restart' input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  PRIOR_TF = `echo $MYSTRING[2] | tr '[:upper:]' '[:lower:]'`
set  POSTE_TF = `echo $MYSTRING[3] | tr '[:upper:]' '[:lower:]'`

# its a little tricky to remove both styles of quotes from the string.

set  MYSTRING = `grep 'inf_in_file_name' input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set  PRIOR_INF_IFNAME = $MYSTRING[2]
set  POSTE_INF_IFNAME = $MYSTRING[3]

set  MYSTRING = `grep 'inf_out_file_name' input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set  PRIOR_INF_OFNAME = $MYSTRING[2]
set  POSTE_INF_OFNAME = $MYSTRING[3]

set  MYSTRING = `grep 'inf_diag_file_name' input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set  PRIOR_INF_DIAG = $MYSTRING[2]
set  POSTE_INF_DIAG = $MYSTRING[3]

# IFF we want PRIOR inflation:

if ( $PRIOR_INF > 0 ) then

   if ($PRIOR_TF == false) then
      # we are not using an existing inflation file.
      echo "inf_flavor(1) = $PRIOR_INF, using namelist values."

   else if ( -e ../ice_inflation_cookie ) then
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
      (ls -rt1 ../ice_${PRIOR_INF_OFNAME}.* | tail -n 1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      # If one exists, use it as input for this assimilation
      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest ${PRIOR_INF_IFNAME}
      else
         echo "ERROR: Requested PRIOR inflation but specified no incoming inflation file."
         echo "ERROR: expected something like ../ice_${PRIOR_INF_OFNAME}.YYYY-MM-DD-SSSSS"
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

   else if ( -e ../ice_inflation_cookie ) then
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
      (ls -rt1 ../ice_${POSTE_INF_OFNAME}.* | tail -n 1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      # If one exists, use it as input for this assimilation
      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest ${POSTE_INF_IFNAME}
      else
         echo "ERROR: Requested POSTERIOR inflation but specified no incoming inflation file."
         echo "ERROR: expected something like ../ice_${POSTE_INF_OFNAME}.YYYY-MM-DD-SSSSS"
         exit -5
      endif
   endif
else
   echo "Posterior Inflation       not requested for this assimilation."
endif

# Eat the cookie regardless
${REMOVE} ../ice_inflation_cookie

#=========================================================================
# Block 4: filter has the ability to directly modify the cice restart files -
#          i.e. it creates the posterior IN-PLACE..
#          We usually want a prior estimate, so we have to save a copy of the
#          input files before we feed them to filter. 
#
# At the end of the block, we have DART initial condition files  filter_ics.[1-N]
# that icee from pointer files ../rpointer.ice_[1-N]
#
# REQUIRED DART namelist settings:
# &filter_nml:           restart_list_file       = 'cice_restarts.txt'
#                        direct_netcdf_read      = .true.
#                        direct_netcdf_write     = .true.
#                        restart_out_file_name   = 'cice_out.r'
#
## TODO FIXME ... this whole section
#
# &ensemble_manager_nml: single_restart_file_in  = '.false.'
# &ice_to_dart_nml:      ice_to_dart_output_file = 'dart_ics',
# &dart_to_ice_nml:      dart_to_ice_input_file  = 'dart_restart',
#                        advance_time_present    = .false.
#
# NOTE: when starting an OSSE by perturbing a single file, use
# &filter_nml
#   start_from_restart        = .false.,
#   output_restart            = .true.,
#   restart_in_file_name      = "filter_ics.0001",
#   restart_out_file_name     = "filter_restart", 
# &ensemble_manager_nml
#   single_restart_file_in    = .true.,
#   single_restart_file_out   = .false.,
# &model_nml
#   pert_names          = 'T'
#   pert_sd             = 1.0e-11,
#   pert_base_vals      = -888888.0d0,
#=========================================================================

echo "`date` -- BEGIN ICE-TO-DART"

# Check to make sure we are running what we are supporting

if ( $ICE_COMPONENT == 'cice' ) then

else
  # TODO error out if we
endif

cat ../rpointer.ice_* >!  cice_restarts.txt

# TODO error out if the number of files in cice_restart.txt does not match the
# ensemble size


# CP the priors since everything in cice_restart.txt will be overwritten by filter

set member = 1
while ( ${member} <= ${ensemble_size} )

   set ICE_FILENAME = `head -n $member cice_restart.txt | tail -n 1`
   
   set DART_FILENAME = `printf cice_out.r.%04d     ${member}`

   ${COPY} ${ICE_FILENAME} ${DART_FILENAME} &

   @ member++
end

wait

echo "`date` -- END ICE-TO-DART for all ${ensemble_size} members."

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

echo "`date` -- BEGIN FILTER"
${LAUNCHCMD} ${EXEROOT}/filter || exit -7
echo "`date` -- END FILTER"

# ${MOVE} Prior_Diag.nc      ../ice_Prior_Diag.${ICE_DATE_EXT}.nc
# ${MOVE} Posterior_Diag.nc  ../ice_Posterior_Diag.${ICE_DATE_EXT}.nc
# ${MOVE} obs_seq.final      ../ice_obs_seq.${ICE_DATE_EXT}.final
# ${MOVE} dart_log.out       ../ice_dart_log.${ICE_DATE_EXT}.out

# Copy obs_seq.final files to a place that won't be archived,
# so that they don't need to be retrieved from the HPSS.
if (! -d ../../Obs_seqs) mkdir ../../Obs_seqs
${COPY} ../ice_obs_seq.${ICE_DATE_EXT}.final ../../Obs_seqs &

# Accommodate any possible inflation files
# 1) rename file to reflect current date
# 2) move to RUNDIR so the DART INFLATION BLOCK works next time and
#    that they can get archived.

foreach FILE ( ${PRIOR_INF_OFNAME} ${POSTE_INF_OFNAME} ${PRIOR_INF_DIAG} ${POSTE_INF_DIAG} )
   if ( -e ${FILE} ) then
      ${MOVE} ${FILE} ../ice_${FILE}.${ICE_DATE_EXT}
   else
      echo "No ${FILE} for ${ICE_DATE_EXT}"
   endif
end

# Handle localization_diagnostics_files
set MYSTRING = `grep 'localization_diagnostics_file' input.nml`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set loc_diag = $MYSTRING[2]
if (-f $loc_diag) then
   $MOVE $loc_diag ../ice_${loc_diag}.${ICE_DATE_EXT}
endif

# Handle regression diagnostics
set MYSTRING = `grep 'reg_diagnostics_file' input.nml`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set reg_diag = $MYSTRING[2]
if (-f $reg_diag) then
   $MOVE $reg_diag ../ice_${reg_diag}.${ICE_DATE_EXT}
endif

# 
#=========================================================================
# Block 6: Update the ice restart files ... simultaneously ...
#
# Each member will do its job in its own directory.
# Block 7: The ice files have now been updated, move them into position.
#=========================================================================

echo "`date` -- BEGIN DART-TO-ICE"
set member = 1
while ( $member <= $ensemble_size )

   cd member_${member}

   ${REMOVE} output.${member}.dart_to_ice

   echo "starting dart_to_ice for member ${member} at "`date`
   ${EXEROOT}/dart_to_ice >! output.${member}.dart_to_ice &

   set inst_string = `printf _%04d $member`

   set ICE_INITIAL_FILENAME = ${CASE}.ice${inst_string}.i.${ICE_DATE_EXT}.nc

   ${LINK} ${ICE_INITIAL_FILENAME} dart_restart.nc || exit -9

   cd ..

   @ member++
end

wait

set nsuccess = `fgrep 'Finished ... at YYYY' member*/output.[0-9]*.dart_to_ice | wc -l`
if (${nsuccess} != ${ensemble_size}) then
   echo "ERROR ... DART died in 'dart_to_ice' ... ERROR"
   echo "ERROR ... DART died in 'dart_to_ice' ... ERROR"
   exit -8
endif

echo "`date` -- END DART-TO-ICE for all ${ensemble_size} members."

#

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

echo "`date` -- END ICE_ASSIMILATE"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

