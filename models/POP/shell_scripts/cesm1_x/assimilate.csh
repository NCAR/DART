#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

echo "`date` -- BEGIN POP_ASSIMILATE"

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

      set BASEOBSDIR = /glade/proj3/image/Observations/WOD09
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

      set BASEOBSDIR = /glade/p/image/Observations/WOD09
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

set ensemble_size = ${NINST_OCN}

#-------------------------------------------------------------------------
# Determine time of model state ... from file name of first member
# of the form "./${CASE}.pop_${ensemble_member}.r.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#-------------------------------------------------------------------------

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

#-----------------------------------------------------------------------------
# Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of POP.
#-----------------------------------------------------------------------------

set YYYYMMDD = `printf %04d%02d%02d ${OCN_YEAR} ${OCN_MONTH} ${OCN_DAY}`
set YYYYMM   = `printf %04d%02d     ${OCN_YEAR} ${OCN_MONTH}`
set OBSFNAME = obs_seq.0Z.${YYYYMMDD}
set OBS_FILE = ${BASEOBSDIR}/${YYYYMM}/${OBSFNAME}

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
else
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   exit -2
endif

echo "`date` -- END COPY BLOCK"

# If possible, use the round-robin approach to deal out the tasks.

if ($?TASKS_PER_NODE) then
   if ($#TASKS_PER_NODE > 0) then
      ${COPY} input.nml input.nml.$$
      sed -e "s#layout.*#layout = 2#" \
          -e "s#tasks_per_node.*#tasks_per_node = $TASKS_PER_NODE#" input.nml.$$ >! input.nml
      ${REMOVE} input.nml.$$
   endif
endif

#=========================================================================
# Block 2: DART INFLATION
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
         exit -4
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

# Eat the cookie regardless
${REMOVE} pop_inflation_cookie latestfile

#=========================================================================
# Block 3: Actually run the assimilation. 
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

${REMOVE} restarts_in.txt restarts_out.txt

foreach FILE ( rpointer.ocn_????.restart )
   head -n 1 ${FILE} >> restarts_in.txt
end

# WARNING: this is the part where the files just get overwritten
${COPY} restarts_in.txt restarts_out.txt

# POP always needs a pop_in and a pop.r.nc to start.
# Lots of ways to get the filename

set OCN_RESTART_FILENAME = `head -n 1 rpointer.ocn_0001.restart`

${LINK} $OCN_RESTART_FILENAME pop.r.nc
${LINK} pop2_in_0001          pop_in

# On yellowstone, you can explore task layouts with the following:
if ( $?LSB_PJL_TASK_GEOMETRY ) then
   setenv ORIGINAL_LAYOUT "${LSB_PJL_TASK_GEOMETRY}"

   # setenv GEOMETRY_32_1NODE \
   #    "{(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31)}";
   # setenv LSB_PJL_TASK_GEOMETRY "${GEOMETRY_32_1NODE}"
endif

echo "`date` -- BEGIN FILTER"
${LAUNCHCMD} ${EXEROOT}/filter || exit -7
echo "`date` -- END FILTER"

if ( $?LSB_PJL_TASK_GEOMETRY ) then
   setenv LSB_PJL_TASK_GEOMETRY "${ORIGINAL_LAYOUT}"
endif

#========================================================================
# Block 4: Tag the output with the valid time of the model state
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

echo "`date` -- END POP_ASSIMILATE"

exit 0


