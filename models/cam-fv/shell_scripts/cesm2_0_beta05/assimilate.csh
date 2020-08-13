#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# This script is designed to interface cesm2_0_beta05 or later
# and $dart/rma_trunk v11###.

# See 'RMA' for places where there are questions about RMA,
# especially file naming.



#=========================================================================
# Block 0: Set command environment
#=========================================================================
# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

echo "`date` -- BEGIN CAM_ASSIMILATE"
pwd

set nonomatch       # suppress "rm" warnings if wildcard does not match anything

setenv CASEROOT $1
# Python uses C indexing on loops; cycle = [0,....,$DATA_ASSIMILATION_CYCLES - 1]
# "Fix" that here, so the rest of the script isn't confusing.
@ cycle = $2 + 1

# In CESM1_4 xmlquery must be executed in $CASEROOT.
cd ${CASEROOT}
setenv CASE           $CASEROOT:t
setenv ensemble_size  `./xmlquery NINST_ATM   --value`
setenv CAM_DYCORE     `./xmlquery CAM_DYCORE  --value`
setenv EXEROOT        `./xmlquery EXEROOT     --value`
setenv RUNDIR         `./xmlquery RUNDIR      --value`
setenv archive        `./xmlquery DOUT_S_ROOT --value`
setenv TOTALPES       `./xmlquery TOTALPES    --value`
setenv DATA_ASSIMILATION_CYCLES        `./xmlquery DATA_ASSIMILATION_CYCLES --value`
cd $RUNDIR

setenv save_all_inf TRUE
# A switch to signal how often to save the stages' ensemble members: NONE, RESTART_TIMES, ALL
setenv save_stages RESTART_TIMES

set BASEOBSDIR = BOGUSBASEOBSDIR

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following
switch ($HOSTNAME)
   case ys*:
      # NCAR "yellowstone"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'
      set TASKS_PER_NODE = `echo $LSB_SUB_RES_REQ | sed -ne '/ptile/s#.*\[ptile=\([0-9][0-9]*\)]#\1#p'`
      setenv MP_DEBUG_NOTIMEOUT yes
      set  LAUNCHCMD = mpirun.lsf
   case ch*:
      # NCAR "cheyenne"   -v removed because it doesn't work on (my) cheyenne.
      set   MOVE = 'mv -f'
      set   COPY = 'cp -f --preserve=timestamps'
      set   LINK = 'ln -fs'
      set REMOVE = 'rm -fr'
      # echo "Trying to set TASKS_PER_NODE using PBS_NUM_PPN $PBS_NUM_PPN"
      # Unavailable for some reason:  set TASKS_PER_NODE = $PBS_NUM_PPN 
      set TASKS_PER_NODE = 36
      setenv MP_DEBUG_NOTIMEOUT yes
      set  LAUNCHCMD = mpiexec_mpt
   breaksw

   case linux_system_with_utils_in_other_dirs*:
      # example of pointing this script at a different set of basic commands
      set   MOVE = '/usr/local/bin/mv -fv'
      set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
      set   LINK = '/usr/local/bin/ln -fvs'
      set REMOVE = '/usr/local/bin/rm -fr'
      set LAUNCHCMD  = mpirun.lsf
   breaksw

   default:
      # NERSC "hopper"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'
      set LAUNCHCMD  = "aprun -n $TOTALPES"

   breaksw
endsw

#=========================================================================
# Block 1: Populate a run-time directory with the input needed to run DART.
#=========================================================================

echo "`date` -- BEGIN COPY BLOCK"

if (  -e   ${CASEROOT}/input.nml ) then
   # ${COPY} ${CASEROOT}/input.nml .
   # Put a pared down copy (no comments) of input.nml in this assimilate_cam directory.
   sed -e "/#/d;/^\!/d;/^[ ]*\!/d" ${CASEROOT}/input.nml >! input.nml  || exit 39
else
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   exit -2
endif

echo "`date` -- END COPY BLOCK"

# If possible, use the round-robin approach to deal out the tasks.

if ($?TASKS_PER_NODE) then
   if ($#TASKS_PER_NODE > 0) then
      ${MOVE} input.nml input.nml.$$
      sed -e "s#layout.*#layout = 2#" \
          -e "s#tasks_per_node.*#tasks_per_node = $TASKS_PER_NODE#" \
          input.nml.$$ >! input.nml || exit 40
      ${REMOVE} input.nml.$$
   endif
endif

#=========================================================================
# Block 2: Identify requested output stages
#=========================================================================
# 
set MYSTRING = `grep stages_to_write input.nml`
set MYSTRING = (`echo $MYSTRING | sed -e "s#[=,'\.]# #g"`)
set STAGE_input     = TRUE
set STAGE_forecast  = FALSE
set STAGE_preassim  = FALSE
set STAGE_postassim = FALSE
set STAGE_analysis  = FALSE
set STAGE_output    = TRUE
# Assemble list of stages to write out, which are not the 'output' stage.
# Stages 'input' and 'output' are always written out.
set stages_except_output = "{input"
set stage = 2
while ($stage <= $#MYSTRING) 
   if ($MYSTRING[$stage] == 'forecast')  then
      set STAGE_forecast = TRUE
      set stages_except_output = "${stages_except_output},forecast"
   endif
   if ($MYSTRING[$stage] == 'preassim')  then
      set STAGE_preassim = TRUE
      set stages_except_output = "${stages_except_output},preassim"
   endif
   if ($MYSTRING[$stage] == 'postassim') then
      set STAGE_postassim = TRUE
      set stages_except_output = "${stages_except_output},postassim"
   endif
   if ($MYSTRING[$stage] == 'analysis')  then
      set STAGE_analysis = TRUE
      set stages_except_output = "${stages_except_output},analysis"
   endif
   @ stage++
end
set stages_all = "${stages_except_output},output}"
set stages_except_output = "${stages_except_output}}"
# Checking
echo "stages_except_output = $stages_except_output"
echo "stages_all = $stages_all"


#=========================================================================
# Block 3: Preliminary clean up, which can run in the background.
#=========================================================================
# CESM2_0's new archiver has a mechanism for removing restart file sets,
# which we don't need, but it runs only after the (multicycle) job finishes.  
# We'd like to remove unneeded restarts as the job progresses, allowing more
# cycles to run before needing to stop to archive data.  So clean them out of
# RUNDIR, and st_archive will never have to deal with them.
#-------------------------------------------------------------------------

# Move any hidden restart sets back into the run directory so they can be used or purged.
ls -d ../Hide*
if ($status == 0) then
   echo 'Moving files from ../Hide* to $rundir'
   $MOVE ../Hide*/* .
   rmdir ../Hide*
endif

# Make a place to store inflation restarts to protect from purging until
# st_archive can make a home for them.
if (! -d $archive/esp/rest && $save_all_inf =~ TRUE) mkdir -p $archive/esp/rest

# Cwd is currently RUNDIR.
set log_list = `ls -t cesm.log.*`
echo "log_list is $log_list"

# For safety, leave the most recent *2* restart sets in place.
# Prevents catastrophe if the last restart set is partially written before a crash.
# Add 1 more because the restart set used to start this will be counted:
# there will be 3 restarts when there are only 2 cesm.log files,
# which caused all the files to be deleted.
if ($#log_list >= 3) then

   # List of potential restart sets to remove
   set re_list = `ls -t *cpl.r.*`
   if ($#re_list < 3) then
      echo "Not enough restart sets, even though there are $#log_list cesm.log files"
      exit 14
   endif
   # Member restarts to remove
   set rm_date = `echo $re_list[3] | sed -e "s/-/ /g;s/\./ /g;"`
   @ day = $#rm_date - 2
   @ sec = $#rm_date - 1
   # Fix the misinterpretation of day_o_month = 08 as an octal number
   # by stripping any leading '0's.
   set day_o_month = `echo $rm_date[$day] | bc`
   set sec_o_day   = $rm_date[$sec]
   set day_time = ${day_o_month}-${sec_o_day}

   # Identify log files to be removed or moved.
   # [3] means the 3rd oldest restart set is being (re)moved.
   set rm_log = `echo $log_list[3] | sed -e "s/\./ /g;"`
   set rm_slot = $#rm_log
   if ($rm_log[$#rm_log] == 'gz') @ rm_slot--
   echo '$rm_log['$rm_slot']='$rm_log[$rm_slot]

   if ( $sec_o_day !~ '00000' || \
       ($sec_o_day =~ '00000' && $day_o_month % BOGUS_save_every_Mth != 0) ) then
      echo "Removing unneeded restart file set from RUNDIR: "
      echo "    ${CASE}"'*.{r,rs,rh0,h0,i}.*'${day_time}
      # Optionally save inflation restarts, even if it's not a 'save restart' time.
      if ($save_all_inf =~ TRUE) $MOVE ${CASE}*inf*${day_time}*  $archive/esp/rest

      # Remove intermediate member restarts (but not DART means, sd, obs_seq, inflation restarts output)
      # Note that *cpl.ha.* is retained, and any h#, #>0.
      #        $CASE                          DD          -SSSSS
      $REMOVE  ${CASE}*.{r,rs,rh0,h0}.*${day_time}* &
      # Handle .i. separately to avoid sweeping up .input_{mean,sd,...} files.
      $REMOVE ${CASE}*.i.[^io]*${day_time}*  &
      if ($save_stages =~ NONE || $save_stages =~ RESTART_TIMES) then
         $REMOVE  ${CASE}.*[0-9].${stages_except_output}*${day_time}* &
      endif
   else
      # Optionally COPY inflation restarts to the same place as the other inflation restarts.
      if ($save_all_inf =~ TRUE) $COPY ${CASE}*inf*${day_time}*  $archive/esp/rest

      # Optionally REMOVE stages' ensemble members (not means and sds).
      if ($save_stages =~ NONE ) then
         $REMOVE  ${CASE}.*[0-9].${stages_except_output}*${day_time}* &
      endif

      # Save the restart set to archive/rest/$datename, where it will be safe
      # from removes of $component/rest and ${case}.locked/archive.
      set save_date = `echo $re_list[3] | sed -e "s/\./ /g;"`
      @ piece = $#save_date - 1
      set save_root = $archive/rest/$save_date[$piece]
      if (! -d $save_root) then
         mkdir -p $save_root
         echo "Moving restart to $save_root"
         $MOVE ${CASE}*.{r,rs,rh0,h0}.*${day_time}*  $save_root &
         # Handle .i. separately to avoid sweeping up .{in,out}put_{mean,sd,...} files.
         $MOVE ${CASE}*.i.[^io]*${day_time}* $save_root &
         $COPY            *inf*${day_time}*  $save_root &
         # Save a few log files
         echo "Moving logs to $save_root"
         $MOVE *0001*$rm_log[$rm_slot]*      $save_root 

      else
         echo "$save_root already exists.  Not archiving restart there."
      endif
   endif
   # Remove log files: *YYMMDD-HHMMSS*.  Except not da.log files
   $REMOVE  [^d]*$rm_log[$rm_slot]*  &

   # I'd like to remove the CAM .r. files, since we always use the .i. files to do a hybrid start,
   # but apparently CESM needs them to be there, even though it doesn't read fields from them.
   # $REMOVE  ${CASE}.cam*.r.*${day_time}.nc &


endif


#=========================================================================
# Block 4: Determine time of model state 
#=========================================================================
# ... from file name of first member
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

#-----------------------------------------------------------------------------
# Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of CAM.
#-----------------------------------------------------------------------------
# Make sure the file name structure matches the obs you will be using.
# PERFECT model obs output appends .perfect to the filenames

set YYYYMM   = `printf %04d%02d                ${ATM_YEAR} ${ATM_MONTH}`
if (! -d ${BASEOBSDIR}/${YYYYMM}_6H_CESM) then
   echo "CESM+DART requires 6 hourly obs_seq files in directories of the form YYYYMM_6H_CESM"
   echo "The directory ${BASEOBSDIR}/${YYYYMM}_6H_CESM is not found.  Exiting"
   exit -10
endif

set OBSFNAME = `printf obs_seq.%04d-%02d-%02d-%05d ${ATM_YEAR} ${ATM_MONTH} ${ATM_DAY} ${ATM_SECONDS}`

set OBS_FILE = ${BASEOBSDIR}/${YYYYMM}_6H_CESM/${OBSFNAME}
echo "OBS_FILE = $OBS_FILE"

# RMA new obs_seq names?
if (  -e   ${OBS_FILE} ) then
   ${LINK} ${OBS_FILE} obs_seq.out
else
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... no observation file $OBS_FILE"
   exit -1
endif

#=========================================================================
# Block 5: Stage the files needed for SAMPLING ERROR CORRECTION
#
# The sampling error correction is a lookup table.
# The tables were originally in the DART distribution, but should
# have been staged to $CASEROOT at setup time.
# RMA:
# There's a single file which has a table for each ensemble size 2,...,100
# It is only needed if
# input.nml:&assim_tools_nml:sampling_error_correction = .true.,
# which is the default.
#=========================================================================

set  MYSTRING = `grep sampling_error_correction input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set SECSTRING = `echo $MYSTRING[2] | tr '[:upper:]' '[:lower:]'`

# RMA change file name (and location?).
if ( $SECSTRING == true ) then
   set SAMP_ERR_FILE = ${CASEROOT}/final_full.${ensemble_size}
   # eventually: 
   # set SAMP_ERR_FILE = ${CASEROOT}/sampling_error_correction_table.nc
   if (  -e   ${SAMP_ERR_FILE} ) then
      ${COPY} ${SAMP_ERR_FILE} .
   else
      echo "ERROR: no sampling error correction file."
      echo "ERROR: looking for ${SAMP_ERR_FILE}"
      exit -3
   endif
else
   echo "Sampling Error Correction not requested for this assimilation."
endif

#=========================================================================
# Block 6: DART INFLATION
# This stages the files that contain the inflation values.
# The inflation values change through time and should be archived.
#
# This file is only relevant if 'inflation' is turned on -
# i.e. if inf_flavor(:) /= 0 AND inf_initial_from_restart = .TRUE.
#
# filter_nml
# inf_flavor                  = 2,                       0,
# inf_initial_from_restart    = .true.,                  .false.,
# inf_diag_file_name          = 'output_priorobsinfl',   'output_postobsinfl',
#
# NOTICE: the archiving scripts will look for inflation files which have
# names in the CESM style, but with new file types (where .i., .r.,... go).
#
# The inflation file is essentially a duplicate of the DART model state ...
# For the purpose of this script, they are the output of a previous assimilation,
# so they should be named something like input_priorinf.YYYY-MM-DD-SSSSS
#
# RMA; did Tim and I figure out that we can always have an inflation file(s)
#      ready before the first assimilation?  -> no cookie?
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
#    should die. The DART_config script automatically creates a cookie
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

if ($PRIOR_TF == FALSE && $STAGE_forecast == TRUE && $STAGE_preassim == TRUE) then
   echo " "
   echo "WARNING ! ! prior inflation is OFF, "
   echo "            but output at stages 'forecast' and 'preassim' is requested."
   echo "            This is redundant output."
   echo " "
endif
if ($POSTE_TF == FALSE && $STAGE_postassim == TRUE && $STAGE_analysis == TRUE) then
   echo " "
   echo "WARNING ! ! posterior inflation is OFF, "
   echo "            but output at stages 'postassim' and 'analysis' is requested."
   echo "            This is redundant output."
   echo " "
endif
# its a little tricky to remove both styles of quotes from the string.

# IFF we want PRIOR inflation:

if ( $PRIOR_INF > 0 ) then

   if ($PRIOR_TF == false) then
      # we are not using an existing inflation file.
      echo "inf_flavor(1) = $PRIOR_INF, using namelist values."

   else if ( -e cam_inflation_cookie ) then
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
      # RMA; file 'type' output_priorinf is hardwired according to DART2.0 naming convention.
      # Yes, we want the 'output' or 'postassim' version of the prior inflation,
      # because the 'preassim' version has not been updated by this cycle's assimilation.

      # If inflation files exists, use them as input for this assimilation
      # Must be separate commands because the 'order' that means and sds 
      # are finished being written out varies from cycle to cycle.
      # Leaving $CASE off of these ls allows it to find inflation files from ref_case.
      (ls -rt1 *.cam.e.output_priorinf_mean* | tail -n 1 >! latestfile) > & /dev/null
      (ls -rt1 *.cam.e.output_priorinf_sd*   | tail -n 1 >> latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`
      if ( $nfiles > 0 ) then
         set latest_mean = `head -n1 latestfile`
         set latest_sd   = `tail -n1 latestfile`
         # These file names are in use as of (r10786 2016-12-14).
         # filter/filter_mod.f90:
         #    call set_file_metadata(file_info, PRIOR_INF_MEAN, stage, 'priorinf_mean', 'prior inflation mean')
         # Need to COPY instead of link because of short-term archiver and disk management.
         ${COPY} $latest_mean input_priorinf_mean.nc
         ${COPY} $latest_sd   input_priorinf_sd.nc
      else
         echo "ERROR: Requested PRIOR inflation restart, "
         echo "       but there are no ouput_priorinf_* files for the next cycle to use."
         ls -l *inf*
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

   else if ( -e cam_inflation_cookie ) then
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
      # Look for the output from the previous assimilation.
      # (The only stage after posterior inflation.)
      (ls -rt1 *.cam.e.output_postinf_mean* | tail -n 1 >! latestfile) > & /dev/null
      (ls -rt1 *.cam.e.output_postinf_sd*   | tail -n 1 >> latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      # If one exists, use it as input for this assimilation
      if ( $nfiles > 0 ) then
         set latest_mean = `head -n1 latestfile`
         set latest_sd   = `tail -n1 latestfile`
         ${LINK} $latest_mean input_postinf_mean.nc
         ${LINK} $latest_sd   input_postinf_sd.nc
      else
         # Stage 'output' will always be written, so the postinf_* files will be written.
         echo "ERROR: Requested POSTERIOR inflation restart, " 
         echo "       but there are no output_postinf_* files for this cycle to use."
         exit -5
      endif
   endif
else
   echo "Posterior Inflation       not requested for this assimilation."
endif

# Eat the cookie regardless
${REMOVE} cam_inflation_cookie

#=========================================================================
# Block 7: Actually run the assimilation. 
# WARNING: this version may overwrite the input, depending on {in,out}put_state_file_list.
#
# DART namelist settings required (or highly recommended):
# &filter_nml:           async                   = 0,
# &filter_nml:           adv_ens_command         = "no_CESM_advance_script",
# &filter_nml:           stages_to_write         = 'preassim','analysis'
# &filter_nml:           distributed_state       = .true.
# &filter_nml:           input_state_file_list   = 'cam_init_files'
# &filter_nml:           output_state_file_list  = 'cam_init_files',
# &filter_nml:           single_file_in          = .false.
# &filter_nml:           single_file_out         = .false.
# &filter_nml:           obs_sequence_in_name    = 'obs_seq.out'
# &filter_nml:           obs_sequence_out_name   = 'obs_seq.final'
# &filter_nml:           init_time_days          = -1,
# &filter_nml:           init_time_seconds       = -1,
# &filter_nml:           first_obs_days          = -1,
# &filter_nml:           first_obs_seconds       = -1,
# &filter_nml:           last_obs_days           = -1,
# &filter_nml:           last_obs_seconds        = -1,
#
#=========================================================================

# CAM:static_init_model() always needs a caminput.nc and a cam_phis.nc
# for geometry information, etc.

set ATM_INITIAL_FILENAME = ${CASE}.cam_0001.i.${ATM_DATE_EXT}.nc
set ATM_HISTORY_FILENAME = ${CASE}.cam_0001.h0.${ATM_DATE_EXT}.nc
${LINK} $ATM_INITIAL_FILENAME caminput.nc
${LINK} $ATM_HISTORY_FILENAME cam_phis.nc

# RMA
# Put the names of the CAM initial files, from which filter will get the model state(s),
# into a text file, whose name is provided to filter in filter_nml:input_state_file_list.
# If filter will create an ensemble from a single state, it's fine (and convenient)
# to put the whole list of files in input_state_file_list.  Filter will just use the first.
# NOTE: if the files in input_state_file_list are CESM format (all vars and all meta data),
#       then they may end up with a different structure than the preassim and postassim stage
#       output written by filter.  This can be prevented (at the cost of more disk space)
#       by copying the CESM format input files into the names filter will use, e.g.
#       $case.cam_0001.i.$date.nc --> preassim_member_0001.nc.  
#       Filter will replace the state variables with updated versions, but leave the other
#       variables and all metadata unchanged.
set line = `grep input_state_file_list input.nml | sed -e "s#[=,'\.]# #g"`
echo "$line"
set input_file_list = $line[2]

ls -1 ${CASE}.cam_[0-9][0-9][0-9][0-9].i.${ATM_DATE_EXT}.nc >! $input_file_list

# If the file names in $output_state_file_list = names in $input_state_file_list,
# then the restart file contents will be overwritten with the states updated by DART.
# This is the behavior from DART1.0.
set line = `grep output_state_file_list input.nml | sed -e "s#[=,'\.]# #g"` 
set output_file_list = $line[2]

if ($input_file_list != $output_file_list) then
   # Replace the filetype with something showing it's after the assimilation.
   # This file_type name is consistent with RMA's file naming convention.
   sed -e "s#\.i\.#.i_output.#" $input_file_list >! $output_file_list
   # If we want the soon-to-be-updated state to be put into an initial file,
   # the initial file must exist first.  Copy the existing i files into a new ones.
   set i_input  = (`cat $input_file_list`)
   set i_output = (`cat $output_file_list`)
   set f = 1
   while ($f <= $#i_input)
      ${COPY} $i_input[$f] $i_output[$f] &
      @ f++
   end
   wait
endif

echo "`date` -- BEGIN FILTER"
${LAUNCHCMD} ${EXEROOT}/filter || exit -7
echo "`date` -- END FILTER"

#========================================================================
# Block 8: Tag the output with the valid time of the model state
#=========================================================================

# Tag the state output
# Rename filter's output files with CESM format:  
#    ${case}.cam{_$inst}.${filetype}{.dart_file}.${date}.nc 
#    These should all be named with $scomp = "cam" to distinguish
#    them from the same output from other components in multi-component assims.
#    The file type "i" or "e.$dart_file" will distinguish CAM files from DART.
# If output_state_file_list is filled with custom (CESM) filenames,
# then 'output' ensemble members will not appear with filter's default,
# hard-wired names.  But file types output_{mean,sd} will appear and be
# renamed here.
foreach FILE (`ls ${stages_all}_{mean,sd}*.nc`)
   set parts = `echo $FILE | sed -e "s#\.# #g"`

   set type = "e"
   echo $FILE | grep "put"
   if ($status == 0) set type = "i"
   mv $FILE ${CASE}.cam.${type}.$parts[1].${ATM_DATE_EXT}.nc
   # Checking
   echo "moving $FILE ${CASE}.cam.${type}.$parts[1].${ATM_DATE_EXT}.nc"
end

# Files with instance/member number handled differently from others.
foreach FILE (`ls ${stages_all}_member_*.nc`)
   # split off the .nc
   set parts = `echo $FILE | sed -e "s#\.# #g"`
   # separate the pieces of the remainder
   set list = `echo $parts[1]  | sed -e "s#_# #g"`
   # grab all but the trailing 'member' and #### parts.
   @ last = $#list - 2
   # and join them back together
   set dart_file = `echo $list[1-$last] | sed -e "s# #_#g"`

   set type = "e"
   echo $FILE | grep "put"
   if ($status == 0) set type = "i"
   mv $FILE ${CASE}.cam_$list[$#list].${type}.${dart_file}.${ATM_DATE_EXT}.nc
   # Checking
   echo "moving $FILE ${CASE}.cam_$list[$#list].${type}.${dart_file}.${ATM_DATE_EXT}.nc"
end

# Tag the observation file and run-time output

${MOVE} obs_seq.final             ${CASE}.cam.e.obs_seq_final.${ATM_DATE_EXT}
${MOVE} dart_log.out              cam_dart_log.${ATM_DATE_EXT}.out

# Copy obs_seq.final files to a place that won't be archived,
# so that they don't need to be retrieved from the HPSS.
if (! -d ../Obs_seqs) mkdir ../Obs_seqs
cp ${CASE}.cam.e.obs_seq_final.${ATM_DATE_EXT}  ../Obs_seqs &

# Tag the inflation files

# Accommodate any possible inflation files
# 1) rename file to reflect current date
# 2) move to RUNDIR so the DART INFLATION BLOCK works next time and
#    that they can get archived.

foreach FILE ( `ls ${stages_all}_{prior,post}inf_*`)
   set parts = `echo $FILE | sed -e "s#\.# #g"`
   ${MOVE} $FILE $CASE.cam.e.$parts[1].${ATM_DATE_EXT}.nc
   echo "Moved ${FILE} to $CASE.cam.e.$parts[1].${ATM_DATE_EXT}.nc"
end

# RMA; do these files have new names?
# Handle localization_diagnostics_files
set MYSTRING = `grep 'localization_diagnostics_file' input.nml`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set loc_diag = $MYSTRING[2]
if (-f $loc_diag) then
   $MOVE $loc_diag cam_${loc_diag}.${ATM_DATE_EXT}
endif

# Handle regression diagnostics
set MYSTRING = `grep 'reg_diagnostics_file' input.nml`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set reg_diag = $MYSTRING[2]
if (-f $reg_diag) then
   $MOVE $reg_diag cam_${reg_diag}.${ATM_DATE_EXT}
endif


# RMA
# Then this script will need to feed the files in output_restart_list_file
# to the next model advance.  
set line = `grep 0001 $output_file_list | sed -e "s#[\.]# #g"` 
set l = 1
while ($l < $#line)
   if ($line[$l] =~ cam_0001) then
      @ l++
      set file_type = $line[$l]
      break
   endif
   @ l++
end

set member = 1
while ( ${member} <= ${ensemble_size} )

   set inst_string = `printf _%04d $member`
   set ATM_INITIAL_FILENAME = ${CASE}.cam${inst_string}.${file_type}.${ATM_DATE_EXT}.nc

   ${LINK} ${ATM_INITIAL_FILENAME} cam_initial${inst_string}.nc || exit -9

   @ member++

end

if ($cycle == $DATA_ASSIMILATION_CYCLES) then
   if ($#log_list >= 3) then
      # During the last cycle, hide the 2nd newest restart set 
      # so that it's not archived, but is available for debugging.
      # This is assuming that DATA_ASSIMILATION_CYCLES has not been changed in env_run.xml
      # since the start (not submission) of this job.
      # (Requested by Karspeck for coupled assims, which need to keep 4 atmospheric
      #  cycles for each ocean cycle.)
      set hide_date = `echo $re_list[2] | sed -e "s/-/ /g;s/\./ /g;"`
      @ day = $#hide_date - 2
      @ sec = $#hide_date - 1
      set day_o_month = `echo $hide_date[$day] | bc`
      set sec_o_day   = $hide_date[$sec]
      set day_time = ${day_o_month}-${sec_o_day}
      set hidedir = ../Hide_${day_time}
      mkdir $hidedir
      # Optionally COPY the 2nd to last inflation restarts to the archive directory.
      # Try copying 2nd-to-last and last (to protect last from st_archive putting them in exp/hist)
       
      if ($save_all_inf =~ TRUE) $MOVE \
            ${CASE}*.${stages_except_output}*inf*  $archive/esp/rest
      $COPY ${CASE}*.output*inf*${day_time}*       $archive/esp/rest
  
      # input*inf at this last cycle will be linked to the output of the previous cycle,
      # which is being hidden here.
      # output*inf must be copied because it needs to be in rundir when st_archive runs
      # to save the results of the following assim (number $DATA_ASSIMILATION_CYCLES).
      $MOVE           ${CASE}*inf*${day_time}*    $hidedir
      $COPY  $hidedir/${CASE}*.output*inf*${day_time}*    .
 
      $MOVE ${CASE}*.[rhi]*${day_time}*    $hidedir
      # Move DART files (but not output inflation) back here for possible archiving.
      $MOVE ${hidedir}/*[^f]_{mean,sd}* .
      $MOVE ${hidedir}/*.${stages_except_output}*inf_{mean,sd}* .
      # $MOVE ${hidedir}/*.[ifpa]*inf_{mean,sd}* .

      # Move log files: *YYMMDD-HHMMSS.  [2] means the 2nd newest restart set is being moved.
      set rm_log = `echo $log_list[2] | sed -e "s/\./ /g;"`
      # -1 skips the gz at the end of the names.
      set rm_slot = $#rm_log
      if ($rm_log[$#rm_log] =~ gz) @ rm_slot--
      echo "Moving log files into $hidedir;"' $rm_log['$rm_slot']='$rm_log[$rm_slot]
      $MOVE  *$rm_log[$rm_slot]*  $hidedir
   endif

   # DEBUG st_archive by making a shadow copy of this directory.
   mkdir ../run_shadow
   foreach f (`ls`)
      ls -l $f > ../run_shadow/$f
   end
endif   

echo "`date` -- END CAM_ASSIMILATE"

# Be sure that the removal of unneeded restart sets and copy of obs_seq.final are finished.
wait

exit 0



