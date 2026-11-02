#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#

# ------------------------------------------------------------------------------
# Purpose: assimilate with a CESM ensemble and perform advanced archiving
#          and compression in support of multiple assimilation cycles in a
#          single CESM job.
#
# The (resulting) assimilate.csh script is called by CESM with two arguments:
# 1) the CASEROOT, and
# 2) the assimilation cycle number in this CESM job
# ------------------------------------------------------------------------------
# This template is lightly modified by the setup scripts to be appropriate
# for specific hardware and other configurations. The modified result is
# then given execute permission and is appropriate to use for an assimilation.
# All of this is automatically performed by the DART-supplied setup scripts.
#
# Tag DART's state output with names using CESM's convention:
#    ${case}.${scomp}[_$inst].${filetype}[.$dart_file].${date}.nc
#    These should all be named with $scomp = "cam" to distinguish
#    them from the same output from other components in multi-component assims.
#
# This script also has logic in it to manage disk space in a way that allows
# for more assimilation cycles to be performed before archiving without losing
# critical restart capability. The same logic is also useful for assimilations
# that may require multiple timesteps to be available.
#
# As a specific example, consider the case when 3 assimilation cycles have been
# performed: 6Z, 12Z, 18Z.
# If we want to keep a restart set and a backup
# restart set, we only need the 18Z and 12Z, so the 6Z set can be removed.
# Let's also say that its the last cycle of job - which automatically kicks off
# the short-term archiver. If we did 'nothing', the 12Z and 18Z get archived
# and the 18Z gets restaged

# machine-specific dereferencing

if ($?SLURM_JOB_ID) then

   # SLURM environment variables:
   # env | grep SLURM | sort

   setenv ORIGINALDIR $SLURM_SUBMIT_DIR
   setenv     JOBNAME $SLURM_JOB_NAME
   setenv       JOBID $SLURM_JOBID
   setenv     MYQUEUE $SLURM_JOB_PARTITION
   setenv   NODENAMES $SLURM_NODELIST
   setenv LAUNCHCMD "mpirun -np $SLURM_NTASKS -bind-to core"
#  untested method for determining NUMTASKS_PERNODE with SLURM
#  set ANY_OLD_NODE = `head -n 1 $SLURM_NODELIST`
#  setenv NUMTASKS_PERNODE `grep $ANY_OLD_NODE $SLURM_NODELIST | wc -l`

else if ($?PBS_NODEFILE) then

   # PBS environment variables:
   # env | grep PBS | sort

   setenv ORIGINALDIR $PBS_O_WORKDIR
   setenv     JOBNAME $PBS_JOBNAME
   setenv       JOBID $PBS_JOBID
   setenv     MYQUEUE $PBS_O_QUEUE
   setenv     NUMCPUS $NCPUS
   setenv    NUMTASKS `cat  $PBS_NODEFILE | wc -l`
   setenv    NUMNODES `uniq $PBS_NODEFILE | wc -l`
   set ANY_OLD_NODE = `head -n 1 $PBS_NODEFILE`
   setenv    NUMTASKS_PERNODE `grep $ANY_OLD_NODE $PBS_NODEFILE | wc -l`
   setenv  MPIEXEC_MPT_DEBUG 0
   setenv MP_DEBUG_NOTIMEOUT yes
 >>>  This has no -n, so the default 1 process is launched.
      CLM launches TOTALPES, which isn't used here
   setenv          LAUNCHCMD mpiexec

   echo "jobname        : $JOBNAME"
   echo "numcpus        : $NUMCPUS"
   echo "numtasks       : $NUMTASKS"
   echo "numnodes       : $NUMNODES"
   echo "tasks_per_node : $NUMTASKS_PERNODE"
   echo " "

else if ($?LSB_HOSTS) then

   # LSF environment variables:
   # env | grep LS | grep -v LS_COLORS | sort

   setenv ORIGINALDIR $LS_SUBCWD
   setenv     JOBNAME $LSB_OUTPUTFILE:ar
   setenv       JOBID $LSB_JOBID
   setenv     MYQUEUE $LSB_QUEUE
   setenv   NODENAMES ${LSB_HOSTS}
   setenv MP_DEBUG_NOTIMEOUT yes
   setenv LAUNCHCMD mpirun.lsf
#  untested method for determining NUMTASKS_PERNODE with LSF
#  setenv NUMTASKS_PERNODE \
#     `echo $LSB_SUB_RES_REQ | sed -ne '/ptile/s#.*\[ptile=\([0-9][0-9]*\)]#\1#p'`

endif

>>> Rename blocks
# ==============================================================================
# Block 0: Set command environment
# ==============================================================================
# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

echo "`date` -- BEGIN CESM_ASSIMILATE"

set nonomatch      # suppress "rm" warnings if wildcard does not match anything
# As of CESM2.0, the assimilate.csh is called by CESM - and has
# two arguments: the CASEROOT and the DATA_ASSIMILATION_CYCLE
setenv CASEROOT $1

# CESM uses C indexing on loops; cycle = [0,....,$DATA_ASSIMILATION_CYCLES - 1]
# That is "fixed" in the assimilation.cshs, so don't fix it here.
@ cycle = $2

# Fetch CESM job parameters from a shared file
source ${CASEROOT}/DART_params || exit 1

cd ${CASEROOT}

setenv scomp                     `./xmlquery COMP_ATM      --value`
setenv CASE                      `./xmlquery CASE          --value`
setenv ENSEMBLE_SIZE             `./xmlquery NINST_ATM     --value`
# setenv CAM_DYCORE                `./xmlquery CAM_DYCORE    --value`
setenv EXEROOT                   `./xmlquery EXEROOT       --value`
setenv RUNDIR                    `./xmlquery RUNDIR        --value`
setenv ARCHIVE                   `./xmlquery DOUT_S_ROOT   --value`
setenv TOTALPES                  `./xmlquery TOTALPES      --value`
setenv CONT_RUN                  `./xmlquery CONTINUE_RUN  --value`
setenv CHECK_TIMING              `./xmlquery CHECK_TIMING  --value`
setenv DATA_ASSIMILATION_CYCLES  `./xmlquery DATA_ASSIMILATION_CYCLES --value`

# Switch CESM's timer script off for the rest of the forecasts of this job.
# The timer takes a significant amount of time, which grows by ~15 s
# for each cycle.  This can double the cycle time in a 2 week job.

./xmlchange CHECK_TIMING=FALSE

cd ${RUNDIR}

# A switch to save all the inflation files
setenv save_all_inf TRUE

# This may be needed before the short-term archiver has been run.
if (! -d ${ARCHIVE}/esp/hist) mkdir -p ${ARCHIVE}/esp/hist

# ==============================================================================
# Block 1: Determine time of current model state from file name of member 1
# These are of the form "${CASE}.cam_${ensemble_member}.i.2000-01-06-00000.nc"
# ==============================================================================

# Piping stuff through 'bc' strips off any preceeding zeros.

# CESM3: the date is in the file name; will be easier to extract
set FILE = `head -n 1 rpointer.cpl_0001`
set FILE = $FILE:r
set DATE_EXT = $FILE:e
set DATE     = `echo $FILE:e | sed -e "s#-# #g"`
set YEAR     = `echo $DATE[1] | bc`
set MONTH    = `echo $DATE[2] | bc`
set DAY      = `echo $DATE[3] | bc`
set SECONDS  = `echo $DATE[4] | bc`
set HOUR     = `echo $DATE[4] / 3600 | bc`

echo "valid time of model is $YEAR $MONTH $DAY $SECONDS (seconds)"
echo "valid time of model is $YEAR $MONTH $DAY $HOUR (hours)"

# Move the hidden restart set back into $rundir so that it is processed properly.

${LIST} -d ../Hide*
if ($status == 0) then
   echo 'Moving hidden restarts into the run directory so they can be used or purged.'
   ${MOVE} ../Hide*/* .
   rmdir   ../Hide*
endif

# ==============================================================================
# Block 2: Populate a run-time directory with the input needed to run DART.
# ==============================================================================

echo "`date` -- BEGIN COPY BLOCK"

echo "`date` -- END COPY BLOCK"

# StrCpl; each assimilator can do this.
# # If possible, use the round-robin approach to deal out the tasks.
# # This facilitates using multiple nodes for the simultaneous I/O operations.
# if ($?NUMTASKS_PERNODE) then
# if ($#NUMTASKS_PERNODE > 0) then
#    sed -i -e "s#layout.*#layout = 2#" \
#           -e "s#tasks_per_node.*#tasks_per_node = $NUMTASKS_PERNODE#" \
#        input.nml || exit 21
# endif
# endif

# ==============================================================================
# Block 3: Identify requested output stages, warn about redundant output.
# ==============================================================================

set MYSTRING = `grep stages_to_write input.nml`
set MYSTRING = (`echo $MYSTRING | sed -e "s#[=,'\.]# #g"`)
set STAGE_input     = FALSE
set STAGE_forecast  = FALSE
set STAGE_preassim  = FALSE
set STAGE_postassim = FALSE
set STAGE_analysis  = FALSE
set STAGE_output    = FALSE

# Assemble lists of stages to write out, which are not the 'output' stage.

set stages_except_output = "{"
@ stage = 2
while ($stage <= $#MYSTRING)
   if ($MYSTRING[$stage] == 'input') then
      set STAGE_input = TRUE
      if ($stage > 2) set stages_except_output = "${stages_except_output},"
      set stages_except_output = "${stages_except_output}input"
   endif
   if ($MYSTRING[$stage] == 'forecast') then
      set STAGE_forecast = TRUE
      if ($stage > 2) set stages_except_output = "${stages_except_output},"
      set stages_except_output = "${stages_except_output}forecast"
   endif
   if ($MYSTRING[$stage] == 'preassim') then
      set STAGE_preassim = TRUE
      if ($stage > 2) set stages_except_output = "${stages_except_output},"
      set stages_except_output = "${stages_except_output}preassim"
   endif
   if ($MYSTRING[$stage] == 'postassim') then
      set STAGE_postassim = TRUE
      if ($stage > 2) set stages_except_output = "${stages_except_output},"
      set stages_except_output = "${stages_except_output}postassim"
   endif
   if ($MYSTRING[$stage] == 'analysis') then
      set STAGE_analysis = TRUE
      if ($stage > 2) set stages_except_output = "${stages_except_output},"
      set stages_except_output = "${stages_except_output}analysis"
   endif
   if ($stage == $#MYSTRING) then
      set stages_all = "${stages_except_output}"
      if ($MYSTRING[$stage] == 'output') then
         set STAGE_output = TRUE
         set stages_all = "${stages_all},output"
      endif
   endif
   @ stage++
end

# Add the closing }
set stages_all = "${stages_all}}"
set stages_except_output = "${stages_except_output}}"

# Checking
echo "stages_except_output = $stages_except_output"
echo "stages_all = $stages_all"
if ($STAGE_output != TRUE) then
   echo "ERROR: assimilate.csh requires that input.nml:filter_nml:stages_to_write includes stage 'output'"
   exit 40
endif

# ==============================================================================
# Block 4: Preliminary clean up, which can run in the background.
# ==============================================================================
# CESM2_0's new archiver has a mechanism for removing restart file sets,
# which we don't need, but it runs only after the (multicycle) job finishes.
# We'd like to remove unneeded restarts as the job progresses, allowing more
# cycles to run before needing to stop to archive data.  So clean them out of
# RUNDIR, and st_archive will never have to deal with them.
# ------------------------------------------------------------------------------

# We need to know the names of the current cesm.log files - one log file is created
# by each CESM model advance.

set log_list = `${LIST} -t cesm.log.*`

echo "most recent log is $log_list[1]"
echo "oldest      log is $log_list[$#log_list]"
echo "entire log list is $log_list"
echo " "


# For safety, leave the most recent *2* restart sets in place.
# Prevents catastrophe if the last restart set is partially written before a crash.
# Add 1 more because the restart set used to start this will be counted:
# there will be 3 restarts when there are only 2 cesm.log files,
# which caused all the files to be deleted.

if ($#log_list >= 3) then

   # List of potential restart sets to remove. 
   # Assimilations should use CESM's multi-driver for good performance,
   # so the coupler restart files will have an 'instance' string in them. 

   set re_list = `${LIST} -t *cpl_0001.r.*`

   if ($#re_list < 3) then
      echo "ERROR: Too many cesm.log files ($#log_list) for the $#re_list restart sets."
      echo "       Clean out the cesm.log files from failed cycles."
      exit 50
   endif

   # Find the date of the oldest restart set from filenames like:
   # setup_test.cpl_0001.r.2016-12-11-21600.nc   
   #
   # Grab the root of the filename (removes the .nc 'extension')
   # and then the extension is the bit we need.
   # Want the YYYY-MM-DD-SSSSS part as well as 'DD-SSSSS'

   set FILE = $re_list[3]
   set FILE = $FILE:r
   if ($FILE:e == 'nc') set FILE = $FILE:r
   set rm_date = $FILE:e
   set RM_DATE_PARTS = `echo $rm_date | sed -e "s#-# #g"`
   set day_o_month = $RM_DATE_PARTS[3]
   set sec_o_day   = $RM_DATE_PARTS[4]
   set day_time    = ${day_o_month}-${sec_o_day}

   # Decide whether to purge restart files at this date and time.
   set save_rest_freq = $save_restart_interval

   set purge = 'true'
   # Learn whether save_rest_freq a string or a number.
   # Character strings must be tested outside of the 'if' statement.
   echo $save_rest_freq | grep '[a-z]'
   if ($status == 0) then
      set purge_date = $RM_DATE_PARTS[1]-$RM_DATE_PARTS[2]-$RM_DATE_PARTS[3]
      set weekday = `date --date="$purge_date" +%A`
      if ($weekday == $save_rest_freq) set purge = 'false'
 
   # Numbers can be tested inside the 'if' statement.
   else if (`echo $save_rest_freq | grep '[0-9]'`) then
      if (${day_o_month} % ${save_rest_freq} == 0) set purge = 'false'
 
   endif

   # Identify log files to be removed or moved.
   # [3] means the 3rd oldest restart set is being (re)moved.
   set rm_log = `echo $log_list[3] | sed -e "s/\./ /g;"`
   set rm_slot = $#rm_log
   if ($rm_log[$#rm_log] == 'gz') @ rm_slot--
   echo 'oldest restart set is from job tagged $rm_log['$rm_slot']='$rm_log[$rm_slot]

   # This first half of the statement removes unwanted restarts.
   # The 'else' block preserves the restarts in the archive directory.

   if ( $sec_o_day !~ '00000' || \
       ($sec_o_day =~ '00000' && $purge == 'true') ) then

      # Optionally save inflation restarts, even if it's not a 'save restart' time.
      if ($save_all_inf =~ TRUE) ${MOVE} ${CASE}*inf*${day_time}*  ${ARCHIVE}/esp/hist

      # Remove intermediate member restarts,
      # but not DART means, sd, obs_seq, inflation restarts output.
      # Note that *cpl.h[ar]* are retained, and any h#, #>0.

      echo "Removing unneeded restart file set (DD_SSSSS ${day_time}) from RUNDIR: "
      echo "     ${CASE}"'*.{r,rs,rs1,rh0,h0}.*'"${day_time}"
      ${REMOVE}  ${CASE}*.{r,rs,rs1,rh0,h0}.*${day_time}* &

      # Handle .i. separately to avoid sweeping up .i.${scomp}_{in,out}put_{mean,sd,...} files.
      echo "     ${CASE}"'*.i.[0-9]*'"${day_time}"
      ${REMOVE} ${CASE}*.i.[0-9]*${day_time}*  &

      if ($save_stages_freq =~ NONE || $save_stages_freq =~ RESTART_TIMES) then
         # 'output' will have been renamed by the time the purging happens.
         echo "     ${CASE}"'*'[0-9].${stages_except_output}'*'${day_time}
         ${REMOVE}  ${CASE}.*[0-9].${stages_except_output}*${day_time}* &
      endif
   else

      echo "Preserving (compressed) restart file set (DD_SSSSS ${day_time})"

      # Optionally COPY inflation restarts to the same place as the other inflation restarts.
      if ($save_all_inf =~ TRUE) then
          ${COPY} ${CASE}*inf*${day_time}*  ${ARCHIVE}/esp/hist &
      endif

      # Optionally REMOVE stages' ensemble members (not means and sds).
      if ($save_stages_freq =~ NONE ) then
         echo "Removing unneeded stages' ensemble members (DD_SSSSS ${day_time}) from RUNDIR: "
         echo "     ${CASE}"'*'[0-9].${stages_except_output}'*'${day_time}
         ${REMOVE}  ${CASE}.*[0-9].${stages_except_output}*${day_time}* &
      endif

      wait

      # The list of components determines which restarts are compressed by this call.
      # List the large files first (more efficient and reliable).
      # There is another call farther down to compress the DART files every cycle.
      echo "compress.csh started at `date`"
>>> Construct the list of components'' file name parts in DART_params,
    based on the list of active components (from the compset).
      ${CASEROOT}/compress.csh $CASE ${rm_date} $ENSEMBLE_SIZE "clm2 cpl cam cice" "$stages_all"
      if ($status != 0) then
         echo "compress.csh failed at `date`"
         exit 55
      endif
      echo "compress.csh finished at `date`"

      # Save the restart set to archive/rest/$datename,
      # where it will be safe from removes of $component/rest.
      # There is an implicit assumption that some sort of inflation will be used.

      set save_root = ${ARCHIVE}/rest/${rm_date}
      if (! -d $save_root) then
         mkdir -p $save_root
         (${MOVE} ${CASE}*.{r,rs,rs1,rh0,h0}.*${day_time}*  $save_root || exit 60) &
         (${MOVE} ${CASE}*.i.[0-9]*${day_time}*             $save_root || exit 61) &
         (${COPY} *.output*inf*${day_time}*                 $save_root || exit 62) &
         (${MOVE} *0001*${rm_log[$rm_slot]}*                $save_root || exit 63) &
         (${MOVE} cesm*${rm_log[$rm_slot]}*                 $save_root || exit 64) &
      else
         echo "WARNING: $save_root already exists.  Did st_archive make it?"
#          exit 65
      endif
   endif

   # Remove log files: *YYMMDD-HHMMSS*.  Except not da.log files
   ${REMOVE}  [^d]*${rm_log[$rm_slot]}*  &

   # I'd like to remove the CAM .r. files, since we always use the .i. files to do a hybrid start,
   # but apparently CESM needs them to be there, even though it doesn't read fields from them.
   # ${REMOVE}  ${CASE}.cam*.r.*${day_time}.nc &

endif

>>> Generalize renaming, archiving, etc., to handle domains
# ==============================================================================
# Block 7: Actually run the assimilations.
#
# Assimilate obs of one component into all of the (chosen) assimilators.
set ocount = 1
while ($ocount <= $#assimilators)
   # Short name for the component providing the obs
   set obs_comp = $assimilators[$ocount]
   # Used by the cross-assims, not by the self assim:
   set obs_dir = ../../${obs_comp}/run
   if (! -d $obs_dir) then
      echo "ERROR: $obs_dir assimilator directory doesn't exist.  Exiting"
      exit
   endif

   cd $obs_dir
   # Link files needed for assimilation here to preserve minimize changes 
   # to assimilate.csh to preserve 'no coupling' mode.
   # State and history files:
   #   h0 is needed for CAM's PHIS field.
   #   h2 is needed for CLM's domains.
   #   rh is needed by all components for inflation restart files.
   # comp_file_bits limits the links to this component (DART_params).
   ln -s ../../CESM/run/*$comp_file_bits[$ocount]*.{r,i,rh,h0,h2}.*.${DATE_EXT}.nc .
   # Inflation files: not used (yet) in cross-assims.
   # Not needed if assimilate.csh gets 4 args:
   #   ln -s ../../rpointer* .

   # Assimilate matching obs into a component.  Note matching indices.
   # Obs come from archive of real obs files.  File name will be built and added to path
   # by the assimilate.cshs.
   # In these calls the CASEROOT (first arg) names the assimilator (= obs_comps),
   # that is, the path to the assimilate.csh.
   # not the CASENAME of the calling (CESM) job.
   # Obs come from the observation archives (not obs_dir).
   $exproot/$obs_comp/assimilate.csh \
      $exproot/$obs_comp $cycle $DATE_EXT $obspaths[$ocount]
   cd -

   set acount = 1
   while ($acount <= $#assimilators)
      # Short name for the assimila'tor'
      set tor = $assimilators[$acount]
      # Already did self-assim, so skip these matches
      if ($tor !~ $obs_comp ) then
         # Look for $pair in $cross_assims.
         set pair = "$tor<$obs_comp"
         @ p = 1
         while ($p <= $#cross_assims)
            if ($pair =~ $cross_assims[$p]) then
               cd $obs_comp

               # Do the assimilation.  Note that obs_dir is the location of the obs,
               # in contrast to obspaths used in the self-assims.
               ${exproot}/${tor}/assimilate.csh \
                  ${exproot}/${tor}  $cycle  $DATE_EXT  ${obs_dir}
               
               cd -
               break
            endif
            @ p++
         end
      endif
      @ acount++
   end
   @ ocount++
end

# ==============================================================================
# Block 8: Rename the output using the CESM file-naming convention.
# ==============================================================================

# Most of this will be done by the assimilators.
# All of the files will (continue to) have the "CESM" root,
# which is the CASENAME of this case.
# For now (2026-10, development) they will be archived with the same root,
# but in an $EXP directory.
# Later we may want to replace the "CESM" with $EXP, to differentiate between
# output from separate coupled assimilations.

# Rename the observation file and run-time output
# StrCpl: ? combine obs_seq files into 1?  With(out) increments?

# Inflation files have been renamed into .rh. files by the assimilators.
# In the end-of-job wrap-up (next) the list of them is encoded into the dart.r. file.


if ($cycle == $DATA_ASSIMILATION_CYCLES) then
   echo "`date` -- BEGIN (NON-RESTART) ARCHIVING LOGIC"

   if ($#log_list >= 3) then

      # During the last cycle, hide the previous restart set
      # so that it's not archived, but is available.
      # (Coupled assimilations may need to keep multiple atmospheric
      #  cycles for each ocean cycle.)

      set FILE = $re_list[2]
      set FILE = $FILE:r
      if ($FILE:e == 'nc') set FILE = $FILE:r
      set hide_date = $FILE:e
      set HIDE_DATE_PARTS = `echo $hide_date | sed -e "s#-# #g"`
      set day_o_month = $HIDE_DATE_PARTS[3]
      set sec_o_day   = $HIDE_DATE_PARTS[4]
      set day_time    = ${day_o_month}-${sec_o_day}

      set hidedir = ../Hide_${day_time}
      mkdir $hidedir

      if ($save_all_inf =~ TRUE) then
         # Put the previous and current inflation restarts in the archive directory.
         # (to protect last from st_archive putting them in esp/hist)
         ${MOVE}   ${CASE}*${stages_except_output}*inf*  ${ARCHIVE}/esp/rest

         # Don't need previous inf restarts now, but want them to be archived later.
         # COPY instead of LINK because they'll be moved or used later.
         ${COPY}   ${CASE}*output*inf* ${ARCHIVE}/esp/rest
      else
         # output*inf must be copied back because it needs to be in ${RUNDIR}
         # when st_archive runs to save the results of the following assim
         ${MOVE}   ${CASE}*inf*${day_time}*  $hidedir

         # Don't need previous inf restarts now, but want them to be archived later.
         ${COPY}   $hidedir/${CASE}*output*inf*${day_time}* .
      endif

      # Hide the CESM 'restart' files from the previous cycle (day_time) from the archiver.
      ${MOVE}           ${CASE}*.{r,rs,rs1,rh0,h0,i}.*${day_time}*    $hidedir

      # Move log files: *YYMMDD-HHMMSS.  [2] means the previous restart set is being moved.
      set rm_log = `echo $log_list[2] | sed -e "s/\./ /g;"`
      # -- (decrement by one) skips the gz at the end of the names.
      set rm_slot = $#rm_log
      if ($rm_log[$#rm_log] =~ gz) @ rm_slot--
      ${MOVE}  *$rm_log[$rm_slot]*  $hidedir
   endif

   # Restore CESM's timing logic for the first cycle of the next job.
   cd ${CASEROOT}
   ./xmlchange CHECK_TIMING=${CHECK_TIMING}
   cd ${RUNDIR}

   # Strcpl; should this be done more often than the end of the job?

   # Create a netCDF file which contains the names of DART inflation restart files.
   # This is needed in order to use the CESM st_archive mechanisms for keeping, 
   # in $RUNDIR, history files which are needed for restarts.
   # These special files must be labeled with '.rh.'.
   # St_archive looks in a .r. restart file for the names of these 'restart history' files.
   # DART's inflation files fit the definition of restart history files, so we put .rh. 
   # in their names.  Those file names must be found in a dart.r. file, which is created here.
   # Inflation restart file names for all components will be in this one restart file,
   # since the inflation restart files have the component names in them.

   set inf_list = `ls *output_{prior,post}inf_*.${DATE_EXT}.nc`
   set file_list = 'restart_hist = "./'$inf_list[1]\"
   set i = 2
   while ($i <= $#inf_list)
      set file_list = (${file_list}\, \"./$inf_list[$i]\")
      @ i++
   end
   cat << ___EndOfText >! inf_restart_list.cdl
       netcdf template {  // CDL file which ncgen will use to make a DART restart file
                          // containing just the names of the needed inflation restart files.
       dimensions:
            num_files = $#inf_list;
       variables:
            string  restart_hist(num_files);
            restart_hist:long_name = "DART restart history file names";
       data:
            $file_list;
       }
___EndOfText

   # I chose 'cpl' to replace ${scomp} because all of the components' inflation files
   # go into this file, so the most suitable CESM 'component' is the coupler/driver.
   # Coincidentally (hopefully not confusingly) it also will indicate that these files
   # are from a coupled assimilation.
   ncgen -k netCDF-4 -o ${CASE}.dart.r.cpl.${DATE_EXT}.nc inf_restart_list.cdl
   if ($status == 0) ${REMOVE} inf_restart_list.cdl

   echo "`date` -- END   ARCHIVING LOGIC"

   # DEBUG st_archive by making a shadow copy of this directory.
   module load nco

   if (-d ../run_shadow) ${REMOVE} -f ../run_shadow
   mkdir ../run_shadow

   set comps = ('cam_'  'clm2_'  'mosart_' 'dart')
   set vars  = ('nhfil' 'locfnh' 'locfnh'  'restart_hist')

   # The [1] is in case $LIST has some arguments that are not wanted here.
   foreach f (`$LIST[1]`)
      @ gr_stat = 1
      echo $f | grep '\.r\.'
      if ($status == 0) then
         set c = 1
         while ($c <= $#comps) 
            echo $f | grep $comps[$c]
            if ($status == 0) then
               echo "c = $c for $f"
#                set echo verbose
               set gr_stat = 0
               ncks -O -v $vars[$c] $f ../run_shadow/$f
               break
            endif
            @ c++
         end
      endif
      if ($gr_stat == 1) then
         ${LIST} -l $f >! ../run_shadow/$f
      endif
   end

endif

# ==============================================================================
# Compress the large coupler history files and DART files.
# ==============================================================================

echo "STARTING: compressing coupler history files and DART files at `date`"

# StrCpl: Don't need now, but may/will, so comment it out.
# Ensure the removal of unneeded restart sets and copy of obs_seq.final are finished.
wait

exit 0

#----------------------------------------------------------------------------------
${CASEROOT}/compress.csh $CASE $DATE_EXT $ENSEMBLE_SIZE "hist dart" "$stages_all"
if ($status != 0) then
   echo "ERROR: Compression of coupler history files and DART files failed at `date`"
   # Ensure the removal of unneeded restart sets and copy of obs_seq.final are finished.
   wait
   exit 250
endif

echo "FINISHED: compressing coupler history files and DART files at `date`"
echo "`date` -- END CESM_ASSIMILATE"

# Ensure the removal of unneeded restart sets and copy of obs_seq.final are finished.
wait

exit 0

