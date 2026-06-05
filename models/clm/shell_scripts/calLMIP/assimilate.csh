#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# This script performs an assimilation by directly reading and writing to
# the CLM restart file.
#
# NOTE: 'dart_to_clm' does not currently support updating the 
# prognostic snow variables based on posterior SWE values.
# Consequently, snow DA is not currently supported.
# Implementing snow DA is high on our list of priorities. 

#=========================================================================
# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.
#=========================================================================

echo "`date` -- BEGIN CLM_ASSIMILATE"

# As of CESM2.0, the assimilate.csh is called by CESM - and has
# two arguments: the CASEROOT and the DATA_ASSIMILATION_CYCLE

setenv CASEROOT $1
setenv ASSIMILATION_CYCLE $2

source ${CASEROOT}/DART_params_tower_assim.csh || exit 1

# Python uses C indexing on loops; cycle = [0,....,$DATA_ASSIMILATION_CYCLES - 1]
# "Fix" that here, so the rest of the script isn't confusing.
@ cycle = $ASSIMILATION_CYCLE + 1

# xmlquery must be executed in $CASEROOT.
cd ${CASEROOT}
setenv CASE           `./xmlquery CASE        --value`
setenv ENSEMBLE_SIZE  `./xmlquery NINST_LND   --value`
setenv EXEROOT        `./xmlquery EXEROOT     --value`
setenv RUNDIR         `./xmlquery RUNDIR      --value`
setenv ARCHIVE        `./xmlquery DOUT_S_ROOT --value`
setenv TOTALPES       `./xmlquery TOTALPES    --value`
setenv STOP_N         `./xmlquery STOP_N      --value`
setenv DATA_ASSIMILATION_CYCLES `./xmlquery DATA_ASSIMILATION_CYCLES --value`
setenv TASKS_PER_NODE `./xmlquery MAX_TASKS_PER_NODE --value`

# Most of this syntax can be determined from CASEROOT  ./preview_run
#setenv MPI_RUN_COMMAND "mpiexec -n $TOTALPES"
# Cannot have more processors than model size which is = 78
# Hard code this for single site assimilation of LAI and SOMC
setenv MPI_RUN_COMMAND "mpiexec -n 78"

cd ${RUNDIR}

#=========================================================================
# Block 1: Determine time of model state ... from file name of first member
# of the form "./${CASE}.clm2_${ensemble_member}.r.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#=========================================================================
set RPFILE  = ( `ls -1tr rpointer.lnd_0001.*` )

# Archiving support
#set RPLFILE = ( `ls -1tr rpointer.lnd_*` )
#set RPAFILE = ( `ls -1tr rpointer.atm_*` )
#set RPCFILE = ( `ls -1tr rpointer.cpl_*` )
#set RLFILE   = ( `ls -1tr ${CASE}.clm2*.r.*nc` )
#set RAFILE   = ( `ls -1tr ${CASE}.datm*.r.*nc` )
#set RCFILE  = ( `ls -1tr ${CASE}.cpl*.r.*nc` )
#set RH0FILE = ( `ls -1tr ${CASE}.clm2*.rh0*nc` )
#set RH1FILE = ( `ls -1tr ${CASE}.clm2*.rh1*nc` )
#set RH2FILE = ( `ls -1tr ${CASE}.clm2*.rh2*nc` )
#set RH3FILE = ( `ls -1tr ${CASE}.clm2*.rh3*nc` )
#set H0FILE = ( `ls -1tr ${CASE}.clm2*.h0*nc` )
#set H1FILE = ( `ls -1tr ${CASE}.clm2*.h1*nc` )
#set LATMFILE = ( `ls -lrt atm_*.log.*`)
#set LDRVFILE = ( `ls -lrt drv_*.log.*`)
#set LMEDFILE = ( `ls -lrt med_*.log.*`)
#set LLNDFILE = ( `ls -lrt lnd_*.log.*`)
#set LCESMFILE = ( `ls -lrt cesm.log.*`)  


if ( $#RPFILE < 1 ) then
   echo "ERROR: no rpointer.lnd_0001.<date> file found"
   exit 1
endif

set FILE = `head -n 1 $RPFILE[$#RPFILE]`
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
echo " "
#echo "Outputting additional archiving variables for troubleshooting"
#echo "The accumulated number of assimilation time steps is rpointer.lnd_0001.*=  $#RPFILE"
#echo "The rpointer.lnd_0001 file list is: $RPFILE"
#echo " "
#echo "The accumulated number of rpointer.lnd files are=  $#RPLFILE"
#echo "The rpointer.lnd file list is: $RPLFILE"
#echo " "
#echo "The accumulated number of rpointer.atm files are=  $#RPAFILE"
#echo "The rpointer.atm file list is: $RPAFILE"
#echo " "
#echo "The accumulated number of rpointer.cpl files are=  $#RPCFILE"
#echo "The rpointer.cpl file list is: $RPCFILE"
#echo " "
#echo "The accumulated number of restart files are=  $#RFILE"
#echo "The restart file list is: $RFILE"
#echo " "
#=========================================================================
# Block 2: Get observation sequence file ... or die right away.
#=========================================================================

# The observation file names have a time that matches the stopping time of CLM.
#
# The CLM observations are stored in two sets of directories.
# If you are stopping every 24 hours or more, the obs are in directories like YYYYMM.
# In all other situations the observations come from directories like YYYYMM_6H.
# The only ugly part here is if the first advance and subsequent advances are
# not the same length. The observations _may_ come from different directories.
#
# The contents of the file must match the history file contents if one is using
# the obs_def_tower_mod or could be the 'traditional' +/- 12Z ... or both.
# Since the history file contains the previous days' history ... so must the obs file.

set OBSDIR = `printf %04d%02d    ${LND_YEAR} ${LND_MONTH}`

set OBS_FILE = ${baseobsdir}/${OBSDIR}/obs_seq.${LND_DATE_EXT}

${REMOVE} obs_seq.out

if (  -e   ${OBS_FILE} ) then
   ${LINK} ${OBS_FILE} obs_seq.out || exit 2
else
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... no observation file $OBS_FILE"
   exit 2
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
   exit 3
endif

echo "`date` -- END COPY BLOCK"

#=========================================================================
# Consistency check: verify that estimate_params in DART_params_tower_assim.csh
# (shell layer) matches estimate_params in input.nml (Fortran layer).
# A mismatch means the shell scripts and Fortran executables will disagree
# on whether to perform parameter DA, which will cause errors or silent
# incorrect behavior.  Warn loudly so the user can fix it before filter runs.
#=========================================================================

set NML_PARAM_TRUE = `grep -i 'estimate_params' input.nml | grep -i 'true'  | wc -l`
set NML_PARAM_FALSE = `grep -i 'estimate_params' input.nml | grep -i 'false' | wc -l`

if ( $estimate_params == TRUE && $NML_PARAM_TRUE == 0 ) then
   echo " "
   echo "WARNING ================================================================"
   echo "WARNING: estimate_params = TRUE in DART_params_tower_assim.csh"
   echo "WARNING: but estimate_params = .false. (or missing) in input.nml"
   echo "WARNING: Shell scripting will attempt param DA but Fortran executables"
   echo "WARNING: will NOT -- this will cause errors. Fix before proceeding."
   echo "WARNING: Set estimate_params = .true. in input.nml:model_nml to resolve."
   echo "WARNING ================================================================"
   echo " "
endif

if ( $estimate_params == FALSE && $NML_PARAM_TRUE > 0 ) then
   echo " "
   echo "WARNING ================================================================"
   echo "WARNING: estimate_params = FALSE in DART_params_tower_assim.csh"
   echo "WARNING: but estimate_params = .true. in input.nml"
   echo "WARNING: Fortran executables will attempt param DA but shell scripting"
   echo "WARNING: will NOT -- param files will be missing and filter will fail."
   echo "WARNING: Set estimate_params = FALSE in input.nml:model_nml to resolve."
   echo "WARNING ================================================================"
   echo " "
endif

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
# Block 4: DART INFLATION
# IF we are doing inflation, we must take the output inflation files from
# the previous cycle and rename them for input to the current cycle.
# The inflation values change through time and should be archived.
#
# If we need to run fill_inflation_restart,
# we need the links to the input files. So this has to come pretty early.
#
# Every variable in the DART vector needs an inflation value if we
# run with any of the temporally- or spatially-adaptive inflation schemes.
# This means that the variables marked 'NO_COPY_BACK' must still have
# inflation values. This is achieved by running fill_inflation_restart
# and copying those input inflation files into the output files, which
# filter will update. By continually copying the input inflation files
# to the output inflation files before filter runs, every variable in the
# DART vector will have an inflation value.
#==========================================================================
# Make sure this is synced with the setup script
# history domain: h0a -->echo "hist_fincl1 = 'LEAFC','TLAI','TOTSOCM_1m'" >> ${fname}
# vector history domain: h1a -->echo "hist_fincl2 = 'TLAI'" >> ${fname}


set     LND_RESTART_FILENAME = ${CASE}.clm2_0001.r.${LND_DATE_EXT}.nc
set     LND_HISTORY_FILENAME = ${CASE}.clm2_0001.h0a.${LND_DATE_EXT}.nc
set LND_VEC_HISTORY_FILENAME = ${CASE}.clm2_0001.h1a.${LND_DATE_EXT}.nc

# remove any potentially pre-existing links
unlink clm_restart.nc
unlink clm_history.nc
unlink clm_vector_history.nc

${LINK} ${LND_RESTART_FILENAME} clm_restart.nc || exit 4
${LINK} ${LND_HISTORY_FILENAME} clm_history.nc || exit 4
if (  -s   ${LND_VEC_HISTORY_FILENAME} ) then
   ${LINK} ${LND_VEC_HISTORY_FILENAME} clm_vector_history.nc || exit 4
endif

# fill_inflation_restart creates files for all the domains in play,
# with names like input_priorinf_[mean,sd]_d0?.nc These should be renamed
# to be similar to what is created during the cycling. fill_inflation_restart
# only takes a second and only runs once.

if ( -e clm_inflation_cookie ) then

   # If parameter estimation is active, expand member 0001's param file now so
   # fill_inflation_restart can determine domain 4 structure and create d04 inflation files.
   # clm_to_dart requires clm.nc (restart) to exist even when only param expansion is needed.
   if ( $estimate_params == TRUE ) then
      set PARAM_0001 = `printf ${param_file_base}%04d.nc 1`
      if ( ! -e $PARAM_0001 ) then
         echo "ERROR: $PARAM_0001 not found in RUNDIR. Check param file staging in setup script."
         exit 4
      endif
      ${COPY} ${LND_RESTART_FILENAME} clm.nc
      ${LINK} $PARAM_0001 clm_params.nc
      ${EXEROOT}/clm_to_dart >& /dev/null
      unlink clm_params.nc
      ${REMOVE} clm.nc
      # clm_params_expanded.nc now exists for fill_inflation_restart domain 4
   endif

   ${EXEROOT}/fill_inflation_restart || exit  4

   foreach FILE ( input_priorinf_*.nc )
      set NEWBASE = `echo $FILE:r | sed -e "s#input#output#"`
      ${MOVE} ${FILE} clm_${NEWBASE}.1601-01-01-00000.nc
   end

   # Make sure this only happens once. Eat the cookie.
   ${REMOVE} clm_inflation_cookie

   # To help keep track of the most recent inflation file,
   # create a 'pointer file' to hold the name of the most recent.

   @ domaincount = 0
   foreach FILE ( clm_output_priorinf_mean*.nc )

      @ domaincount ++

      set POINTERFILE = `printf priorinf_pointer_d%02d.txt $domaincount`

      set SDFILE = `echo $FILE | sed -e "s#mean#sd#"`

      echo $FILE   >! $POINTERFILE
      echo $SDFILE >> $POINTERFILE

   end

   # Not supporting posterior inflation at this time.
   ${REMOVE} input_postinf*nc

endif

# We have to potentially deal with files like:
# clm_output_priorinf_mean_d01.${LND_DATE_EXT}.nc
# clm_output_priorinf_mean_d02.${LND_DATE_EXT}.nc
# clm_output_priorinf_mean_d03.${LND_DATE_EXT}.nc
# clm_output_priorinf_sd_d01.${LND_DATE_EXT}.nc
# clm_output_priorinf_sd_d02.${LND_DATE_EXT}.nc
# clm_output_priorinf_sd_d03.${LND_DATE_EXT}.nc


# Check to see if inflation is being used.

set  MYSTRING = `grep inf_flavor input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  PRIOR_INF = $MYSTRING[2]
set  POSTE_INF = $MYSTRING[3]

if ( $PRIOR_INF != 0 ) then

   # CLM always has at least two domains, but may sometimes have three.
   # Link to the new expected name, if the file does not exist, filter will
   # die and issue a very explicit death message.

   ${REMOVE} input_priorinf_mean*.nc input_priorinf_sd*.nc

   @ domaincount = 1

   foreach POINTERFILE ( priorinf_pointer*.txt )

      set DOMAIN = `printf _d%02d $domaincount`
      set INPUT  =  input_priorinf_mean_${DOMAIN}
      set OUTPUT = output_priorinf_mean_${DOMAIN}

      set latest_mean = `head -n 1 $POINTERFILE`
      set latest_sd   = `tail -n 1 $POINTERFILE`

      # Create the expected output inflation file.
      # The NO_COPY_BACK variables that are part of the DART vector
      # need to have inflation values. 
      ${COPY} ${latest_mean} output_priorinf_mean${DOMAIN}.nc
      ${COPY} ${latest_sd}   output_priorinf_sd${DOMAIN}.nc

      ${LINK} ${latest_mean} input_priorinf_mean${DOMAIN}.nc
      ${LINK} ${latest_sd}   input_priorinf_sd${DOMAIN}.nc

      @ domaincount ++

   end

endif

if ( $POSTE_INF != 0 ) then
   echo "ERROR: assimilate.csh not configured to cycle with posterior inflation."
   exit 4
endif

#=========================================================================
# Block 5: REQUIRED DART namelist settings
#
# "restart_files.txt" is mandatory.
# "history_files.txt" and "history_vector_files.txt" are only needed if
# variables from these files are specified as part of the desired DART state.
# It is an error to specify them if they are not required.
#
# model_nml "clm_restart_filename" and "clm_history_filename" are mandatory
# and are used to determine the domain metadata and *shape* of the variables.
# "clm_vector_history_filename" is used to determine the shape of the
# variables required to be read from the vector history file. If there are no
# vector-based history variables, 'clm_vector_history_filename' is not used.
#
# &filter_nml
#     async                   = 0,
#     obs_sequence_in_name    = 'obs_seq.out'
#     obs_sequence_out_name   = 'obs_seq.final'
#     init_time_days          = -1,
#     init_time_seconds       = -1,
#     first_obs_days          = -1,
#     first_obs_seconds       = -1,
#     last_obs_days           = -1,
#     last_obs_seconds        = -1,
#     input_state_file_list   = "restart_files.txt",
#                               "history_files.txt",
#                               "vector_files.txt"
#     output_state_file_list  = "restart_files.txt",
#                               "history_files.txt",
#                               "vector_files.txt"
# &model_nml
#     clm_restart_filename        = 'clm_restart.nc'
#     clm_history_filename        = 'clm_history.nc'
#     clm_vector_history_filename = 'clm_vector_history.nc'
# &ensemble_manager_nml
#     single_restart_file_in  = .false.
#     single_restart_file_out = .false.
#=========================================================================
# clm always needs a clm_restart.nc, clm_history.nc for geometry information, etc.
# it may or may not need a vector-format history file - depends on user input

${REMOVE} restart_files.txt history_files.txt vector_files.txt param_files.txt

# enscount tracks the ensemble member number for param file naming.
@ enscount = 1

foreach FILE ( ${CASE}.clm2_*.r.${LND_DATE_EXT}.nc )

   # Create unique output filename for the copy that has the indeterminate
   # values replaced by the _FillValue. The copies are the files that will
   # be used to construct the DART vector.
   set OUTPUT = `echo $FILE | sed -e "s/${CASE}.//"`

   ${COPY} $FILE clm.nc

   # If parameter estimation is active, link the working param file for this member
   # so clm_to_dart performs spatial expansion alongside the snow layer fix.
   # The working copy in RUNDIR is used -- the original paramdir files are never touched.
   if ( $estimate_params == TRUE ) then
      set WORK_PARAM = `printf ${param_file_base}%04d.nc $enscount`
      if ( ! -e $WORK_PARAM ) then
         echo "ERROR: $WORK_PARAM not found in RUNDIR. Check param file staging."
         exit 5
      endif
      ${LINK} $WORK_PARAM clm_params.nc
   endif

   ${EXEROOT}/clm_to_dart >& /dev/null

   ${MOVE} clm.nc $OUTPUT

   # Move the per-member expanded param file to a uniquely named file and
   # append its name to param_files.txt for filter to read.
   if ( $estimate_params == TRUE ) then
      set EXP_PARAM = `printf clm_params_expanded_%04d.nc $enscount`
      ${MOVE} clm_params_expanded.nc $EXP_PARAM
      unlink clm_params.nc
      echo $EXP_PARAM >> param_files.txt
   endif

   @ enscount ++

end

ls -1         clm2_*.r.${LND_DATE_EXT}.nc  >! restart_files.txt
ls -1 ${CASE}.clm2_*.h0a.${LND_DATE_EXT}.nc >! history_files.txt
ls -1 ${CASE}.clm2_*.h1a.${LND_DATE_EXT}.nc >! vector_files.txt

#=========================================================================
# Block 6: Actually run the assimilation.
#=========================================================================

echo "`date` -- BEGIN FILTER"
${MPI_RUN_COMMAND} ${EXEROOT}/filter || exit 6
echo "`date` -- END FILTER"

#=========================================================================
# Block 7: Put the DART posterior into the CLM restart file. The CLM
# restart file is also the prior for the next forecast.
#=========================================================================
# Unlink any potentially pre-existing links
unlink clm_restart.nc
unlink dart_posterior.nc

# Identify if SWE re-partitioning is necessary
set  REPARTITION = `grep repartition_swe input.nml`
set  REPARTITION = `echo $REPARTITION | sed -e "s/repartition_swe//g"`
set  REPARTITION = `echo $REPARTITION | sed -e "s/=//g"`


# Reset enscount for Block 7 loops (was incremented in Block 5).
@ enscount = 1

if ($REPARTITION != 0) then
unlink clm_vector_history

   foreach RESTART ( ${CASE}.clm2_*.r.${LND_DATE_EXT}.nc )

     set POSTERIOR_RESTART = `echo $RESTART | sed -e "s/${CASE}.//"`
     set POSTERIOR_VECTOR  = `printf analysis_member_00%02d_d03.nc $enscount`
     set CLM_VECTOR        = `printf ${CASE}.clm2_00%02d.h2.${LND_DATE_EXT}.nc $enscount`

     # Confirm that H2OSNO prior/posterior files exist
     if (! -e $POSTERIOR_VECTOR || ! -e $CLM_VECTOR) then
        echo "ERROR: assimilate.csh could not find $POSTERIOR_VECTOR or $CLM_VECTOR"
        echo "When SWE re-partitioning is enabled H2OSNO must be"
        echo "within vector history file (h2). Also the analysis"
        echo "stage must be output in 'stages_to_write' within filter_nml"
        exit 7
     endif

     ${LINK} $POSTERIOR_RESTART dart_posterior.nc
     ${LINK} $POSTERIOR_VECTOR  dart_posterior_vector.nc
     ${LINK} $RESTART           clm_restart.nc
     ${LINK} $CLM_VECTOR        clm_vector_history.nc

     # Link param files for dart_to_clm parameter averaging (domain 4).
     if ( $estimate_params == TRUE ) then
        set EXP_PARAM  = `printf clm_params_expanded_%04d.nc $enscount`
        set WORK_PARAM = `printf ${param_file_base}%04d.nc $enscount`
        ${LINK} $EXP_PARAM  clm_params_expanded.nc
        ${LINK} $WORK_PARAM clm_params.nc
     endif

     ${EXEROOT}/dart_to_clm >& /dev/null

     if ($status != 0) then
        echo "ERROR: dart_to_clm failed for $RESTART"
        exit 7
     endif

     foreach LIST ( clm_restart.nc clm_vector_history.nc \
                    dart_posterior.nc dart_posterior_vector.nc )
        unlink $LIST
     end

     if ( $estimate_params == TRUE ) then
        unlink clm_params_expanded.nc
        unlink clm_params.nc
     endif

     @ enscount ++
  end


else

foreach RESTART ( ${CASE}.clm2_*.r.${LND_DATE_EXT}.nc )

   set POSTERIOR = `echo $RESTART | sed -e "s/${CASE}.//"`

   ${LINK} $POSTERIOR dart_posterior.nc
   ${LINK} $RESTART   clm_restart.nc

   # Link param files for dart_to_clm parameter averaging (domain 4).
   if ( $estimate_params == TRUE ) then
      set EXP_PARAM  = `printf clm_params_expanded_%04d.nc $enscount`
      set WORK_PARAM = `printf ${param_file_base}%04d.nc $enscount`
      ${LINK} $EXP_PARAM  clm_params_expanded.nc
      ${LINK} $WORK_PARAM clm_params.nc
   endif

   ${EXEROOT}/dart_to_clm >& /dev/null

   if ($status != 0) then
      echo "ERROR: dart_to_clm failed for $RESTART"
      exit 8
   endif

   unlink dart_posterior.nc
   unlink clm_restart.nc

   if ( $estimate_params == TRUE ) then
      unlink clm_params_expanded.nc
      unlink clm_params.nc
   endif

   @ enscount ++
end

endif

# Remove the copies that we no longer need. The posterior values are
# in the DART diagnostic files for the appropriate 'stage'.
\rm -f clm2_*.r.${LND_DATE_EXT}.nc
\rm -f clm2_*.h0a.${LND_DATE_EXT}.nc
\rm -f clm2_*.h1a.${LND_DATE_EXT}.nc

# Remove the expanded param files -- they are intermediate 2D representations
# only needed between Block 5 (clm_to_dart) and Block 7 (dart_to_clm).
# The updated compact param files (calLMIP_precalibratedXXXX.nc) remain in RUNDIR
# as the prior for the next forecast cycle.
if ( $estimate_params == TRUE ) then
   \rm -f clm_params_expanded_*.nc
endif

#=========================================================================
# Block 8: Archive the results and prepare pointer inflation files for
# next cycle. Tag the output with the valid time of the model state.
#=========================================================================

# TODO could move each ensemble-member file to the respective member dir.

foreach FILE ( input*mean*nc      input*sd*nc     input_member*nc \
            forecast*mean*nc   forecast*sd*nc  forecast_member*nc \
            preassim*mean*nc   preassim*sd*nc  preassim_member*nc \
           postassim*mean*nc  postassim*sd*nc postassim_member*nc \
            analysis*mean*nc   analysis*sd*nc  analysis_member*nc \
              output*mean*nc     output*sd*nc )

   if (  -e $FILE ) then
      set FEXT  = $FILE:e
      set FBASE = $FILE:r
      ${MOVE} $FILE clm_${FBASE}.${LND_DATE_EXT}.${FEXT}
   endif
end

# Tag the DART observation file with the valid time of the model state.

${MOVE} obs_seq.final    clm_obs_seq.${LND_DATE_EXT}.final
${MOVE} dart_log.out     clm_dart_log.${LND_DATE_EXT}.out

echo "Updating inflation pointer files."

@ domaincount = 0
foreach FILE ( clm_output_priorinf_mean*.${LND_DATE_EXT}.nc )
   @ domaincount ++
   set POINTERFILE = `printf priorinf_pointer_d%02d.txt $domaincount`
   set SDFILE = `echo $FILE | sed -e "s#mean#sd#"`
   echo $FILE   >! $POINTERFILE
   echo $SDFILE >> $POINTERFILE
end

#-------------------------------------------------------------------------
# Cleanup and archiving
#-------------------------------------------------------------------------
# Periodically archive CESM log, restart and history files
# For rpointer and restart files leave behind current time
#
#if ( $#RPFILE > 1 ) then

#   echo "WARNING, WARNING  Starting archiving of  CESM files"

#   if ( ! -e ${archdir} ) then
#      mkdir -p ${archdir}/rest
#      mkdir -p ${archdir}/logs
#      mkdir -p ${archdir}/lnd/hist
#   endif

 #  mkdir -p ${archdir}/rest/${LND_DATE_EXT}

  # @ nfiles = $#RPLFILE - $num_instances
  # echo "   "
  # echo "Within archive loop"
  # echo "Loop counter is nfiles = #RPLFILE - num_instances"
  # echo "nfiles = $nfiles"
  # echo "#RPLFILE= $#RPLFILE"
  # echo "num_instances= $num_instances"
  # echo " "
  # if ( $nfiles > 0 ) then
  #    @ i = 1
  #    while ( $i <= $nfiles )
         #${MOVE} $RPLFILE[$i] ${archdir}/rest/${LND_DATE_EXT}
         #${MOVE} $RPCFILE[$i] ${archdir}/rest/${LND_DATE_EXT}
         #${MOVE} $RPAFILE[$i] ${archdir}/rest/${LND_DATE_EXT}
         #${MOVE} $RLFILE[$i]   ${archdir}/rest/${LND_DATE_EXT}
         #${MOVE} $RAFILE[$i]   ${archdir}/rest/${LND_DATE_EXT}
         #${MOVE} $RCFILE[$i]   ${archdir}/rest/${LND_DATE_EXT}
         #${MOVE} $RH0FILE[$i]  ${archdir}/rest/${LND_DATE_EXT}
         #${MOVE} $RH1FILE[$i]  ${archdir}/rest/${LND_DATE_EXT}
         #${MOVE} $RH2FILE[$i]  ${archdir}/rest/${LND_DATE_EXT}
         #${MOVE} $RH3FILE[$i]  ${archdir}/rest/${LND_DATE_EXT}
         #${MOVE} $H0FILE[$i]  ${CASE}.clm2_*.h0*nc ${archdir}/lnd/hist
         #${MOVE} $H1FILE[$i]  ${CASE}.clm2_*.h1*nc ${archdir}/lnd/hist
         #${MOVE} $LATMFILE[$i]  ${archdir}/logs
         #${MOVE} $LDRVFILE[$i]  ${archdir}/logs
         #${MOVE} $LMEDFILE[$i]  ${archdir}/logs
         #${MOVE} $LLNDFILE[$i]  ${archdir}/logs
         #${MOVE} $LCESMFILE[$i] ${archdir}/logs  
   #    @ i++
   #   end
   #endif



   #${MOVE} ${CASE}.clm2_*.h2*nc ${archdir}/lnd/hist
   #${MOVE} ${CASE}.clm2_*.h3*nc ${archdir}/lnd/hist


#endif


echo "`date` -- END CLM_ASSIMILATE"

exit 0

