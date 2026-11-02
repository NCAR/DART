#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#

# ------------------------------------------------------------------------------
# Purpose: assimilate with a CAM ensemble in the coupled assimilation context.
#
# This script is called by DART's CESM's multi_assim.csh with 4 arguments:
# 1) the CASEROOT, 
# 2) the assimilation cycle number in this CESM job
# 3) the assimilation date+time (CESM's date format)
# 4) the location of the obs_seq file to assimilate 
#    (may be a "normal" archived obs_seq.out or an obs_seq.final containing increments
#    from a different component.
# ------------------------------------------------------------------------------
# This script imports parameters from a DART_params file in models/CESM;
# experiment, specific hardware, and other configurations. 
#
# Tag DART's output with names using CESM's convention:
#    ${case}.${scomp}[_$inst].${filetype}[.$dart_file].${date}.nc
#    These should all be named with $scomp = "cam" to distinguish
#    them from the same output from other components in multi-component assims.
#    Dart_files are generally DART's generic file names.

# Not from CESM's DART_params; that specifies what some job params WILL BE,
# not necessarily current values.
# SLURM and LSF were removed for now.

if ($?PBS_NODEFILE) then

   # PBS environment variables:
   # env | grep PBS | sort

#    setenv ORIGINALDIR $PBS_O_WORKDIR
   setenv       JOBID $PBS_JOBID
   setenv     MYQUEUE $PBS_O_QUEUE

   setenv     JOBNAME $PBS_JOBNAME
   setenv     NUMCPUS $NCPUS
   setenv    NUMTASKS `cat  $PBS_NODEFILE | wc -l`
   setenv    NUMNODES `uniq $PBS_NODEFILE | wc -l`
   setenv MPI_RUN_COMMAND   "mpiexec -n $NUMCPUS"

   echo "jobname        : $JOBNAME"
   echo "numcpus        : $NUMCPUS"
   echo "numtasks       : $NUMTASKS"
   echo "numnodes       : $NUMNODES"
   echo "tasks_per_node : $use_tasks_per_node"
   echo " "

endif

echo "`date` -- BEGIN CAM_ASSIMILATE"

# ==============================================================================
# Block 0: Set command environment
# ==============================================================================
# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

# suppress "rm" warnings if wildcard does not match anything
set nonomatch

# ----------------------------------------
# This script is called from the CESM CASEROOT/multi_assim.csh in CESM's CASEROOT.
# Import env variables from the CESM CASE (which is where this script is called).
# (But the CASEROOT here is the assimilator's.)

source ./DART_params || exit 1

# ----------------------------------------
# Import variables from arguments
# The assimilator CASEROOT (has the assimilate.csh and other scripts it needs).
# (Will be either a full pathname
setenv CASEROOT $1
cd ${CASEROOT}

# CESM uses C indexing on loops; cycle = [0,....,$DATA_ASSIMILATION_CYCLES - 1]
# "Fix" that here, so the rest of the script isn't confusing.

@ cycle = $2 + 1

setenv scomp                     `./xmlquery COMP_ATM      --value`
setenv CAM_DYCORE                `./xmlquery CAM_DYCORE    --value`
setenv EXEROOT                   `./xmlquery EXEROOT       --value`
setenv RUNDIR                    `./xmlquery RUNDIR        --value`

# Determine what kind of coupling is involved with this assimilation;
#    self-assimilation
#      2 parameters from CESM call + date+time + obs location from multi_assim.csh
#      Inflation should be applied.
#    cross-assimilation
#      same as self-assimilation but inflation is NOT done (until shown to be needed).
if ($#argv == 4) then 
   setenv obs_path $4
   
   set obs_parts = (`echo $obs_path | sed -e "s#\/# #g"`)

   setenv cpld 'self'
   setenv tor_obs "$scomp"

   if ($obs_parts[1] == '..' && $obs_parts[2] == '..') then
      setenv cpld 'cross'
      set obs_comp = "$obs_parts[3]"

      # Cross_assims need to have both components in their names.
      # Put 'AxO' between the Assimilator and the Obs source
      setenv tor_obs "${scomp}AxO${obs_comp}"

   endif

   setenv CASE                      CESM
   setenv archive                   $archdir
   setenv ensemble_size             $num_instances
#  Not used (anymore)  setenv TOTALPES                  `./xmlquery TOTALPES      --value`
#  Need CONTINUE_RUN from the CESM CASEROOT for this cycle, 
#  not from the CASE's first cycle (DART_params).
   cd ../CESM
   setenv CONT_RUN                  `./xmlquery CONTINUE_RUN  --value`
   setenv DATA_ASSIMILATION_CYCLES  `./xmlquery DATA_ASSIMILATION_CYCLES --value`
   cd -

   set hindcast_path = '../../CESM/run/'

else
   echo "ERROR: cam-fv's assimilate.csh must receive 4 arguments"
   exit 7
endif
   
cd ${RUNDIR}

# ==============================================================================
# Block 1: Determine time of current model state from file name of member 1
# These are of the form "${CASE}.cam_${ensemble_member}.i.2000-01-06-00000.nc"
# ==============================================================================

# Piping stuff through 'bc' strips off any preceeding zeros.
set ATM_DATE_EXT = "$3"
set ATM_DATE     = `echo $ATM_DATE_EXT | sed -e "s#-# #g"`
set ATM_YEAR     = `echo $ATM_DATE[1] | bc`
set ATM_MONTH    = `echo $ATM_DATE[2] | bc`
set ATM_DAY      = `echo $ATM_DATE[3] | bc`
set ATM_SECONDS  = `echo $ATM_DATE[4] | bc`
set ATM_HOUR     = `echo $ATM_DATE[4] / 3600 | bc`

echo "valid time of model is $ATM_YEAR $ATM_MONTH $ATM_DAY $ATM_SECONDS (seconds)"
echo "valid time of model is $ATM_YEAR $ATM_MONTH $ATM_DAY $ATM_HOUR (hours)"

# ==============================================================================
# Block 2: Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of CAM.
#
# Make sure the file name structure matches the obs you will be using.
# PERFECT model obs output appends .perfect to the filenames
#
# Coupled assims will get the date and obs root as arguments #3 and 4
# If the assimilator is the same as the source of the obs,
# then archived obs_seq.out (not final) files will be assimilated.
# The location is listed in CESM_params:obspaths, passed in through arg #4.
# Otherwise use the obs_seq.final (which has sequential increments, if needed)
# in the path fed in through arg #4, which will start with ../../{obs_component}/run.
# The file name will be translated into CESM format,
# but also linked by obs_seq.final.

# ==============================================================================

if ($cpld =~ 'self') then 
   set YYYYMM = `printf %04d%02d ${ATM_YEAR} ${ATM_MONTH}`
   
   if (! -d ${obs_path}/${YYYYMM}_6H_CESM) then
      echo "CESM+DART requires 6 hourly obs_seq files in directories of the form YYYYMM_6H_CESM"
      echo "The directory ${obs_path}/${YYYYMM}_6H_CESM is not found.  Exiting"
      exit 70
   endif
   
   set OBSFNAME = `printf obs_seq.%04d-%02d-%02d-%05d ${ATM_YEAR} ${ATM_MONTH} ${ATM_DAY} ${ATM_SECONDS}`
   
   set OBS_FILE = "${obs_path}/${YYYYMM}_6H_CESM/${OBSFNAME}"

else if ($cpld =~ 'cross') then
   # ?>>  This was not Jeff''s naming convention, but I think this should work.
   # It does not label the file as a .sequential_prior_$obs_component.
   # Does that matter, if everything is consistent with this system?
   
   # obs_seq.final is renamed as
   # "../../CESM/run/${CASE}.dart.e.${tor_obs}_obs_seq_final.${ATM_DATE_EXT}"
   # tor_obs = $scomp for 'self' and = ${tor}AxO${obs} for cross assims.
   # Only the self-assim obs_seq_final are linked back into its rundir,
   # for use by the cross-assims.
   set OBS_FILE = "${obs_path}/obs_seq.final"
   endif
else 
   echo "ERROR ... assimilate.csh: cpld has illegal value $cpld"
   exit 75 
endif
echo "OBS_FILE = $OBS_FILE"

${REMOVE} obs_seq.out
if (  -e ${OBS_FILE} ) then
   ${LINK} ${OBS_FILE} obs_seq.out || exit 80
else
   echo "ERROR ... no observation file ${OBS_FILE}"
   echo "ERROR ... no observation file ${OBS_FILE}"
   exit 81
endif

# ==============================================================================
# Block 3: Populate a run-time directory with the input.nml needed to run DART.
# ==============================================================================

echo "`date` -- BEGIN COPY BLOCK"

# Put a pared down copy (no comments) of input.nml in this RUNDIR.
# The contents change from one coupled assimilation phase to another,
# and maybe from one cycle to the next, so always start from
# the known configuration in the CASEROOT directory.
# That version should have had preprocessor_nml and obs_kind_nml imported from
# CESM/work/input.nml so that cam-fv can handle obs from all components.

if (  -e   ${CASEROOT}/input.nml ) then

   sed -e "/#/d;/^\!/d;/^[ ]*\!/d" \
       -e "1,1i\WARNING: Changes to this file will be ignored. \n Edit \$CASEROOT/input.nml instead.\n\n\n" \
       ${CASEROOT}/input.nml >! input.nml  || exit 10

   # Self-assimilations need to write out the obs increments, while
   # cross-assims need to read them.  Location_mod needs to know whether it's a cross-assim.
   # Cross-assims do not use inflation.
   # `sed -i` means "edit in place", not "input".
   if ($cpld =~ 'self') then 
      sed -e "s#output_sequential_prior_post.*#output_sequential_prior_post = .true.#" \
          -e "s#use_sequential_prior_post.*#use_sequential_prior_post = .false. #" \
          -e "s#strongly_coupled.*#strongly_coupled = .false.#" \
          -i input.nml  || exit 12
   else if ($cpld =~ 'cross') then 
      sed -e "s#output_sequential_prior_post.*#output_sequential_prior_post = .false.#" \
          -e "s#use_sequential_prior_post.*#use_sequential_prior_post = .true. #" \
          -e "s#strongly_coupled.*#strongly_coupled = .true.#" \
          -e "s#inf_flavor.*#inflation_flavor = 0, 0#" \
          -i input.nml  || exit 13
   endif
else
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   exit 15
endif

echo "`date` -- END COPY BLOCK"

# If possible, use the round-robin approach to deal out the tasks.
# This facilitates using multiple nodes for the simultaneous I/O operations.

if ($?use_tasks_per_node) then
   if ($#use_tasks_per_node > 0) then
      sed -e "s#layout.*#layout = 2#" \
          -e "s#tasks_per_node.*#tasks_per_node = $use_tasks_per_node#" \
          -i input.nml || exit 21
   endif
endif

# ==============================================================================
# Block 6: DART INFLATION
# This block is only relevant if 'inflation' is turned on AND
# inflation values change through time:
# filter_nml
#    inf_flavor(:)  = 2  (or 5)
#    inf_initial_from_restart    = .TRUE.
#    inf_sd_initial_from_restart = .TRUE.

# This block stages the files that contain the inflation values.
# The inflation files are essentially duplicates of the DART model state,
# which have names in the CESM style, something like
#    ${case}.dart.rh.${scomp}_output_priorinf_{mean,sd}.YYYY-MM-DD-SSSSS.nc
# The strategy is to use the latest such files in ${RUNDIR}.
# If those don't exist at the start of an assimilation,
# this block creates them with 'fill_inflation_restart'.
# If they don't exist AFTER the first cycle, the script will exit
# because they should have been available from a previous cycle.
# The script does NOT check the model date of the files for consistency
# with the current forecast time, so check that the inflation mean
# files are evolving as expected.

# CESM's st_archive should archive the inflation restart files
# like any other "restart history" (.rh.) files; copying the latest files
# to the archive directory, and moving all of the older ones.
# ==============================================================================

# If we need to run fill_inflation_restart, CAM:static_init_model()
# always needs a caminput.nc and a cam_phis.nc for geometry information, etc.

set MYSTRING = `grep cam_template_filename input.nml`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set CAMINPUT = "$MYSTRING[2]"
${REMOVE} ${CAMINPUT}
${LINK} ${CASE}.cam_0001.i.${ATM_DATE_EXT}.nc ${CAMINPUT} || exit 90

# All of the .h0. files contain the same PHIS field, so we can link to any of them.
set MYSTRING = `grep cam_phis_filename input.nml`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set CAM_PHIS = "$MYSTRING[2]"
# Listings are slow in productions runs; only do it if necessary
if (! -f ${CAM_PHIS}) then
   set hists = `${LIST} ${CASE}.cam_0001.h0.*.nc`
   ${LINK} $hists[1] ${CAM_PHIS} || exit 100
endif

# Now, actually check the inflation settings

set  MYSTRING = `grep inf_flavor input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  PRIOR_INF = "$MYSTRING[2]"
set  POSTE_INF = "$MYSTRING[3]"

set  MYSTRING = `grep inf_initial_from_restart input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`

# If no inflation is requested, the inflation restart source is ignored

if ( $PRIOR_INF == 0 ) then
   set  PRIOR_INFLATION_FROM_RESTART = "ignored"
   set  USING_PRIOR_INFLATION = "false"
else
   set  PRIOR_INFLATION_FROM_RESTART = `echo $MYSTRING[2] | tr '[:upper:]' '[:lower:]'`
   set  USING_PRIOR_INFLATION = "true"
endif

if ( $POSTE_INF == 0 ) then
   set  POSTE_INFLATION_FROM_RESTART = ignored
   set  USING_POSTE_INFLATION = "false"
else
   set  POSTE_INFLATION_FROM_RESTART = `echo $MYSTRING[3] | tr '[:upper:]' '[:lower:]'`
   set  USING_POSTE_INFLATION = "true"
endif

# StrCpl; can't do this here (yet) because STAGE_* are not available.
#         (They could be, by reimporting the extraction of stages from input.nml)
#         And can't do this test in multi_assim.csh because each component
#         may have different inflation.
# if ($USING_PRIOR_INFLATION == false ) then
#    @ stages_requested = 0
#    if ( $STAGE_input    == TRUE ) @ stages_requested++
#    if ( $STAGE_forecast == TRUE ) @ stages_requested++
#    if ( $STAGE_preassim == TRUE ) @ stages_requested++
#    if ( $stages_requested > 1 ) then
#       echo " "
#       echo "WARNING ! ! Redundant output is requested at multiple stages before assimilation."
#       echo "            Stages 'input' and 'forecast' are always redundant."
#       echo "            Prior inflation is OFF, so stage 'preassim' is also redundant. "
#       echo "            We recommend requesting just 'preassim'."
#       echo " "
#    endif
# endif

# if ($USING_POSTE_INFLATION == false ) then
#    @ stages_requested = 0
#    if ( $STAGE_postassim == TRUE ) @ stages_requested++
#    if ( $STAGE_analysis  == TRUE ) @ stages_requested++
#    if ( $STAGE_output    == TRUE ) @ stages_requested++
#    if ( $stages_requested > 1 ) then
#       echo " "
#       echo "WARNING ! ! Redundant output is requested at multiple stages after assimilation."
#       echo "            Stages 'output' and 'analysis' are always redundant."
#       echo "            Posterior inflation is OFF, so stage 'postassim' is also redundant. "
#       echo "            We recommend requesting just 'output'."
#       echo " "
#    endif
# endif

# IFF we want PRIOR inflation:

if ($USING_PRIOR_INFLATION == true && $cpld == self) then
   if ($PRIOR_INFLATION_FROM_RESTART == false) then

      echo "inf_flavor(1) = $PRIOR_INF, using namelist values."

   else
      # Look for the output from the previous assimilation (or fill_inflation_restart)
      # If inflation files exists, use them as input for this assimilation
      (${LIST} -rt1 *.dart.rh.${scomp}_output_priorinf_mean* | tail -n 1 >! latestfile) > & /dev/null
      (${LIST} -rt1 *.dart.rh.${scomp}_output_priorinf_sd*   | tail -n 1 >> latestfile) > & /dev/null
      @ nfiles = `cat latestfile | wc -l`

      if ( $nfiles > 0 ) then

         set latest_mean = `head -n 1 latestfile`
         set latest_sd   = `tail -n 1 latestfile`
         # Need to COPY instead of link because of short-term archiver and disk management.
         ${COPY} $latest_mean input_priorinf_mean.nc
         ${COPY} $latest_sd   input_priorinf_sd.nc

      else if ($CONT_RUN == FALSE) then

         # It's the first assimilation; try to find some inflation restart files
         # or make them using fill_inflation_restart.
         # Fill_inflation_restart needs caminput.nc and cam_phis.nc for static_model_init,
         # so this staging is done in assimilate.csh (after a forecast) instead of stage_cesm_files.

         if (-x ${EXEROOT}/fill_inflation_restart) then

            ${EXEROOT}/fill_inflation_restart

         else
            echo "ERROR: Requested PRIOR inflation restart for the first cycle."
            echo "       There are no existing inflation files available "
            echo "       and ${EXEROOT}/fill_inflation_restart is missing."
            echo "EXITING"
            exit 112
         endif

      else
         echo "ERROR: Requested PRIOR inflation restart, "
         echo '       but files *.dart.rh.${scomp}_output_priorinf_* do not exist in the ${RUNDIR}.'
         echo '       If you are changing from cam_no_assimilate.csh to assimilate.csh,'
         echo '       you might be able to continue by changing CONTINUE_RUN = FALSE for this cycle,'
         echo '       and restaging the initial ensemble.'
         ${LIST} -l *inf*
         echo "EXITING"
         exit 115
      endif
   endif
else
   echo "Prior Inflation not requested for this assimilation."
endif

# POSTERIOR: We look for the 'newest' and use it - IFF we need it.

if ($USING_POSTE_INFLATION == true && $cpld == self) then
   if ($POSTE_INFLATION_FROM_RESTART == false) then

      # we are not using an existing inflation file.
      echo "inf_flavor(2) = $POSTE_INF, using namelist values."

   else
      # Look for the output from the previous assimilation (or fill_inflation_restart).
      # (The only stage after posterior inflation.)
      (${LIST} -rt1 *.dart.rh.${scomp}_output_postinf_mean* | tail -n 1 >! latestfile) > & /dev/null
      (${LIST} -rt1 *.dart.rh.${scomp}_output_postinf_sd*   | tail -n 1 >> latestfile) > & /dev/null
      @ nfiles = `cat latestfile | wc -l`

      # If one exists, use it as input for this assimilation
      if ( $nfiles > 0 ) then

         set latest_mean = `head -n 1 latestfile`
         set latest_sd   = `tail -n 1 latestfile`
         ${LINK} $latest_mean input_postinf_mean.nc || exit 120
         ${LINK} $latest_sd   input_postinf_sd.nc   || exit 121

      else if ($CONT_RUN == FALSE) then
         # It's the first assimilation; try to find some inflation restart files
         # or make them using fill_inflation_restart.
         # Fill_inflation_restart needs caminput.nc and cam_phis.nc for static_model_init,
         # so this staging is done in assimilate.csh (after a forecast) instead of stage_cesm_files.

         if (-x ${EXEROOT}/fill_inflation_restart) then
            ${EXEROOT}/fill_inflation_restart
            ${MOVE} prior_inflation_mean.nc input_postinf_mean.nc || exit 125
            ${MOVE} prior_inflation_sd.nc   input_postinf_sd.nc   || exit 126

         else
            echo "ERROR: Requested POSTERIOR inflation restart for the first cycle."
            echo "       There are no existing inflation files available "
            echo "       and ${EXEROOT}/fill_inflation_restart is missing."
            echo "EXITING"
            exit 127
         endif

      else
         echo "ERROR: Requested POSTERIOR inflation restart, "
         echo '       but files *.dart.rh.${scomp}_output_postinf_* do not exist in the ${RUNDIR}.'
         ${LIST} -l *inf*
         echo "EXITING"
         exit 128
      endif
   endif
else
   echo "Posterior Inflation not requested for this assimilation."
endif

# ==============================================================================
# Block 7: Actually run the assimilation.
#
# DART namelist settings required:
# &filter_nml
#    adv_ens_command         = "no_CESM_advance_script",
#    obs_sequence_in_name    = 'obs_seq.out'
#    obs_sequence_out_name   = 'obs_seq.final'
#    single_file_in          = .false.,
#    single_file_out         = .false.,
#    stages_to_write         = stages you want + ,'output'
#    input_state_file_list   = 'cam_init_files'
#    output_state_file_list  = 'cam_init_files',
#
# WARNING: the default mode of this script assumes that
#          input_state_file_list = output_state_file_list, so that
#          the CAM initial files used as input to filter will be overwritten.
#          The input model states can be preserved by requesting that stage
#          'forecast' be output.
#
# ==============================================================================

# In the default mode of CAM assimilations, filter gets the model state(s)
# from CAM initial files.  This section puts the names of those files into a text file.
# The name of the text file is provided to filter in filter_nml:input_state_file_list.

# NOTE:
# If the files in input_state_file_list are CESM initial files (all vars and
# all meta data), then they will end up with a different structure than
# the non-'output', stage output written by filter ('preassim', 'postassim', etc.).
# This can be prevented (at the cost of more disk space) by copying
# the CESM format initial files into the names filter will use for preassim, etc.:
#    > cp $case.cam_0001.i.$date.nc  preassim_member_0001.nc.
#    > ... for all members
# Filter will replace the state variables in preassim_member* with updated versions,
# but leave the other variables and all metadata unchanged.

# If filter will create an ensemble from a single state,
#    filter_nml: perturb_from_single_instance = .true.
# it's fine (and convenient) to put the whole list of files in input_state_file_list.
# Filter will just use the first as the base to perturb.

set line = `grep input_state_file_list input.nml | sed -e "s#[=,'\.]# #g"`
set input_file_list_name = "$line[2]"

# If the file names in $output_state_file_list = names in $input_state_file_list,
# then the restart file contents will be overwritten with the states updated by DART.

set line = `grep output_state_file_list input.nml | sed -e "s#[=,'\.]# #g"`
set output_file_list_name = "$line[2]"

if ($input_file_list_name != $output_file_list_name) then
   echo "ERROR: assimilate.csh requires that input_file_list = output_file_list"
   echo "       You can probably find the data you want in stage 'forecast'."
   echo "       If you truly require separate copies of CAM's initial files"
   echo "       before and after the assimilation, see revision 12603, and note that"
   echo "       it requires changing the linking to cam_initial_####.nc, below."
   exit 130
endif

${LIST} -1 ${hindcast_path}${CASE}.cam_[0-9][0-9][0-9][0-9].i.${ATM_DATE_EXT}.nc \
        >! $input_file_list_name

echo "`date` -- BEGIN FILTER"
${MPI_RUN_COMMAND} ${EXEROOT}/filter || exit 140
echo "`date` -- END FILTER"

# ==============================================================================
# Block 8: Rename the output using the CESM file-naming convention.
# ==============================================================================

# If output_state_file_list is filled with custom (CESM) filenames,
# then 'output' ensemble members will not appear with filter's default,
# hard-wired names.  But file types output_{mean,sd} will appear and be
# renamed here.
#
# We don't know the exact set of files which will be written,
# so loop over all possibilities: use LIST in the foreach.
# LIST will expand the variables and wildcards, only existing files will be
# in the foreach loop. (If the input.nml has num_output_state_members = 0,
# there will be no output_member_xxxx.nc even though the 'output' stage
# may be requested - for the mean and sd) 
#
# Handle files with instance numbers first.
#    split off the .nc
#    separate the pieces of the remainder
#    grab all but the trailing 'member' and #### parts.
#    and join them back together

echo "`date` -- BEGIN FILE RENAMING"

# The short-term archiver archives files depending on pieces of their names.
# '_####.i.' files are CESM initial files.
# '.dart.i.' files are ensemble statistics (mean, sd) of just the state variables 
#            in the initial files.
# '.e.'      designates a file as something from the 'external system processing ESP', e.g. DART.

# DART 'output_member_****.nc' files are actually linked to cam input files

foreach FILE (`${LIST} ${stages_all}_member_*.nc`)

#    set parts = `echo $FILE | sed -e "s#\.# #g"`
#    set list = `echo $parts[1]  | sed -e "s#_# #g"`
   set list = `echo $FILE:r | sed -e "s#_# #g"`
   # DART file type doesn't need 'member' and inst number is moved to $scomp.
   # Keep the stage
#    @ last = $#list - 2
#    set dart_file = `echo $list[1-$last] | sed -e "s# #_#g"`
   set dart_file = "$list[1]"
   set inst = "$list[$#list]"
   set type = "e"
   echo $FILE | grep "put"
   if ($status == 0) set type = "i"

   ${MOVE} $FILE \
       ${hindcast_path}${CASE}.${tor_obs}_${inst}.${type}.${dart_file}.${ATM_DATE_EXT}.nc || exit 150
end

# Files without instance numbers need to have the scomp part of their names = "dart".
# This is because in st_archive, all files with  scomp = "cam"
# (= compname in env_archive.xml) will be st_archived using a pattern
# which has the instance number added onto it.  {mean,sd} files don't have 
# instance numbers, so they need to be archived by the "dart" section of env_archive.xml.
# But they still need to be different for each component, so include $scomp in the
# ".dart_file" part of the file name.  Somewhat awkward and inconsistent, but effective.

# Means and standard deviation files (except for inflation).
foreach FILE (`${LIST} ${stages_all}_{mean,sd}*.nc`)

   set parts = $FILE:r
   set type = "e"
   echo $FILE | grep "put"
   if ($status == 0) set type = "i"

   ${MOVE} $FILE \
      ${hindcast_path}${CASE}.dart.${type}.${tor_obs}_${parts}.${ATM_DATE_EXT}.nc || exit 160
end

# Rename files so that st_archive will find them.

set obs_seq_final = "${CASE}.dart.e.${tor_obs}_obs_seq_final.${ATM_DATE_EXT}"
${MOVE} obs_seq.final ${hindcast_path}$obs_seq_final || exit 170

# There will be multiple obs_seq.final files created by this assimilator,
# but only the self-assim one should be linked to a local generic file name.
if ($cpld =~ 'self') then 
   ${LINK} ${hindcast_path}$obs_seq_final obs_seq.final
endif

${MOVE} dart_log.out ${hindcast_path}${tor_obs}_dart_log.${ATM_DATE_EXT}.out || exit 171
gzip input.nml
${MOVE} input.nml.gz ${hindcast_path}${tor_obs}input.nml.log.${JOBID}.gz

# Rename the inflation files and designate them as 'rh' files - which get
# reinstated in the run directory by the short-term archiver and are then
# available for the next assimilation cycle.
#
# Accommodate any possible inflation files.
# The .${scomp}_ part is needed by DART to distinguish
# between inflation files from separate components in coupled assims.

# StrCpl; for now inflation is only applied to self-assims and 'no coupling' cases.
#         So the existing scripting should be sufficient.
#         But don't manage the cross-assim phase files.

if ($cpld =~ 'self') then 
   foreach FILE (`${LIST} ${stages_all}_{prior,post}inf_*`)
   
      set parts = `echo $FILE | sed -e "s#\.# #g"`
      set DEST = "${hindcast_path}${CASE}.dart.rh.${scomp}_$parts[1].${ATM_DATE_EXT}.nc "
      ${MOVE} $FILE $DEST || exit 180
      ${LINK} $DEST .
   
   end
endif

# Handle localization_diagnostics_files
set MYSTRING = `grep 'localization_diagnostics_file' input.nml`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set loc_diag = "$MYSTRING[2]"
if (-f $loc_diag) then
   ${MOVE} $loc_diag \
      ${hindcast_path}${tor_obs}_${loc_diag}.dart.e.${ATM_DATE_EXT} || exit 190
endif

# Handle regression diagnostics
set MYSTRING = `grep 'reg_diagnostics_file' input.nml`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set reg_diag = "$MYSTRING[2]"
if (-f $reg_diag) then
   ${MOVE} $reg_diag \
      ${hindcast_path}${tor_obs}_${reg_diag}.dart.e.${ATM_DATE_EXT} || exit 200
endif

# Then this script will need to feed the files in output_restart_list_file
# to the next model advance.   This is done in . (single-comp assims)
# or in ../../CESM/run in coupled assims (hindcast_path handles the situation).
cd $hindcast_path
# This gets the .i. or .r. piece from the CESM format file name.
set line = `grep 0001 $output_file_list_name | sed -e "s#[\.]# #g"`
@ l = 1
while ($l < $#line)
   if ($line[$l] =~ ${scomp}_0001) then
      @ l++
      set file_type = "$line[$l]"
      break
   endif
   @ l++
end

@ member = 1
while ( ${member} <= ${ensemble_size} )

   set inst_string = `printf _%04d $member`
   set ATM_INITIAL_FILENAME = "${CASE}.${scomp}${inst_string}.${file_type}.${ATM_DATE_EXT}.nc"

   # This is the generic name of the initial file in the user_nl_cam_$member file.
   ${REMOVE} ${scomp}_initial${inst_string}.nc
   ${LINK} $ATM_INITIAL_FILENAME ${scomp}_initial${inst_string}.nc || exit 210

   @ member++

end
cd -
echo "Should be in cam-fv/run"
pwd

echo "`date` -- END   FILE RENAMING"

exit 0
