#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# This script performs an assimilation by directly reading and writing to
# the CLM restart file.
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
#
# Tag DART's state output with names using CESM's convention:
#    ${case}.${scomp}[_$inst].${filetype}[.$dart_file].${date}.nc
#    These should all be named with $scomp = "clm" to distinguish
#    them from the same output from other components in multi-component assims.

#=========================================================================

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
   setenv MPI_RUN_COMMAND "mpiexec -n $NUMCPUS"

   echo "jobname        : $JOBNAME"
   echo "numcpus        : $NUMCPUS"
   echo "numtasks       : $NUMTASKS"
   echo "numnodes       : $NUMNODES"
   echo "tasks_per_node : $use_tasks_per_node"
   echo " "

endif

# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.
#=========================================================================

echo "`date` -- BEGIN CLM_ASSIMILATE"
# suppress "rm" warnings if wildcard does not match anything
set nonomatch

# Import env variables from the CESM CASE (which is where this script is called).
# This is called from the CESM CASEROOT/multi_assim.csh in CESM's CASEROOT.
source ./DART_params || exit 1

# As of CESM2.0, the assimilate.csh is called by CESM - and has
# two arguments: the CASEROOT and the DATA_ASSIMILATION_CYCLE
# StrCpl; when this is called by a strongly-coupled multi_assim.csh
#         it gets 4 arguments; the orginal 2, the date string,
#         and the path to the observations, which may be archived obs_seq.outs
#         or obs_seq.finals (with increments) from another component's assimilation.

setenv CASEROOT $1


# Python uses C indexing on loops; cycle = [0,....,$DATA_ASSIMILATION_CYCLES - 1]
# "Fix" that here, so the rest of the script isn't confusing.
@ cycle = $2 + 1

# xmlquery must be executed in $CASEROOT.
cd ${CASEROOT}

setenv scomp      `./xmlquery COMP_LND    --value`
setenv EXEROOT    `./xmlquery EXEROOT     --value`
setenv RUNDIR     `./xmlquery RUNDIR      --value`

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
      setenv obs_comp "$obs_parts[3]"

      # Cross_assims need to have both components in their names.
      # Put 'AxO' between the Assimilator and the Obs source
      setenv tor_obs "${scomp}AxO${obs_comp}"

   endif


   setenv CASE                      CESM
   setenv ARCHIVE        $archdir
   setenv ENSEMBLE_SIZE  $num_instances
   setenv STOP_N         `./xmlquery STOP_N      --value`
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

#=========================================================================
# Block 1: Determine time of model state ... from file name of first member
# of the form "./${CASE}.clm2_${ensemble_member}.r.2000-01-06-00000.nc"

# Piping stuff through 'bc' strips off any preceeding zeros.
set LND_DATE_EXT = "$3"
set LND_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
set LND_YEAR     = `echo $LND_DATE[1] | bc`
set LND_MONTH    = `echo $LND_DATE[2] | bc`
set LND_DAY      = `echo $LND_DATE[3] | bc`
set LND_SECONDS  = `echo $LND_DATE[4] | bc`
set LND_HOUR     = `echo $LND_DATE[4] / 3600 | bc`

echo "valid time of model is $LND_YEAR $LND_MONTH $LND_DAY $LND_SECONDS (seconds)"
echo "valid time of model is $LND_YEAR $LND_MONTH $LND_DAY $LND_HOUR (hours)"

#=========================================================================
# Block 2: Get observation sequence file ... or die right away.
#=========================================================================

# The observation file names have a time that matches the stopping time of CLM.
#
# Coupled assims will get the date and obs root as arguments #3 and 4
# If the assimilator is the same as the source of the obs,
# then archived obs_seq.out (not final) files will be assimilated.
# The location is listed in CESM_params:obspaths, passed in through arg #4.
# Otherwise use the obs_seq.final (which has sequential increments, if needed)
# in the path fed in through arg #4, which will start with ../../{obs_component}/run.
# The file name will be translated into CESM format,
# but also linked by obs_seq.final.

# The CLM observations are stored in two sets of directories.
# If you are stopping every 24 hours or more, the obs are in directories like YYYYMM.
# In all other situations the observations come from directories like YYYYMM_6H.
# The only ugly part here is if the first advance and subsequent advances are
# not the same length. The observations _may_ come from different directories.
#
# The contents of the file must match the history file contents if one is using
# the obs_def_tower_mod or could be the 'traditional' +/- 12Z ... or both.
# Since the history file contains the previous days' history ... so must the obs file.

if ($cpld =~ 'self') then 
   if ($STOP_N >= 24) then
      set OBSDIR = `printf %04d%02d    ${LND_YEAR} ${LND_MONTH}`
   else
      set OBSDIR = `printf %04d%02d_6H ${LND_YEAR} ${LND_MONTH}`
   endif

   set OBS_FILE = ${obs_path}/${OBSDIR}/obs_seq.${LND_DATE_EXT}
else if ($cpld =~ 'cross') then
   set OBS_FILE = "${obs_path}/obs_seq.final"
   endif
endif
echo "OBS_FILE = $OBS_FILE"

${REMOVE} obs_seq.out

if (  -e   ${OBS_FILE} ) then
   ${LINK} ${OBS_FILE} obs_seq.out || exit 2
else
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... no observation file $OBS_FILE"
   exit 2
endif

#=========================================================================
# Block 3: Populate a run-time directory with the input.nml needed to run DART.
#=========================================================================

# REQUIRED DART namelist settings
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
# &fill_inflation_restart_nml
#     input_state_files = 'clm_restart.nc','clm_history.nc','clm_vector_history.nc'
# &clm_to_dart_nml
#    clm_restart_file =  'clm.nc'

# &dart_to_clm_nml
#    dart_to_clm_input_file  = 'dart_posterior.nc'
#    dart_to_clm_output_file = 'clm_restart.nc'
#    repartition_vhist_file = 'clm_vector_history.nc'
#    repartition_analysis_file = 'dart_posterior_vector.nc'



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
# DART vector will continue to have an inflation value.

# Inflation files from a previous cycle have names in the CESM style, something like
#    ${case}.dart.rh.${scomp}_output_priorinf_{mean,sd}.YYYY-MM-DD-SSSSS.nc

# The script does NOT check the model date of the files for consistency
# with the current forecast time, so check that the inflation mean
# files are evolving as expected.

# CESM's st_archive should archive the inflation restart files
# like any other "restart history" (.rh.) files; copying the latest files
# to the archive directory, and moving all of the older ones.

#==========================================================================

# fill_inflation_restart creates files for all the domains in play,
# with names like input_priorinf_[mean,sd]_d0?.nc These should be renamed
# to be similar to what is created during the cycling. fill_inflation_restart
# only takes a second and only runs once.

if ( -e clm_inflation_cookie ) then
   # Set up template files for fill_inflation restart.
   set     LND_RESTART_FILENAME = ${CASE}.clm2_0001.r.${LND_DATE_EXT}.nc
   set     LND_HISTORY_FILENAME = ${CASE}.clm2_0001.h0.${LND_DATE_EXT}.nc
   set LND_VEC_HISTORY_FILENAME = ${CASE}.clm2_0001.h2.${LND_DATE_EXT}.nc
   
   # remove any potentially pre-existing links
   # These links are expected by fill_inflation_restart and by
   ${REMOVE} clm_restart.nc
   ${REMOVE} clm_history.nc
   ${REMOVE} clm_vector_history.nc
   
   ${LINK} ${LND_RESTART_FILENAME} clm_restart.nc || exit 4
   ${LINK} ${LND_HISTORY_FILENAME} clm_history.nc || exit 4
   if (  -s   ${LND_VEC_HISTORY_FILENAME} ) then
      ${LINK} ${LND_VEC_HISTORY_FILENAME} clm_vector_history.nc || exit 4
   endif

   ${EXEROOT}/fill_inflation_restart || exit  4

   foreach FILE ( input_priorinf_*.nc )
      set NEWBASE = `echo $FILE:r | sed -e "s#input#output#"`
      ${MOVE} ${FILE} clm_${NEWBASE}.1601-01-01-00000.nc
   end

   # Make sure this only happens once. Eat the cookie.
   ${REMOVE} clm_inflation_cookie

   # We have to potentially deal with files like:
   # clm_output_priorinf_mean_d01.${LND_DATE_EXT}.nc
   # clm_output_priorinf_mean_d02.${LND_DATE_EXT}.nc
   # clm_output_priorinf_mean_d03.${LND_DATE_EXT}.nc
   # clm_output_priorinf_sd_d01.${LND_DATE_EXT}.nc
   # clm_output_priorinf_sd_d02.${LND_DATE_EXT}.nc
   # clm_output_priorinf_sd_d03.${LND_DATE_EXT}.nc

   # To help keep track of the most recent inflation file,
   # create a 'pointer file' to hold the name of the most recent.

   @ domaincount = 0
   foreach FILE ( clm_output_priorinf_mean*.nc )

      @ domaincount ++

      set POINTERFILE = `printf clm_priorinf_pointer_d%02d.txt $domaincount`

      set SDFILE = `echo $FILE | sed -e "s#mean#sd#"`

      echo $FILE   >! $POINTERFILE
      echo $SDFILE >> $POINTERFILE

   end

   # Not supporting posterior inflation at this time.
   ${REMOVE} input_postinf*nc

endif


# Check to see if inflation is being used.

set  MYSTRING = `grep inf_flavor input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  PRIOR_INF = $MYSTRING[2]
set  POSTE_INF = $MYSTRING[3]

if ( $PRIOR_INF != 0 ) then

   # CLM always has at least two domains, but may sometimes has three.
   # Link to the new expected name, if the file does not exist, filter will
   # die and issue a very explicit death message.

   ${REMOVE} input_priorinf_mean*.nc input_priorinf_sd*.nc

   @ domaincount = 1

   foreach POINTERFILE ( clm_priorinf_pointer*.txt )

      set DOMAIN = `printf _d%02d $domaincount`
      set INPUT  =  input_priorinf_mean_${DOMAIN}
      set OUTPUT = output_priorinf_mean_${DOMAIN}

#  Has files clm_output_priorinf_mean*.nc  and ...sd...
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
# clm always needs a clm_restart.nc, clm_history.nc for geometry information, etc.
# it may or may not need a vector-format history file - depends on user input

${REMOVE} restart_files.txt history_files.txt vector_files.txt

# StrCpl; LND_DATE_EXT is extracted from the rpointer.lnd_0001,
#         which is linked here by multi-assim.csh.
#         CASE is 'CESM' for strongly-coupled cases, something else for all others.
foreach FILE ( ${CASE}.clm2_*.r.${LND_DATE_EXT}.nc )

   # create unique output filename for the copy that has the indeterminate
   # values replaced by the _FillValue. The copies are the files that will
   # be used to construct the DART vector.
   set OUTPUT = `echo $FILE | sed -e "s/${CASE}.//"`

# StrCpl; OUTPUT is redefined below here and clm.nc is used only here
#         why COPY?  Safety?
   ${COPY} $FILE clm.nc
   ${EXEROOT}/clm_to_dart >& /dev/null
   ${MOVE} clm.nc $OUTPUT

end

ls -1         clm2_*.r.${LND_DATE_EXT}.nc  >! restart_files.txt
ls -1 ${CASE}.clm2_*.h0.${LND_DATE_EXT}.nc >! history_files.txt
ls -1 ${CASE}.clm2_*.h2.${LND_DATE_EXT}.nc >! vector_files.txt

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
${REMOVE} clm_restart.nc
${REMOVE} dart_posterior.nc

# Identify if SWE re-partitioning is necessary
set  REPARTITION = `grep repartition_swe input.nml \
                    | sed -e "s/repartition_swe//g" -e "s/=//g"`

if ($REPARTITION != 0) then
#    `was unlink clm_vector_history`, which appears only as part of other names
   ${REMOVE} clm_vector_history.nc
   # Track the ensemble count to locate matching vector analysis file
   @ enscount = 1

   foreach RESTART ( ${CASE}.clm2_*.r.${LND_DATE_EXT}.nc )
    
     # ? Why remove $CASE?  What's left is not used to make other file names
     #   and it's not a DART file name.
     set POSTERIOR_RESTART = `echo $RESTART | sed -e "s/${CASE}.//"`
     # The vectors are the third domain(?)
     # This is a DART name format; analysis_member_####_d03.nc
     # ? But this assumes that stage 'analysis' is written?  It's not in the default input.nml.
     # ? Why not 'output'?
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
        
     # Make links named according to dart_to_clm_nml.
     # Target is a CESM file name without the CASE name.
     ${LINK} $POSTERIOR_RESTART  dart_posterior.nc
     # Target is not CESM format and maybe not DART (see above)
     ${LINK} $POSTERIOR_VECTOR  dart_posterior_vector.nc
     # Target is a CESM file name restart.
     ${LINK} $RESTART   clm_restart.nc
     # Target is a CESM h2 file name.
     ${LINK} $CLM_VECTOR   clm_vector_history.nc

     # Does dart_to_clm yield files with CESM names?
     # Yes, that is what the RESTART link does.
     ${EXEROOT}/dart_to_clm >& /dev/null

     if ($status != 0) then
        echo "ERROR: dart_to_clm failed for $RESTART"
        exit 7
     endif
     
     foreach LIST (clm_restart.nc clm_vector_history.nc dart_posterior.nc \
              dart_posterior_vector.nc )

             ${REMOVE} $LIST
     end
     @ enscount ++
  end


else

   foreach RESTART ( ${CASE}.clm2_*.r.${LND_DATE_EXT}.nc )
   
      set POSTERIOR = `echo $RESTART | sed -e "s/${CASE}.//"`
   
      ${LINK} $POSTERIOR dart_posterior.nc
      ${LINK} $RESTART   clm_restart.nc
   
      ${EXEROOT}/dart_to_clm >& /dev/null
   
      if ($status != 0) then
         echo "ERROR: dart_to_clm failed for $RESTART"
           exit 8
      endif
   
      ${REMOVE} dart_posterior.nc
      ${REMOVE} clm_restart.nc
   end

endif

# Remove the copies that we no longer need. The posterior values are
# in the DART diagnostic files for the appropriate 'stage'.
# StrCpl; this removes the local links, but not the target files in ../../CESM/run,
#         Why don't these use $REMOVE (rm -rf), which also doesn't remove soft link targets ?
${REMOVE} clm2_*.{r,h0,h2}.${LND_DATE_EXT}.nc 

#=========================================================================
# Block 8: Rename the output using the CESM file-naming convention.
# and prepare inflation pointer files for next cycle.
#=========================================================================

# The short-term archiver archives files depending on pieces of their names.
# '.dart.i.' files are ensemble statistics (mean, sd) of just the state variables 
#            in the initial files.
# '.e.'      designates a file as something from the 'external system processing ESP', e.g. DART.
# '_d##'     The domain label is not part of the CESM naming, so keep it with the DART file type.

# stages_all not available from CESM's DART_params yet.
# DART 'output_member_****.nc' files are actually linked to cam input files

# foreach FILE (`${LIST} ${stages_all}_member_*.nc`)
foreach FILE (`${LIST} *_member_*.nc`)

   set list = `echo $FILE:r | sed -e "s#_# #g"`
   # DART file type doesn't need 'member' and inst number is moved to $scomp.
   # Keep the stage and domain.
   set dart_file = `echo $list[1] $list[$#list] | sed -e "s# #_#g"`
   @ elem = $#list - 2
   set inst = $list[$elem] 
   set type = "e"

   ${MOVE} $FILE \
       ${hindcast_path}${CASE}.${tor_obs}_${inst}.${type}.${dart_file}.${LND_DATE_EXT}.nc || exit 150
end

# Files without instance numbers need to have the scomp part of their names = "dart".
# This is because in st_archive, all files with  scomp = "cam"
# (= compname in env_archive.xml) will be st_archived using a pattern
# which has the instance number added onto it.  {mean,sd} files don't have 
# instance numbers, so they need to be archived by the "dart" section of env_archive.xml.
# But they still need to be different for each component, so include $scomp (tor_obs) in the
# ".dart_file" part of the file name.  Somewhat awkward and inconsistent, but effective.

# Means and standard deviation files (except for inflation).
foreach FILE (`${LIST} ${stages_all}_{mean,sd}*.nc`)

   set parts = $FILE:r 
   set type = "e"

   ${MOVE} $FILE \
      ${hindcast_path}${CASE}.dart.${type}.${tor_obs}_${parts}.${LND_DATE_EXT}.nc || exit 160
end

set obs_seq_final = "${CASE}.dart.e.${tor_obs}_obs_seq_final.${LND_DATE_EXT}"
${MOVE} obs_seq.final ${hindcast_path}$obs_seq_final || exit 170
${LINK} ${hindcast_path}$obs_seq_final oXTerm2024-02-25.06:56:35bs_seq.final

${MOVE} dart_log.out ${hindcast_path}${tor_obs}_dart_log.${LND_DATE_EXT}.out || exit 171
gzip input.nml
${MOVE} input.nml.gz ${hindcast_path}${tor_obs}input.nml.log.${JOBID}.gz


# Rename the inflation files and designate them as 'rh' files - which get
# reinstated in the run directory by the short-term archiver and are then
# available for the next assimilation cycle.
#        
# Accommodate any possible inflation files, including domains (clm).
# The .${scomp}_ part is needed by DART to distinguish
# between inflation files from separate components in coupled assims.
   
# StrCpl; for now inflation is only applied to self-assims and 'no coupling' cases.
#         So the existing scripting should be sufficient.
#         But don't manage the cross-assim phase files.

if ($cpld =~ 'self') then  
   echo "Renaming inflation files and updating inflation pointer files."
   
   # stages_all not available from CESM's DART_params yet.
   # foreach FILE (`${LIST} ${stages_all}_{prior,post}inf_*`)
   @ domaincount = 0
   foreach FILE (`${LIST} *_{prior,post}inf_mean*`)
      @ domaincount ++
      set POINTERFILE = `printf clm_priorinf_pointer_d%02d.txt $domaincount`

      set DEST = "${hindcast_path}${CASE}.dart.rh.$FILE:r.${LND_DATE_EXT}.nc "
      ${MOVE} $FILE $DEST || exit 180
      ${LINK} $DEST .
      
      set SDFILE = `echo $DEST | sed -e "s#mean#sd#"`
      echo $DEST   >! $POINTERFILE
      echo $SDFILE >> $POINTERFILE
   end

endif 

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

echo "`date` -- END CLM_ASSIMILATE"

exit 0

