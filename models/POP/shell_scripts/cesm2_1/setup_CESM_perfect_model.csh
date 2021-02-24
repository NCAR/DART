#!/bin/csh -f
#
# Copyright 2020 University Corporation for Atmospheric Research
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License.
# Please view the License at http://www.apache.org/licenses/LICENSE-2.0
#
# ==============================================================================
#
# ---------------------
# Purpose
# ---------------------
#
# This script is designed to set up, stage, and build a single-instance run
# of CESM using a G compset where only POP is active and the atm and land
# states are specified by data files. The initial state can come from a
# single multi-instance reference case so a CESM hybrid setup is used.
#
# This script has a counterpart that is a multi-instance setup for either a
# free run or an assimilation experiment. To make it easy to maintain (and
# hopefully understand), the two scripts are intended to parallel each other.
# That means this script performs a lot of manipulation of the 'instance'
# portion of the filenames, which seems unnecessary initially.
#
# This script results in a viable setup for a CESM single instance experiment.
# You are STRONGLY encouraged to run the single instance CESM a few times and
# experiment with different settings BEFORE you try to generate 'perfect'
# observations. You should become comfortable using CESM's restart capability
# to re-stage files in your RUN directory.
#
# ${CASEROOT}/CESM_DART_config.csh will change the xml settings of the CESM
# case to enable data assimilation and set a data assimilation script, in
# this case, perfect_model.csh, which configures DART to harvest synthetic
# observations.
#
# This, and the required setup, can be run at a later date. e.g. you can
# advance an ensemble from 2004-01-01 to 2004-02-01 and then run
# CESM_DART_config.csh to augment the existing run script, modify STOP_N
# to 6 hours, and start harvesting synthetic observations when CESM stops
# at 2004-02-01 06Z.
#
# This script relies heavily on the information in:
# http://esmci.github.io/cime/versions/master/html/index.html
#
# ---------------------
# How to use this script.
# ---------------------
#
# -- You will have to read and understand the script in its entirety.
#    You will have to modify things outside this script.
#    This script sets up a plain CESM single-instance run without DART,
#    intentionally.  Once it is running, calls to DART can be added.
#
# -- Examine the whole script to identify things to change for your experiments.
#
# -- Edit this script in the $DART/models/POP/shell_scripts/cesm2_1 directory
#    or copy it to somewhere where it will be preserved.
#
# -- Locate the initial condition files that CESM will need.
#
# -- Run this script. When it is executed, it will create:
#    1) a CESM 'CASE' directory, where the model will be built,
#    2) a run directory, where each forecast will take place,
#    3) a bld directory for the executables.
#    4) The short term archiver will use a fourth directory for
#       storage of model output until it can be moved to long term storage
#       (campaign storage or the research data archive)
#
# -- To run DART; read, understand, and execute ${CASEROOT}/CESM_DART_config.csh
#
# -- Submit the job using ${CASEROOT}/case.submit -M begin,end
#
# If you want to change something in your case other than the runtime
# settings, it is safest to delete everything and start the run from scratch.
#
# ==============================================================================

# ==============================================================================
# Source DART_params.csh.
#
# Many of the environmental variables are stored in the DART_params.csh file
# and are sourced by the various scripts. This increases the ease of maintenance
# and decreases the errors associated with inconsistent script settings.
# ==============================================================================

echo ""
echo "Source DART_params.csh."
echo ""

source ./DART_params.csh || exit 1

# ==============================================================================
# Set CASE to the combination of the compset, resolution and an optional extra
# string.
#
# extra_string  An optional extra_string to differentiate a case from others
#               that have the same compset and resolution.
# CASE          The value of "case" will be used many ways; directory and file
#               names both locally and on HPSS, and script names; so consider
#               its length and information content.
# ==============================================================================

setenv extra_string ''
setenv CASE         ${case_string}.${resolution}${extra_string}

echo ""
echo "Case name is set to ${CASE}"
echo ""

# ==============================================================================
# Directories:
#
# CASEROOT        Will create the CESM case directory here, where the CESM+DART
#                 configuration files will be stored.  This should probably not
#                 be in scratch (on GLADE, your 'work' partition is suggested).
#                 This script will delete any existing CASEROOT, so this script,
#                 and other useful things should be kept elsewhere.
# RUNDIR          Will create the CESM run directory here.  Will need large
#                 amounts of disk space, generally on a scratch partition.
# EXEROOT         Will create the CESM executable directory here, where the
#                 CESM executables will be built.  Medium amount of space
#                 needed, generally on a scratch partition.
# ARCHDIR         Will create the CESM short-term archive directories here.
#                 Large, generally on a scratch partition.  Files will remain
#                 here until the long-term archiver moves it to permanent storage.
# ==============================================================================

setenv CASEROOT     /glade/work/${USER}/cases/${CASE}
setenv RUNDIR       /glade/scratch/${USER}/${CASE}/run
setenv EXEROOT      /glade/scratch/${USER}/${CASE}/bld
setenv ARCHDIR      /glade/scratch/${USER}/archive/${CASE}

# ==============================================================================
# Create the case - this creates the CASEROOT directory.
#
# For list of the pre-defined component sets: ./create_newcase --list
# To create a variant compset, see the CESM documentation and carefully
# incorporate any needed changes into this script.
#
# The script won't allow the user to make CASEROOT the same dir as where this
# setup script is since the build process removes all files in the ${CASEROOT}
# dir before populating it.
# ==============================================================================

if ( $CASEROOT == `dirname $0` ) then
   echo "ERROR: the setup script should not be located in the CASEROOT"
   echo "directory, because all files in the CASEROOT dir will be removed"
   echo "before creating the new case.  move the script to a safer place."
   exit -1
endif

echo "removing old files from ${CASEROOT}"
echo "removing old files from ${EXEROOT}"
echo "removing old files from ${RUNDIR}"
${REMOVE} ${CASEROOT}
${REMOVE} ${EXEROOT}
${REMOVE} ${RUNDIR}

# ==============================================================================
# Invoking create_newcase.
# The --run-unsupported option is necessary for CESM2 because this compset
# hasn't been tested.
# ==============================================================================

set use_tasks_per_node = 36
set nthreads = 1

echo ""
echo "Invoking create_newcase."
echo ""
echo "${cesmroot}/cime/scripts/create_newcase"
echo "--case ${CASEROOT}"
echo "--mach ${mach}"
echo "--res ${resolution}"
echo "--compset ${compset}"
echo "--run-unsupported"

${cesmroot}/cime/scripts/create_newcase \
--case ${CASEROOT} \
--mach ${mach} \
--res ${resolution} \
--compset ${compset} \
--run-unsupported

if ( $status != 0 ) then
   echo "ERROR: Case could not be created."
   exit -1
endif

# preserve a copy of this script as it was run
set ThisFileName = $0:t
${COPY} $ThisFileName ${CASEROOT}/${ThisFileName}.original
${COPY} DART_params.csh ${CASEROOT}

# ==============================================================================
# CESM_DART_config.csh can be run at some later date if desired, but it presumes
# to be run from a CASEROOT directory. If CESM_DART_config.csh does not exist
# locally, then it better exist in the expected part of the DARTROOT tree.
# ==============================================================================

if ( ! -e CESM_DART_config.csh ) then
   ${COPY} ${DARTROOT}/models/POP/shell_scripts/${cesmtagmajor}/CESM_DART_config.csh .
endif

if (  -e   CESM_DART_config.csh ) then
   ${COPY} CESM_DART_config.csh ${CASEROOT}
else
   echo "WARNING: the script to configure for data assimilation is not available."
   echo "         CESM_DART_config.csh should be present locally or in"
   echo "         ${DARTROOT}/models/POP/shell_scripts/${cesmtagmajor}/"
   echo "         You can stage this script later if you must."
endif

# ==============================================================================
# Configuring the case.
# ==============================================================================

echo ""
echo "Configuring the case."
echo ""

cd ${CASEROOT}

foreach FILE ( *xml )
   if ( ! -e          ${FILE}.original ) then
      ${COPY} ${FILE} ${FILE}.original
   endif
end

# Make sure the land and atmosphere are 'data' components.

setenv COMP_ATM  `./xmlquery COMP_ATM --value`
setenv COMP_LND  `./xmlquery COMP_LND --value`

if ( (${COMP_ATM} != datm) || (${COMP_LND} != slnd) ) then
   echo ""
   echo "ERROR: This setup script is not appropriate for active land or atmosphere compsets."
   echo "ERROR: Please use the models/CESM/shell_scripts examples for these cases."
   echo ""
   exit -3
endif

# This configures each POP model instance to run on a single node. 
# Given the resolution and the machine architecture, you may need to use 
# more nodes per instance.  In the multi-instance configuration, getting the 
# optimal layout is quite problematic. Running the components sequentially 
# has been show to be robust, so we are setting the ROOTPEs to zero 
# (which makes them run sequentially).

@ number_of_threads = 1
@ nodes_per_instance = 1
@ ntasks_active = -1 * $nodes_per_instance
@ ntasks_data   = -1

./xmlchange ROOTPE_ATM=0,NTHRDS_ATM=$nthreads,NTASKS_ATM=$ntasks_active,NINST_ATM=1
./xmlchange ROOTPE_LND=0,NTHRDS_LND=$nthreads,NTASKS_LND=$ntasks_active,NINST_LND=1
./xmlchange ROOTPE_ICE=0,NTHRDS_ICE=$nthreads,NTASKS_ICE=$ntasks_active,NINST_ICE=1
./xmlchange ROOTPE_ROF=0,NTHRDS_ROF=$nthreads,NTASKS_ROF=$ntasks_active,NINST_ROF=1
./xmlchange ROOTPE_OCN=0,NTHRDS_OCN=$nthreads,NTASKS_OCN=$ntasks_active,NINST_OCN=1
./xmlchange ROOTPE_GLC=0,NTHRDS_GLC=$nthreads,NTASKS_GLC=$ntasks_active,NINST_GLC=1
./xmlchange ROOTPE_WAV=0,NTHRDS_WAV=$nthreads,NTASKS_WAV=$ntasks_active,NINST_WAV=1
./xmlchange ROOTPE_CPL=0,NTHRDS_CPL=$nthreads,NTASKS_CPL=$ntasks_active
./xmlchange ROOTPE_ESP=0,NTHRDS_ESP=$nthreads,NTASKS_ESP=$ntasks_data

# ==============================================================================
# Definition of a HYBRID run:
#
# https://esmci.github.io/cime/versions/master/html/users_guide/running-a-case.html
# "A hybrid run indicates that CESM is initialized more like a startup, but uses
# initialization datasets from a previous case. This is somewhat analogous to a
# branch run with relaxed restart constraints. A hybrid run allows users to bring
# together combinations of initial/restart files from a previous case (specified
# by $RUN_REFCASE) at a given model output date (specified by $RUN_REFDATE).
# Unlike a branch run, the starting date of a hybrid run (specified by $RUN_STARTDATE)
# can be modified relative to the reference case. In a hybrid run, the model does not
# continue in a bit-for-bit fashion with respect to the reference case. The resulting
# climate, however, should be continuous provided that no model source code or
# namelists are changed in the hybrid run. In a hybrid initialization, the ocean
# model does not start until the second ocean coupling (normally the second day),
# and the coupler does a "cold start" without a restart file."
#
# A hybrid start is better for POP because the velocities are used rather than
# just T,S. A hybrid start is also more desirable because initial values can be
# specified for ROF - as opposed to just zeros.
# ==============================================================================

./xmlchange RUN_TYPE=hybrid
./xmlchange RUN_STARTDATE=${start_year}-${start_month}-${start_day}
./xmlchange START_TOD=$start_tod
./xmlchange RUN_REFCASE=$refcase
./xmlchange RUN_REFDATE=$refdate
./xmlchange RUN_REFTOD=$reftod
./xmlchange GET_REFCASE=FALSE
./xmlchange EXEROOT=${EXEROOT}
./xmlchange RUNDIR=${RUNDIR}

./xmlchange DATM_MODE=CPLHIST
./xmlchange DATM_CPLHIST_CASE=${CASE}
./xmlchange DATM_CPLHIST_YR_ALIGN=$stream_year_align
./xmlchange DATM_CPLHIST_YR_START=$stream_year_first
./xmlchange DATM_CPLHIST_YR_END=$stream_year_last

# The default for G compsets, e.g. those using the COREv2 and JRA-55
# streamfiles, is to have NO_LEAP calendars by default while the CAM
# reanalysis has a GREGORIAN calendar, so we change it here.
./xmlchange CALENDAR=GREGORIAN

./xmlchange STOP_OPTION=$stop_option
./xmlchange STOP_N=$first_STOP_N
./xmlchange CONTINUE_RUN=FALSE
./xmlchange RESUBMIT=$resubmit

./xmlchange PIO_TYPENAME=pnetcdf

./xmlchange GLC_NCPL=1
./xmlchange OCN_NCPL=1

# These are archiving options that may be used.
# You can turn the short/long term archivers on or off,
# but these settings should be made in either event.

./xmlchange DOUT_S_ROOT=${ARCHDIR}
./xmlchange DOUT_S_SAVE_INTERIM_RESTART_FILES=TRUE

if ($short_term_archiver == 'off') then
   ./xmlchange DOUT_S=FALSE
else
   ./xmlchange DOUT_S=TRUE
endif

# level of debug output, 0=minimum, 1=normal, 2=more, 3=too much, valid values: 0,1,2,3 (integer)

./xmlchange DEBUG=FALSE
./xmlchange INFO_DBUG=0

# ==============================================================================
# Updating source files.
#
# Ideally, using DART would not require any modifications to the model source.
# Until then, this script accesses sourcemods from a hardwired location.
# If you have additional sourcemods, they will need to be merged into any DART
# mods and put in the SourceMods subdirectory found in the 'CASEROOT' directory.
# ==============================================================================

echo ""
echo "Updating source files."
echo ""

if (    -d     ~/${cesmtag}/SourceMods ) then
   ${COPY} -r  ~/${cesmtag}/SourceMods/* ${CASEROOT}/SourceMods/
else
   echo "ERROR - No SourceMods for this case."
   echo "DART requires modifications to several src files that cause POP to"
   echo "recompute the barotropic velocity after assimilation."
   echo "Download the appropriate files for CESM2_1 from:"
   echo "http://www.image.ucar.edu/pub/DART/CESM"
   echo "untar these into your HOME directory - they will create a"
   echo "~/${cesmtag}  directory with the appropriate SourceMods structure."
   exit -4
endif

# ==============================================================================
# Updating the run script to reflect queue and wallclock
# ==============================================================================

echo ""
echo "Updating the run script to set wallclock and queue."
echo ""

./xmlchange --subgroup case.st_archive --id JOB_WALLCLOCK_TIME --val $timewall
./xmlchange --subgroup case.run --id USER_REQUESTED_WALLTIME --val $timewall
./xmlchange --subgroup case.run --id JOB_WALLCLOCK_TIME --val $timewall

# ==============================================================================
# Setting up the case.
#
# This creates the EXEROOT and RUNDIR directories.
# ==============================================================================

echo ""
echo "Setting up the case ..."
echo "This creates the EXEROOT and RUNDIR directories."
echo ""

./case.setup

if ( $status != 0 ) then
   echo "ERROR: Case could not be set up."
   exit -2
endif

# ==============================================================================
# Modifying namelist templates for this instance.
# 
# Even though there is only one instance, we use a loop to keep the pmo and
# multiple instance scripts as similar as possible.
# ==============================================================================

echo ""
echo "Modifying namelist templates for this instance."
echo ""

@ inst = 1
while ($inst <= 1)

   # following the CESM strategy for 'inst_string'
   set inst_string = ''

   # ===========================================================================
   set fname = "user_nl_datm${inst_string}"
   # ===========================================================================
   # DATM Namelist

   echo "dtlimit  = 1.5, 1.5"               >> ${fname}
   echo "fillalgo = 'nn', 'nn'"             >> ${fname}
   echo "fillmask = 'nomask','nomask'"      >> ${fname}
   echo "mapalgo  = 'bilinear','bilinear'"  >> ${fname}
   echo "mapmask  = 'nomask','nomask'"      >> ${fname}
   echo "streams  = 'datm.streams.txt.CPLHISTForcing.nonSolarFlux${inst_string} $stream_year_align $stream_year_first $stream_year_last'," >> ${fname}
   echo "           'datm.streams.txt.CPLHISTForcing.Solar${inst_string}        $stream_year_align $stream_year_first $stream_year_last',"  >> ${fname}
   echo "           'datm.streams.txt.CPLHISTForcing.State1hr${inst_string}     $stream_year_align $stream_year_first $stream_year_last',"  >> ${fname}
   echo "           'datm.streams.txt.CPLHISTForcing.State3hr${inst_string}     $stream_year_align $stream_year_first $stream_year_last'"  >> ${fname}
   echo "taxmode  = 'cycle','cycle'"        >> ${fname}
   echo "tintalgo = 'linear','linear'"      >> ${fname}
   echo "restfils = 'unset'"                >> ${fname}
   echo "restfilm = 'unset'"                >> ${fname}

   # ===========================================================================
   set fname = "user_nl_cice${inst_string}"
   # ===========================================================================
   # CICE Namelist
   # Note well: this is true as of CESM2_1 but may change in the future.
   # In a hybrid case, setting restart files for CICE is a distinct process from
   # setting restart files for POP. The restart file is set in the user_nl_cice
   # file in the ${CASEROOT} directory using the ice_ic (ice initial condition)
   # key in the namelist file instead of using an rpointer file in ${RUNDIR}.

   if ( $SingleInstanceRefcase ) then
      set true_string = ''
   else
      set true_string = `printf _%04d $TRUTHinstance`
   endif

   echo "ice_ic = '${stagedir}/${refcase}.cice${true_string}.r.${reftimestamp}.nc'" >> ${fname}

   # ===========================================================================
   set fname = "user_nl_pop${inst_string}"
   # ===========================================================================
   # POP Namelist
   # init_ts_suboption = 'data_assim'   for non bit-for-bit restarting (assimilation mode)
   # init_ts_suboption = 'rest'         --> default behavior
   #
   # README:
   # Configuring the contents of the history files for POP is best explained in
   # the cesm2 pop2 namelist documentation.
   # http://www.cesm.ucar.edu/models/cesm2/settings/current/pop2_nml.html
   #
   # and the CESM-specific documentation for the tavg output variables in the pop2
   # online documentation:
   # https://ncar.github.io/POP/doc/build/html/users_guide/model-diagnostics-and-output.html#time-averaged-history-files

   echo "init_ts_suboption  = 'data_assim'" >> ${fname}

   @ inst ++
end

# ==============================================================================
# Staging restart files now that the run directory exists
# ==============================================================================

echo ""
echo "Staging restart files now that the run directory exists"
echo ""

set init_time = ${reftimestamp}

cat << EndOfText >! stage_cesm_files
#!/bin/csh -f
# This script can be used to help restart an experiment from any previous step.
# The appropriate files are copied to the RUN directory.
#
# Before running this script:
#  1) be sure CONTINUE_RUN is set correctly in the env_run.xml file in
#     your CASEROOT directory.
#     CONTINUE_RUN=FALSE => you are starting over at the initial time.
#     CONTINUE_RUN=TRUE  => you are starting from a previous step but not
#                           the very first one.
#  2) be sure 'restart_time' is set to the day and time that you want to
#     restart from if not the initial time.

set restart_time = $init_time

# get the settings for this case from the CESM environment

cd ${CASEROOT}

setenv RUNDIR           \`./xmlquery RUNDIR           --value\`
setenv CONTINUE_RUN     \`./xmlquery CONTINUE_RUN     --value\`
setenv ARCHIVED         \`./xmlquery DOUT_S           --value\`
setenv ARCHIVE_DIR      \`./xmlquery DOUT_S_ROOT      --value\`

cd \${RUNDIR}

echo "Copying the required CESM files to the run directory to rerun"
echo "a previous step.  CONTINUE_RUN from env_run.xml is" \${CONTINUE_RUN}
if ( \${CONTINUE_RUN} == TRUE ) then
  echo "so files for some later step than the initial one will be restaged."
  echo "Date to reset files to is: \${restart_time}"
else
  echo "so files for the initial step of this experiment will be restaged."
  echo "Date to reset files to is: ${init_time}"
endif
echo ""

if ( \${CONTINUE_RUN} == TRUE ) then

   #----------------------------------------------------------------------
   # This block copies over a set of restart files from any previous step of
   # the experiment that is NOT the initial step.
   # After running this script resubmit the job to rerun.
   #----------------------------------------------------------------------

   echo "Staging restart files for run date/time: " \${restart_time}

   #  The short term archiver is on, so the files we want should be in one
   #  of the short term archive 'rest' restart directories.  This assumes
   #  the long term archiver has NOT copied these files to the HPSS yet.

   if (  \${ARCHIVED} == TRUE ) then

      # The restarts should be in the short term archive directory.  See
      # https://esmci.github.io/cime/versions/master/html/users_guide/running-a-case.html
      # for more help and information.

      set RESTARTDIR = \${ARCHIVE_DIR}/rest/\${restart_time}

      if ( ! -d \${RESTARTDIR} ) then

         echo ""
         echo "ERROR: restart file directory not found: "
         echo " \${RESTARTDIR}"
         echo "If the long-term archiver is on, you may have to restore this directory first."
         echo "You can also check for either a .sta or a .sta2 hidden subdirectory in"
         echo "\${ARCHIVE_DIR}"
         echo "which may contain the 'rest' directory you need,"
         echo "and then modify RESTARTDIR in this script."
         exit -1

      endif

      ${COPY} \${RESTARTDIR}/* . || exit -1

   else

      # The short term archiver is off, which leaves all the restart files
      # in the run directory.  The rpointer files must still be updated to
      # point to the files with the right day/time.

      @ inst=1
      while (\$inst <= 1)

         set inst_string = ''

         echo "${CASE}.datm\${inst_string}.r.\${restart_time}.nc"    >! rpointer.atm\${inst_string}
         echo "${CASE}.datm\${inst_string}.rs1.\${restart_time}.bin" >> rpointer.atm\${inst_string}
         echo "${CASE}.drof\${inst_string}.r.\${restart_time}.nc"    >! rpointer.rof\${inst_string}
         echo "${CASE}.drof\${inst_string}.rs1.\${restart_time}.bin" >> rpointer.rof\${inst_string}
         echo "${CASE}.cice\${inst_string}.r.\${restart_time}.nc"    >! rpointer.ice\${inst_string}

         echo "${CASE}.pop\${inst_string}.ro.\${restart_time}"    >! rpointer.ocn\${inst_string}.ovf
         echo "${CASE}.pop\${inst_string}.r.\${restart_time}.nc"  >! rpointer.ocn\${inst_string}.restart
         echo "RESTART_FMT=nc"                                    >> rpointer.ocn\${inst_string}.restart

         if ( -e rpointer.ocn\${inst_string}.tavg ) then
            echo "${CASE}.pop\${inst_string}.rh.\${restart_time}.nc" >! rpointer.ocn\${inst_string}.tavg
         endif
         if ( -e rpointer.ocn\${inst_string}.tavg.2 ) then
            echo "${CASE}.pop\${inst_string}.rh.nday1.\${restart_time}.nc" >! rpointer.ocn\${inst_string}.tavg.2
         endif

         @ inst ++
      end

      # The coupler file has no instance string.
      echo "${CASE}.cpl.r.\${restart_time}.nc" >! rpointer.drv

   endif

   echo "All files reset to rerun experiment step for time " \$restart_time

else     # CONTINUE_RUN == FALSE

   #----------------------------------------------------------------------
   # This block links the right files to rerun the initial (very first)
   # step of an experiment.  The names and locations are set during the
   # building of the case; to change them rebuild the case.
   # After running this script resubmit the job to rerun.
   #----------------------------------------------------------------------

   @ inst=1
   while (\$inst <= 1)

      set inst_string = ''

      if ( $SingleInstanceRefcase ) then
         set true_string = ''
         echo "Staging initial files."
      else
         set true_string = `printf _%04d $TRUTHinstance`
         echo "Staging initial files from instance $TRUTHinstance for the truth run."
      endif

      echo "${stagedir}/${refcase}.pop\${true_string}.ro.${init_time}"   >! rpointer.ocn\${inst_string}.ovf
      echo "${stagedir}/${refcase}.pop\${true_string}.r.${init_time}.nc" >! rpointer.ocn\${inst_string}.restart
      echo "RESTART_FMT=nc"                                              >> rpointer.ocn\${inst_string}.restart

      @ inst ++
   end

   echo "All files set to run the FIRST experiment step at time" $init_time

endif
exit 0

EndOfText
chmod 0755 stage_cesm_files

./stage_cesm_files

# ==============================================================================
#
# To create custom stream files we:'
# 1. Use preview_namelists to obtain the contents of the stream txt files
#    in CaseDocs, and then place a copy of the modified stream txt file in
#    ${CASEROOT} with the string user_ prepended, and
# 2. Copy a template stream txt file from this directory:
#    ${DARTROOT}/models/POP/shell_scripts/${cesmtagmajor}
#     and modify one for each instance.
#
# ==============================================================================

echo ""
echo "To create custom stream files we:"
echo "1. Use preview_namelists to obtain the contents of the stream txt files"
echo "   in CaseDocs, and then place a copy of the modified stream txt file in"
echo "   ${CASEROOT} with the string user_ prepended, and"
echo "2. Copy a template stream txt file from this directory:"
echo "   ${DARTROOT}/models/POP/shell_scripts/${cesmtagmajor}"
echo "   and modify one for each instance."
echo ""

./preview_namelists || exit -3

# This gives us a stream txt file for each instance that we can
# modify for our own purpose.

foreach FILE (CaseDocs/*streams*)
   set FNAME = $FILE:t

   switch ( ${FNAME} )
      case *presaero*:
         echo "Using default prescribed aerosol stream.txt file ${FNAME}"
         breaksw
      case *diatren*:
         echo "Using default runoff stream.txt file ${FNAME}"
         breaksw
      case *\.Precip_*:
         echo "Precipitation in user_datm.streams.txt.CPLHISTForcing.nonSolarFlux - not ${FNAME}"
         breaksw
      default:
         ${COPY} $FILE user_${FNAME}
         chmod   644   user_${FNAME}
         breaksw
   endsw

end

# Replace each default stream txt file with one that uses the CAM DATM
# conditions for a default year and modify the instance number.
# The stream files for POP have no leading zeros in the instance number.

echo ""
echo "Replacing each default stream txt file with one that uses CAM DATM"
echo ""

foreach FNAME (user*streams*)
   set name_parse = `echo ${FNAME} | sed 's/\_/ /g'`
   @ filename_index = $#name_parse
   set streamname = $name_parse[$filename_index]
   if ( $SingleInstanceRefcase ) then
      set instance = ''
   else
      set instance = `printf %04d $TRUTHinstance`
   endif

   if (-e $DARTROOT/models/POP/shell_scripts/$cesmtagmajor/user_$streamname*template) then

      echo "Copying DART template for ${FNAME} and changing instance."

      ${COPY} $DARTROOT/models/POP/shell_scripts/$cesmtagmajor/user_$streamname*template ${FNAME}

      sed s/NINST/$instance/g ${FNAME} >! out.$$
      ${MOVE} out.$$ ${FNAME}

   else
      echo "DIED Looking for a DART stream txt template for ${FNAME}"
      echo "DIED Looking for a DART stream txt template for ${FNAME}"
      exit -3
   endif

end

# ==============================================================================
# Building the case.
# ==============================================================================

echo ""
echo "Building the case."
echo ""

${BUILD_WRAPPER} ./case.build

if ( $status != 0 ) then
   echo "ERROR: Case could not be built."
   exit -5
endif

# ==============================================================================
# Checking the case.
# ==============================================================================

cat << EndOfText >! CESM_instructions.txt

-------------------------------------------------------------------------

Checking the case.

1) cd ${RUNDIR}
   and check the compatibility between the namelists/pointer
   files and the files that were staged.

2) cd ${CASEROOT}

3) The case is initially configured to NOT INVOKE ANY DART CODE.
   When you are ready to generate synthetic observations, configure and execute
   the ${CASEROOT}/CESM_DART_config.csh script.

4) The very first CESM advance (i.e. CONTINUE_RUN=FALSE)
   STOP_N must be longer than *AT LEAST 2 TIMES* the coupling
   frequency between the atmosphere and ocean.
   If coupling once a day, the first advance MUST be at least 48 hours.
   If coupling 4 times a day, the first advance MUST be at least 12 hours.
   After that, STOP_N can be as short as a single coupling frequency.

5) Verify the contents of env_run.xml and submit the CESM job:
   ./case.submit -M begin,end

6) After the job has run, check to make sure it worked and that
   a: POP is creating netCDF restart files,
   b: the right restart files exist in the run directory,
   c: (if you're running DART) the archive dart/hist directory has the DART output,
   d: everything is working correctly ...

7) To extend the run in $STOP_N "$stop_option" steps,
   change the env_run.xml variables:

   ./xmlchange CONTINUE_RUN=TRUE
   ./xmlchange RESUBMIT=<number_of_cycles_to_run>
   ./xmlchange STOP_N=$STOP_N

   and submit:

   ./case.submit -M begin,end

Check the streams listed in the streams text files.  If more or different
dates need to be added, change the $CASEROOT/user_*files*
and invoke './preview_namelists' so you can check the information in the
${CASEROOT}/CaseDocs or
${RUNDIR} directories.

-------------------------------------------------------------------------

EndOfText

cat CESM_instructions.txt

exit 0
