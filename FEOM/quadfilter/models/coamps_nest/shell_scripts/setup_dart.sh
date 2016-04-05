#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################
#
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Initialize the DART system to run experiments with the Navy's COAMPS
# model.  This takes the following actions: 
#
# 1. Create the DART experiment folder.  This is DART_WORK_DIR which
#    yields the following directory structure
#       DART_WORK_BASE
#          +-- 00001
#          +-- 00002
#          +-- ...
#    where the DART_WORK_DIR houses links to the executable files, the
#    customized template run scripts, along with some of the basic
#    files that DART needs to run the executables.
#
# 2. Copy data files that the binaries need to properly run from both 
#    the COAMPS directory and the DART template directory.
#
# 3. Copy the shell scripts needed to run the ensemble/model and
#    parallel assimilation and update the paths in them to point to
#    the current working directory structure.
# 
# 4. Build the DART binaries
#
# 5. Link the DART binaries to the DART working directory.
#
# The -c option to the script specifies the location of the path
# configuration file that contains all the setup information
######

######
# BEGIN DEFINITIONS
######

# Function to check the command return status 
function check_return
{
    if [ $1 -ne 0 ]; then
    echo "Call failed - aborting script"
    echo "Error in $2"
    exit $1
    fi
}

# Function to make sure that I'm leaving nothing undefined in the 
# paths.config file that needs to be there.  Note that this does
# *NOT* constrain the value to make sense (i.e. paths existing,
# proper permissions, etc.) - it only ensures that the variable
# has a defined value
function validate_paths_config
{
    validate_var COAMPS_HOME
    validate_var COAMPS_UTIL_HOME
    validate_var DART_HOME
    validate_var COAMPS_DATA
    validate_var DART_ARCH
    validate_var DART_BASE
    validate_var DART_NAMELIST
    validate_var COAMPS_FILES
    validate_var PERFECT_FILES
    validate_var STATE_DEF_FILE
    validate_var ENSEMBLE_SIZE
    validate_var CONCURRENT_RUNS
    validate_var PBS_QUEUE_ENS
    validate_var PBS_QUEUE_TRANS
    validate_var PBS_QUEUE_FILTER
    validate_var FILTER_PROCS
    validate_var FLAT_FILES
#    validate_var ARCMACH
}

# Function customizes the path config file name
function customize_path_config
{
	echo "  ${1}         - PATH_CONFIG=${ABS_CONFIG_FILE}"
	if [ ${1} == 'navdas_preproc_obs.sh' ]; then 
      ${SED} "s|CONFIG\/FILE\/GOES\/HERE|${ABS_CONFIG_FILE}|" ${DART_SCRIPTS_WORK}/${1}
    else
      ${SED} "s|^\(.*\bPATH_CONFIG\s*=\s*\).*|\1${ABS_CONFIG_FILE}|" ${DART_SCRIPTS_WORK}/${1}
    fi
}


# Helper function - see if a variable is defined and abort if it
# is not.  Use indirect referencing to take in the variable name
# as a string and access its value
function validate_var
{
    var_name=$1
    if [[ -z ${!var_name} ]]
    then
        echo "Variable ${var_name} is undefined - aborting script"
        exit 1
    fi
}

# Command Aliases
RM='rm -f -r'
CP='cp -f -p'
LN='ln -sf'
PERL='perl -i -p -e'
SED='sed -i -e'

######
# END DEFINITIONS
######

usage="Usage: `basename $0` -c PATH_CONFIG -b BUILD_DART -s SUBMIT_JOB [t|f] -h"

if [ $# -eq 0 ]
then
    echo $usage
    exit 1
fi

# Parse the options and grab the file name
while getopts ":c:b:s:" option
do
  case $option in
      c  ) PATH_CONFIG_FILE=$OPTARG;;
      b  ) build_dart=$OPTARG;;
      s  ) submit_jobs=$OPTARG;;
      \? ) echo $usage
           exit 65;;
      *  ) echo $usage
       exit 65;;
  esac
done

# Make sure the configuration file exists
if [ ! -e ${PATH_CONFIG_FILE} ]; then
    echo "Configuration file ${PATH_CONFIG_FILE} not found"
    exit 1
fi

# Load the path values from the supplied configuration file
# Use the relative path here since we haven't switched dirs
. $PATH_CONFIG_FILE

# Make sure the hpc configuration file exists
if [ ! -e ${HPC_CONFIG_FILE} ]; then
    echo "HPC configuration file ${HPC_CONFIG_FILE} not found."
    echo "Defined in ${PATH_CONFIG_FILE}."
    exit 1
fi

# DEFINE THE LAYOUT OF THE MACHINE BEING USED
. $HPC_CONFIG_FILE
MACH_LAYOUT

validate_paths_config

# Default value for the queue to submit the pmo job to should be the
# same as filter, since they both rely on waiting for the coamps
# integration to finish
: ${PBS_QUEUE_PMO:=$PBS_QUEUE_FILTER}

echo "------------------------------------------"
echo "  Setting up DART ensemble experiment"
echo "------------------------------------------"
echo "The COAMPS model directory is            : $COAMPS_HOME"
echo "The DART COAMPS directory is             : $DART_HOME"
echo "The COAMPS data for this experiment is   : $COAMPS_DATA"
echo "The DART experiment base directory is    : $DART_BASE"
echo "------------------------------------------"

# Directory Shortcuts
: ${DART_WORK:=$DART_HOME/work}
: ${DART_SCRIPTS:=$DART_HOME/shell_scripts}
: ${DART_SCRIPTS_TEMPLATE:=$DART_HOME/shell_scripts/TEMPLATES}

# Convert the configuration file name to something absolute
echo "Converting configuration to absolute path..."
ABS_CONFIG_FILE=`${DART_SCRIPTS}/rel_to_abs_path.sh ${PATH_CONFIG_FILE}`

: ${DART_SCRIPTS_WORK:=$DART_BASE/scripts}
: ${DART_BIN:=$DART_BASE/bin}

# Build the dart binaries
: ${build_dart:='f'}
if [ $build_dart = 't' ]; then
${DART_SCRIPTS}/dart_build.sh -c ${PATH_CONFIG_FILE} 
check_return $? dart_build.sh
fi

# Create the basic DART directories
echo "Creating DART experiment directory..."
rm -rf $DART_BASE
mkdir -p $DART_BASE
mkdir -p $DART_BIN
mkdir -p $DART_LOG
mkdir -p $DART_SCRIPTS_WORK
#rsh ${ARCMACH} "mkdir -p ${DART_ARCH}"

cd $DART_BASE
${LN} ${ABS_CONFIG_FILE} ${DART_BASE}

# Copy the data files into the DART working directory
echo "Copying data files to DART working directory..."
${DART_SCRIPTS}/populate_data_files.sh -c ${ABS_CONFIG_FILE} > ${DART_LOG}/pop_data.out
check_return $? populate_data_files.sh

# Set some default values
: ${cdtg:=$cdtg_beg}
: ${icycle_loc:=$icycle}
: ${DART_RESTART_FCST:='coamps_state_fcst'}
: ${DART_RESTART_ANAL:='coamps_state_anal'}

#
# Link the DART binaries to the DART working directory
dart_exe=( create_obs_sequence create_fixed_network_seq perfect_model_obs filter
		   trans_coamps_to_dart trans_dart_to_coamps advance_time create_mean_std
		   create_increment obs_diag obs_seq_to_netcdf innov_to_obs_seq 
		   perturb_bndy perturb_init)
echo
echo "Linking DART binaries to working directory..."
for file in ${dart_exe[@]}; do 
  echo " ${file}"
#${LN} $DART_WORK/${file} ${DART_BIN}
  ${CP} $DART_WORK/${file} ${DART_BIN}
done

# Copy over the shell scripts
echo "Copying shell scripts..."
${CP} ${DART_SCRIPTS}/HPC_CONFIG.sh        ${DART_SCRIPTS_WORK}
${CP} ${DART_SCRIPTS}/strip_namelist.pl    ${DART_SCRIPTS_WORK}

flist=( run_filter.sh advance_wrapper.sh update_namelists.sh create_mean_std.sh create_increment.sh forecast_wrapper.sh restart_wrapper.sh 
		scale_coamps_perts.sh archive_coamps_ens.sh )
if [ ${IS_REAL_DATA} == 'T' ]; then
  flist=( ${flist[@]} navdas_preproc_obs.sh ) 
else
  flist=( ${flist[@]} run_pmo.sh advance_perfect.sh update_perfect.sh )
fi

# Copy all the run scirpts to the working director
echo "Customizing scripts..."
for file in ${flist[@]}; do 
  ${CP} ${DART_SCRIPTS_TEMPLATE}/${file} ${DART_SCRIPTS_WORK} 
  customize_path_config ${file}
done

echo "  advance_wrapper.sh    - ic_file_name"
  ${SED} "s/^\(\s*ic_file\s*=\s*.printf\).*\(\.\%04.*\)/\1 ${DART_RESTART_ANAL}\2/" ${DART_SCRIPTS_WORK}/advance_wrapper.sh
echo "  advance_wrapper.sh    - ud_file_name"
  ${SED} "s/^\(\s*ud_file\s*=\s*.printf\).*\(\.\%04.*\)/\1 ${DART_RESTART_FCST}\2/" ${DART_SCRIPTS_WORK}/advance_wrapper.sh

# Also do the namelists
echo "Customizing paths in namelists..."
echo "  COAMPS namelist       - dsnrff"
${PERL} "s|(\s*dsnrff\s*=\s*\').*(\'.*$)|\1${DART_BASE}/ENSDIR/data/\2|" ./namelist
echo "  COAMPS namelist       - dsngff"
${PERL} "s|(\s*dsngff\s*=\s*\').*(\'.*$)|\1\/NOGAPS_ENSEMBLE\2|" ./namelist

echo "  input.nml             - ens_size"
  ${PERL} "s/(ens_size\s*)=\s*\d+/\1= ${ENSEMBLE_SIZE}/" input.nml 
#echo "  input.nml             - num_output_state_members"
#  ${PERL} "s/(num_output_state_members\s*)=\s*\d+/\1= ${ENSEMBLE_SIZE}/" input.nml 
echo "  input.nml             - num_output_obs_members"
  ${PERL} "s/(num_output_obs_members\s*)=\s*\d+/\1= ${ENSEMBLE_SIZE}/" input.nml 
echo "  input.nml             - start_from_restart (only filter_nml)"
  ${SED} "/filter_nml/,/smoother_nml/ s/^\(\s*start_from_restart\s*\)=.*/\1= .true./" input.nml 
echo "  input.nml             - adv_ens_command"
  ${SED} "s|^\(\s*adv_ens_command\s*=\s*.\).*\.\/\(.*\)|\1${DART_SCRIPTS_WORK}\/\2 | " input.nml
echo "  input.nml             - single_restart_file_in"
  ${PERL} "s/^(\s*single_restart_file_in\s*)=.*/\1= .false.,/" input.nml
echo "  input.nml             - single_restart_file_out"
  ${PERL} "s/^(\s*single_restart_file_out\s*)=.*/\1= .false.,/" input.nml
echo "  input.nml             - cdtg"
  ${PERL} "s/^(\s*cdtg\s*=\s*).*/\1\'${cdtg}\',/" input.nml
echo "  input.nml             - restart_in_file_name"
  ${SED} "/filter_nml/,/restart_in_file_name/ s/^\(\s*restart_in_file_name\s*\)=.*/\1= \"${DART_RESTART_FCST}\",/" input.nml
echo "  input.nml             - restart_out_file_name"
  ${SED} "/filter_nml/,/restart_out_file_name/ s/^\(\s*restart_out_file_name\s*\)=.*/\1= \"${DART_RESTART_ANAL}\",/" input.nml

# Get the number of cpus in the namelist.  This is not currently used.
XCPUS=`awk -F= '/ndxnam/{print $2}' namelist | awk -F, '{print $1}'`
YCPUS=`awk -F= '/ndynam/{print $2}' namelist | awk -F, '{print $1}'`
let "COMP_CPUS = XCPUS * YCPUS"

: ${submit_jobs:='t'}

if [ ${submit_jobs} == 't' ]; then
  cd ${DART_BASE}
  if [ ${IS_REAL_DATA} == 'T' ]; then
    ${DART_SCRIPTS_WORK}/run_filter.sh
  else
    RUN_PMO_ID=`${DART_SCRIPTS_WORK}/run_pmo.sh`
    ${DART_SCRIPTS_WORK}/run_filter.sh -W ${RUN_PMO_ID} 
  fi
fi

# Now we should have done everything we need in order to properly run everything!!!! 
echo "Finished"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

