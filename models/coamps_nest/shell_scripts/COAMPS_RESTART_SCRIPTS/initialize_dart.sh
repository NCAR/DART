#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################
#
# SCRIPT:   initialize_dart.sh
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Initialize the DART system to run experiments with the Navy's COAMPS
# model.  This takes the following actions: 
#
# 1. Create the DART experiment folder.  This is DART_WORK_DIR which
#    yields the following directory structure
#       DART_WORK_BASE
#          +-- 001
#          +-- 002
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
    validate_var DART_BASE
    validate_var DART_NAMELIST
    validate_var COAMPS_RESTART_FILE
    validate_var RESTART_DAT
    validate_var ENSEMBLE_SIZE
    validate_var CONCURRENT_RUNS
    validate_var PBS_QUEUE_ENS
    validate_var PBS_QUEUE_FILTER
    validate_var FILTER_PROCS
}


# Helper function - see if a variable is defined and abort if it
# is not.  Use indirect referencing to take in the variable name
# as a string and access its value
function validate_var
{
    var_name=$1
    if [[ -z ${!var_name} ]]
    then
        echo "Variable ${var_name} is undefined - aborting script!"
        exit 1
    fi
}

# Command Aliases
CP='cp -f'
PERL='perl -i -p -e'
SED='sed -i -e'

######
# END DEFINITIONS
######

usage="Usage: `basename $0` -c PATHCONFIG"

if [ $# -eq 0 ]
then
    echo $usage
    exit 1
fi

# Parse the options and grab the file name
while getopts ":c:" option
do
  case $option in
      c  ) PATH_CONFIG_FILE=$OPTARG;;
      \? ) echo $usage
           exit 65;;
      *  ) echo $usage
       exit 65;;
  esac
done

# Make sure the configuration file exists
if [ ! -e ${PATH_CONFIG_FILE} ]; then
    echo "Configuration file ${PATH_CONFIG_FILE} not found!"
    exit 1
fi

# Load the path values from the supplied configuration file
# Use the relative path here since we haven't switched dirs
. $PATH_CONFIG_FILE
validate_paths_config

## Handle non-specified values
# Give default values for a few of the parameters in the configuration
# file.  These have to do with a pre-specified ensemble.
: ${PERFECT_OBS_DATA:=$COAMPS_DATA}
: ${PERFECT_RESTART_FILE:=$COAMPS_RESTART_FILE}

# Default value for the queue to submit the pmo job to should be the
# same as filter, since they both rely on waiting for the coamps
# integration to finish
: ${PBS_QUEUE_PMO:=$PBS_QUEUE_FILTER}

# If we're dealing with a pre-existing ensemble, then we're going to have to use the printf
# syntax for COAMPS data.  No need to specify explicitly if the ensemble pre-exists or not:
# just check for a format character, assuming that % is *NOT* part of a normal path name - 
# change this if you have that.
if [[ `expr index ${COAMPS_DATA} '%'` == 0 ]]
then
    ENSEMBLE_PREEXISTS="no"
else
    ENSEMBLE_PREEXISTS="yes"
fi

# Add the calculated value to the paths.config file
if [[ -z `grep ENSEMBLE_PREEXISTS $PATH_CONFIG_FILE` ]]
then
    # Not found in file - add it
    echo "ENSEMBLE_PREEXISTS=${ENSEMBLE_PREEXISTS}" >> $PATH_CONFIG_FILE
else
    ${PERL} "s/(ENSEMBLE_PREEXISTS=).*/\1${ENSEMBLE_PREEXISTS}/g" $PATH_CONFIG_FILE
fi

echo "------------------------------------------"
echo "  Setting up DART ensemble experiment"
echo "------------------------------------------"
echo "The COAMPS model directory is            : $COAMPS_HOME"
echo "The DART COAMPS directory is             : $DART_HOME"
echo "The COAMPS data for this experiment is   : $COAMPS_DATA"
echo "The DART experiment base directory is    : $DART_BASE"
echo "------------------------------------------"
if [ "$ENSEMBLE_PREEXISTS" != "no" ]; then
    echo "       Using a pre-existing ensemble      "
    echo "The COAMPS data for perfect_model_obs is : $PERFECT_OBS_DATA"
    echo "------------------------------------------"
fi

# Directory Shortcuts
DART_WORK=$DART_HOME/work
DART_SCRIPTS=$DART_HOME/shell_scripts

# Convert the configuration file name to something absolute
echo "Converting configuration to absolute path..."
ABS_CONFIG_FILE=`${DART_SCRIPTS}/rel_to_abs_path.sh ${PATH_CONFIG_FILE}`

# Create the basic DART directories
echo "Creating DART experiment directory..."
rm -rf $DART_BASE
mkdir -p $DART_BASE

# Copy the data files into the DART working directory
echo "Copying data files to DART working directory..."
cd $DART_BASE
${DART_SCRIPTS}/populate_data_files.sh -c ${ABS_CONFIG_FILE}
check_return $? populate_data_files.sh

# Define the restart.vars file based on our data file - this creates
# a "restart.vars" file in the current directory.
$DART_SCRIPTS/populate_restart_vars.pl ${RESTART_DAT} > pop_restart.out
check_return $? populate_restart_vars.pl

# Copy over the shell scripts
echo "Copying shell scripts..."
${CP} ${DART_SCRIPTS}/advance_*.csh .
${CP} ${DART_SCRIPTS}/run_filter.csh .
${CP} ${DART_SCRIPTS}/run_pmo.csh .
${CP} ${DART_SCRIPTS}/initialize_ensemble.sh .
${CP} ${DART_SCRIPTS}/job_setup.csh .

# Change the paths - we only need to do this for advance_model.csh since
# advance_group.csh only uses PBS_O_WORKDIR  
# Do this with perl instead of sed since I can't figure out how to
# make it handle the path name variable interpolation properly
echo "Updating paths in shell scripts..."
echo "  advance_model.csh - COAMPS_BIN_DIR"
${PERL} "s|^(.*\bCOAMPS_BIN_DIR\s*=\s*).*$|\1${COAMPS_HOME}/bin|" advance_model.csh
echo "  advance_model.csh - COAMPS_DART_DIR"
${PERL} "s|^(.*\bCOAMPS_DART_DIR\s*=\s*).*$|\1${DART_HOME}|" advance_model.csh
echo "  advance_model.csh - COAMPS_ENS"
${PERL} "s|^(.*\bCOAMPS_ENS\s*=\s*).*({.*})\"?$|\1${DART_BASE}/\\$\2|" advance_model.csh

# Also do the namelists
echo "Updating paths in namelists..."
echo "  COAMPS namelist - dsnrff"
${PERL} "s|(\s*dsnrff\s*=\s*\').*(\'.*$)|\1${DART_BASE}/ENSDIR/data/\2|" namelist
echo "  initialize_ensemble.sh - COAMPS_DIR"
${PERL} "s|^(\s*COAMPS_DIR=).*$|\1${COAMPS_DATA}|" initialize_ensemble.sh

# Update the ensemble size
echo "Updating ensemble size count to ${ENSEMBLE_SIZE}..."
echo "  initialize_ensemble.sh - NUM_MEMS"
${PERL} "s/(NUM_MEMS)=\d+/\1=${ENSEMBLE_SIZE}/" initialize_ensemble.sh
echo "  run_filter.csh - num_ens"
${PERL} "s/(num_ens)\s*=\s*\d+/\1 = ${ENSEMBLE_SIZE}/" run_filter.csh
echo "  input.nml - ens_size"
${PERL} "s/(ens_size\s*)=\s*\d+/\1= ${ENSEMBLE_SIZE}/" input.nml 
echo "  input.nml - num_output_state_members"
${PERL} "s/(num_output_state_members\s*)=\s*\d+/\1= ${ENSEMBLE_SIZE}/" input.nml 
echo "  input.nml - num_output_obs_members"
${PERL} "s/(num_output_obs_members\s*)=\s*\d+/\1= ${ENSEMBLE_SIZE}/" input.nml 

# Need to be careful here - there are several namelists in input.nml that have
# the start_from_restart parameter
if [ "${ENSEMBLE_PREEXISTS}" != "no" ]
then
    echo "Updating namelists to reflect pre-existing ensemble..."
    echo "  input.nml - start_from_restart (only filter_nml)"
    ${SED} "/filter_nml/,/smoother_nml/ s/^\(\s*start_from_restart\s*\)=.*/\1= .true./" input.nml 
    echo "  input.nml - single_restart_file_in"
    ${PERL} "s/^(\s*single_restart_file_in\s*)=.*/\1= .false.,/" input.nml
else
    echo "Updating namelists to reflect no pre-existing ensemble..."
    echo "  input.nml - start_from_restart (only filter_nml)"
    ${SED} "/filter_nml/,/smoother_nml/ s/^\(\s*start_from_restart\s*\)=.*/\1= .false./" input.nml 
    echo "  input.nml - single_restart_file_in"
    ${PERL} "s/^(\s*single_restart_file_in\s*)=.*/\1= .true./" input.nml
fi

# Update the DART DTG
echo "Updating date-time group"
echo "  input.nml - cdtg"
CDTG=`perl -n -e 'if (m/cdtg\s*=\s*(.\d{10}.)/) {print $1;}' namelist`
${PERL} "s/(cdtg\s*=\s*)'\d{10}'/\1${CDTG}/" input.nml

# It's silly to try and run more members than the ensemble size
# concurrently, so reset it if they are too big.
if [ ${CONCURRENT_RUNS} -gt ${ENSEMBLE_SIZE} ]
then
    CONCURRENT_RUNS=${ENSEMBLE_SIZE}
fi

# Hopefully this should never have to be done
if [ ${CONCURRENT_RUNS} -le 0 ]
then
    CONCURRENT_RUNS=1
fi

# Update the number of CPUs that we need to use - this figures in to
# advance_group.csh and how to handle the queue system.
# Do this by taking the total number of CPUs required by the COAMPS
# namelist for COAMPS to run, multiplying it by the number of
# concurrent runs we'll handle,  then figure out how many 2-CPU nodes
# we need and update the processor count and the queue request.

# Count the X and Y CPUS - need dual awk programs to split up the
# namelist by whitespace, then by commas - assume that we only pay
# attention to the first entry in the list of 1,1,1,1,1,1,1,1,
XCPUS=`awk -F= '/ndxnam/{print $2}' namelist | awk -F, '{print $1}'`
YCPUS=`awk -F= '/ndynam/{print $2}' namelist | awk -F, '{print $1}'`
let "COMP_CPUS = XCPUS * YCPUS"
IOCPUS=`awk -F= '/npr0nam/{print $2}' namelist | awk -F, '{print $1}'`
let "NCPUS = COMP_CPUS + IOCPUS"
let "TOTALCPUS = NCPUS * CONCURRENT_RUNS"
let "NODES = (TOTALCPUS / 2) + (TOTALCPUS % 2)"
let "PMO_NODES = (NCPUS / 2) + (NCPUS % 2)"

# For perfect model obs
echo "Updating perfect model PBS request to use ${PMO_NODES} node(s)..."
echo "  run_pmo.csh - PBS -lnodes"
${PERL} "s/(\#PBS.*nodes)=\d+/\1=${PMO_NODES}/" run_pmo.csh

# For the ensemble integration
echo "Updating ensemble PBS request to use $NODES node(s)..."
echo "  advance_group.csh - PBS -lnodes"
${PERL} "s/(\#PBS.*nodes)=\d+/\1=${NODES}/" advance_group.csh
echo "  advance_group.csh - NPROCS"
${PERL} "s/(.*NPROCS)\s*=\s*\d+/\1 = ${NCPUS}/" advance_group.csh

# For perfect model obs
let "PMO_NODES = (NCPUS / 2) + (NCPUS % 2)"
echo "Updating perfect model PBS request to use ${PMO_NODES} node(s)..."
echo "  run_pmo.csh - PBS -lnodes"
${PERL} "s/(\#PBS.*nodes)=\d+/\1=${PMO_NODES}/" run_pmo.csh

# For the MPI version of filter
let "FILTER_NODES = (FILTER_PROCS / 2) + (FILTER_PROCS % 2)"
echo "Updating filter PBS request to use $FILTER_NODES node(s)..."
echo "  run_filter.csh - PBS -lnodes"
${PERL} "s/(\#PBS.*nodes)=\d+/\1=${FILTER_NODES}/" run_filter.csh
echo "  run_filter.csh - filterprocs"
${PERL} "s/(.*filterprocs)\s*=\s*\d+/\1 = ${FILTER_PROCS}/" run_filter.csh

# Update the queue information
echo "Updating ensemble PBS request to use queue ${PBS_QUEUE_ENS}..."
echo "  advance_group.csh"
${PERL} "s/(\#PBS\s*-q).*/\1 ${PBS_QUEUE_ENS}/" advance_group.csh
echo "Updating perfect model PBS request to use queue ${PBS_QUEUE_PMO}..."
echo "  run_pmo.csh"
${PERL} "s/(\#PBS\s*-q).*/\1 ${PBS_QUEUE_PMO}/" run_pmo.csh
echo "Updating filter PBS request to use queue ${PBS_QUEUE_FILTER}..."
echo "  run_filter.csh"
${PERL} "s/(\#PBS\s*-q).*/\1 ${PBS_QUEUE_FILTER}/" run_filter.csh

# Update the number of concurrent runs
echo "Making $CONCURRENT_RUNS concurrent ensemble runs available..."
echo "  advance_group.csh - MAX_RUNNING"
${PERL} "s/(.*MAX_RUNNING)\s*=\s*\d+/\1 = ${CONCURRENT_RUNS}/" advance_group.csh
echo "  advance_model.csh - pernode parameter"
if [ ${CONCURRENT_RUNS} -lt 2 ]
then
  ${PERL} "s/-pernode//g" advance_model.csh
fi

# Get the DART externals into their proper location
echo "Installing COAMPS external files to the DART tree..."
${DART_SCRIPTS}/install_externals.sh ${DART_HOME}/../../ > ${DART_BASE}/external.log
check_return $? install_externals.sh

# Generate the get_name_info file based on the COAMPS code
echo "Extracting restart file parameters from COAMPS..."
cd ${DART_HOME}
${DART_SCRIPTS}/generate_restart_field_list.sh "${COAMPS_HOME}" ${DART_SCRIPTS}
check_return $? generate_restart_field_list.sh
${DART_SCRIPTS}/generate_get_name_function.pl fields_master_list
check_return $? generate_get_name_function.pl
rm fields_master_list

# Create the module containing the COAMPS utility subroutines
echo "Creating module from COAMPS utility package..."
(cd ${DART_HOME}; ${DART_SCRIPTS}/create_coamps_intrinsic_mod.sh ${COAMPS_UTIL_HOME})

# Build the DART binaries
echo "Building DART binaries..."
cd ${DART_WORK}
cp $DART_BASE/input.nml .
${DART_SCRIPTS}/quickbuild.sh > ${DART_BASE}/make.out
check_return $? quickbuild.sh

# Link the DART binaries to the DART working directory
echo
echo "Linking DART binaries to working directory..."
cd $DART_BASE
echo " create_obs_sequence"
ln -sf $DART_WORK/create_obs_sequence .
echo " create_fixed_network_seq"
ln -sf $DART_WORK/create_fixed_network_seq .
echo " perfect_model_obs"
ln -sf $DART_WORK/perfect_model_obs .
echo " filter"
ln -sf $DART_WORK/filter .
echo " wakeup_filter"
ln -sf $DART_WORK/wakeup_filter .
echo " trans_coamps_to_dart"
ln -sf $DART_WORK/trans_coamps_to_dart .
echo " trans_dart_to_coamps"
ln -sf $DART_WORK/trans_dart_to_coamps .

# Create DART initial condition file
echo "Creating DART initial condition file..."
${DART_SCRIPTS}/convert_coamps_to_dart.sh -c ${ABS_CONFIG_FILE} \
  -o perfect_ics -i "${PERFECT_RESTART_FILE}" > convert.perfect_ics.log 2>&1
check_return $? convert_coamps_to_dart.sh

# Create DART ensemble initial conditions file
if [ $ENSEMBLE_PREEXISTS != "no" ]; then
    echo "Creating DART ensemble initial conditions files..."
    for member in `seq 1 ${ENSEMBLE_SIZE}`
    do
      # Multiple files have a four-digit member extension
      FILTER_IC_NAME=`printf "filter_ics.%04d" ${member}`
      echo " ${FILTER_IC_NAME}"
      
      RESTART_INPUT=`printf ${COAMPS_RESTART_FILE} ${member}`
      ${DART_SCRIPTS}/convert_coamps_to_dart.sh -c ${ABS_CONFIG_FILE} \
      -o ${FILTER_IC_NAME} -i "${RESTART_INPUT}" > convert.${FILTER_IC_NAME}.log 2>&1
      check_return $? convert_coamps_to_dart.sh
    done
else
    echo "Creating DART ensemble initial conditions file..."
    cp perfect_ics filter_ics
fi

# Create and populate the directories for both the individual
# ensemble members (named with the number of the ensemble member)
# and the perfect model/truth member (named with "perfect")
echo "Creating perfect model and ensemble member directories..."
RESTART=`basename "${PERFECT_RESTART_FILE}"`
echo "  Only copying ${RESTART}"
${PERL} "s/(rsync.*)restart\*/\1${RESTART}/" initialize_ensemble.sh


# Copy over all the run scripts and tweak them so they handle the
# special case of a perfect model run, then use the newly modified
# initialization code to copy over the perfect model data
${DART_SCRIPTS}/create_pmo_files.sh $(( NCPUS/2 + NCPUS%2 )) \
                                    ${PERFECT_OBS_DATA}

echo "  Initializing perfect model directory"
./initialize_perfect_model.sh > perfect_initialization.log 

# Initialize the ensemble - modify script so it only points to the
# restart file that we use as an initial condition.
echo "  Initializing ensemble directories"
./initialize_ensemble.sh > ensemble_initialization.log

# Now we should have done everything we need in order to properly run
# everything!!!! 
rm *.bak
echo "Finished!"

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

