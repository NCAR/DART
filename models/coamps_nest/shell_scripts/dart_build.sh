#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################
#
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
    validate_var DART_BASE
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

######
# END DEFINITIONS
######

usage="Usage: `basename $0` -c PATHCONFIG"

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

DART_SCRIPTS=$DART_HOME/shell_scripts
DART_WORK=$DART_HOME/work

build='t'
if [ $build = 't' ]; then
# Get the DART externals into their proper location
echo "Installing COAMPS external files to the DART tree..."
(
${DART_SCRIPTS}/install_externals.sh ${DART_HOME}/../../ > ${DART_WORK}/external.log
check_return $? install_externals.sh
)

# Generate the get_name_info file based on the COAMPS code
echo "Extracting restart file parameters from COAMPS..."
(
cd ${DART_HOME}
${DART_SCRIPTS}/generate_restart_field_list.sh "${COAMPS_HOME}" ${DART_SCRIPTS}
check_return $? generate_restart_field_list.sh
${DART_SCRIPTS}/generate_get_name_function.pl fields_master_list
check_return $? generate_get_name_function.pl
rm -rf fields_master_list
)

# Create the module containing the COAMPS utility subroutines
echo "Creating module from COAMPS utility package..."
(cd ${DART_HOME}; ${DART_SCRIPTS}/create_coamps_intrinsic_mod.sh ${COAMPS_UTIL_HOME})
fi

# Build the DART binaries
echo "Building DART binaries..."
cd ${DART_WORK}
echo "cp $DART_BASE/input.nml ."
cp -f "${DART_NAMELIST}" ./input.nml
echo "${DART_SCRIPTS}/quickbuild.sh > ${DART_WORK}/make.out"
${DART_SCRIPTS}/quickbuild.sh | tee ${DART_WORK}/make.out
check_return $? quickbuild.sh

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

