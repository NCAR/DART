#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################
#
# SCRIPT:   populate_data_files.sh
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Populate the data files necessary to run the DART
# programs.  These include:
#   - input.nml    (DART)
#   - namelist     (COAMPS)
#   - datahd*      (COAMPS domain info)
#   - terrht*      (COAMPS terrain heights)
#   - restart.vars (DART state vector definition)
#   - convert.vars (DART <> COAMPS converter variables)
#
# This script requires the -c option of the path configuration file.
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
if [ ! -e "${PATH_CONFIG_FILE}" ]; then
    echo "Configuration file ${PATH_CONFIG_FILE} not found!"
    exit 1
fi

# Grab the definitions
. "${PATH_CONFIG_FILE}"

# For ease of use, just assume that we will be copying things
# in the current working directory
echo "Populating DART & COAMPS data files here:"
pwd

# Bring in the files that we need
# Get the desired restart file (supplied in the path configuration)
echo "Copying in data files..."
echo "  DART namelist"
cp -f "${DART_NAMELIST}" ./input.nml
echo "  COAMPS namelist"
cp -f "${COAMPS_NAMELIST}" ./namelist

# If we don't have a pre-existing ensemble, the data directory is
# just a single directory
echo "  COAMPS data files"
if [ "${ENSEMBLE_PREEXISTS}" != "no" ]; then
    DATA_PATH="${PERFECT_OBS_DATA}"
else
    DATA_PATH="${COAMPS_DATA}"
fi
cp -f ${DATA_PATH}/datahd* .
cp -f ${DATA_PATH}/terrht* .

cp -f ${DART_HOME}/convert.vars ./convert.vars
echo "  State vector definition"
${DART_HOME}/shell_scripts/populate_restart_vars.pl ${RESTART_DAT} > restart.log

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

