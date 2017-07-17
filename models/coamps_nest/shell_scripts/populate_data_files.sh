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
# Populate the data files necessary to run the DART
# programs.  These include:
#   - input.nml    (DART)
#   - namelist     (COAMPS)
#   - datahd*      (COAMPS domain info)
#   - terrht*      (COAMPS terrain heights)
#   - state.vars (DART state vector definition)
#   - convert.vars (DART <> COAMPS converter variables)
#
# This script requires the -c option of the path configuration file.
######

usage="Usage: `basename $0` -c PATHCONFIG [-m MEMBER]"

if [ $# -eq 0 ]
then
    echo $usage
    exit 1
fi

# Parse the options and grab the file name
while getopts ":c:m:" option
do
  case $option in
      c  ) PATH_CONFIG_FILE=$OPTARG;;
      m  ) MEMBER=$OPTARG;;
      \? ) echo $usage
           exit 65;;
      *  ) echo $usage
       exit 65;;
  esac
done
: ${MEMBER:=-1}

# Make sure the configuration file exists
if [ ! -e "${PATH_CONFIG_FILE}" ]; then
    echo "Configuration file ${PATH_CONFIG_FILE} not found!"
    exit 1
fi

# Grab the definitions
. "${PATH_CONFIG_FILE}"
: ${is_fcp_bndy:='f'}

# For ease of use, just assume that we will be copying things
# in the current working directory
echo "Populating DART & COAMPS data files in `pwd`"

# Bring in the files that we need
# Get the desired restart file (supplied in the path configuration)
echo "Copying in data files..."
echo "  DART namelist"
cp -f "${DART_NAMELIST}" ./input.nml
echo "  COAMPS namelist"
cp -f "${COAMPS_NAMELIST}" ./namelist
echo "  Configure file"
ln -sf "${PATH_CONFIG_FILE}" ./

if [ -e ${COS_FILE} ]; then
  echo "  copying $COS_FILE"
  cp -f "${COS_FILE}" ./coseq.dat
fi

#if [ -e ${CFN_FILE} ]; then
#  echo "  copying $CFN_FILE"
#  cp -f "${CFN_FILE}" ./cfn.dat
#fi

# If we don't have a pre-existing ensemble, the data directory is
# just a single directory
echo "  COAMPS data files"
#if [ "${ENSEMBLE_PREEXISTS}" != "no" ]; then
#    DATA_PATH="${PERFECT_OBS_DATA}"
#else
#    DATA_PATH="${COAMPS_DATA}"
#fi

if [ $MEMBER = 'perfect' ]; then
    DATA_PATH=${PERFECT_DATA} ; MEMBER=1
else
    DATA_PATH=`printf ${COAMPS_DATA} $MEMBER`
fi

if [ $MEMBER -ge 0 ]; then
  cp -f ${DATA_PATH}/datahd* .
  cp -f ${DATA_PATH}/terrht* .
fi

echo ${DART_HOME}/convert.vars
cat ${DART_HOME}/convert.vars

cp -f ${DART_HOME}/convert.vars ./convert.vars
echo "  State vector definition"
${DART_HOME}/shell_scripts/populate_state_vars.pl ${STATE_DEF_FILE} ${FLAT_FILES} > ${DART_LOG}/state_var.log

if [ $is_fcp_bndy == 'true' ]; then
  cp -f ${DART_HOME}/perturb.vars ./perturb.vars
fi

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

