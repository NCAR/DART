#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Convert a COAMPS restart file to a DART restart file.
# 
# The command line argument(s) specify the location of the paths.config
# file that's used to share filename/path information between scripts,
# an optional input file name and an output file name.  Since we're using
# the bash "getopts" routines, these can be in any order.
######

# Argument should supply where we can find the path configuration file
usage="Usage: `basename $0` -c PATHCONFIG -o OUTPUTFILE [-i INPUTFILE]"

if [ $# -eq 0 ]
then
    echo $usage
    exit 1
fi

# Parse the options and grab the file name
# For ease of processing with pre-existing ensembles, allow the
# command line specification of an input filename as well.
while getopts ":c:o:i:" option
do
  case $option in
      c  ) PATH_CONFIG_FILE=$OPTARG;;
      o  ) REL_OUTPUT_FILE=$OPTARG;;
      i  ) REL_INPUT_FILE=$OPTARG;;
      \? ) echo $usage
           exit 65;;
      *  ) echo $usage
       exit 65;;
  esac
done

# Load the path values from the supplied configuration file
. $PATH_CONFIG_FILE

# Directory Shortcuts
DART_WORK=$DART_HOME/work
DART_SCRIPTS=${DART_HOME}/shell_scripts

# If we didn't use the optional argument, the relative input file is
# simply the COAMPS restart file specified in the config file
if [ -n "$REL_INPUT_FILE" ]; then
    echo "Using command-line specified input file..."
else
    echo "Using configuration-file specified input file..."
fi
: ${REL_INPUT_FILE:=$COAMPS_RESTART_FILE}
echo "$REL_INPUT_FILE"

# Make absolute path so we can call other things
INPUT_FILE=`${DART_SCRIPTS}/rel_to_abs_path.sh "${REL_INPUT_FILE}"`
OUTPUT_FILE=`${DART_SCRIPTS}/rel_to_abs_path.sh "${REL_OUTPUT_FILE}"`
ABS_CONFIG_FILE=`${DART_SCRIPTS}/rel_to_abs_path.sh "${PATH_CONFIG_FILE}"`

# Create temporary working directory
echo "Creating temporary directory for conversion..."
CURDIR=`pwd`
CONVERT_TEMP=`mktemp -d ${CURDIR}/convert.XXXXXX`
cd ${CONVERT_TEMP}

# Count the X and Y CPUS - need dual awk programs to split up the
# namelist by whitespace, then by commas - assume that we only pay
# attention to the first entry in the list of 1,1,1,1,1,1,1,1,
XCPUS=`awk -F= '/ndxnam/{print $2}' ${COAMPS_NAMELIST} | awk -F, '{print $1}'`
YCPUS=`awk -F= '/ndynam/{print $2}' ${COAMPS_NAMELIST} | awk -F, '{print $1}'`
IOCPUS=`awk -F= '/npr0nam/{print $2}' ${COAMPS_NAMELIST} | awk -F, '{print $1}'`

# Get the desired restart file (supplied in the path configuration or
# specified on the command line) - if it's from the command line,
# we should use an absolute path
echo "Linking in restart file(s)..."
ln -sf ${INPUT_FILE} .

# Grab the rest of the data files
${DART_HOME}/shell_scripts/populate_data_files.sh -c ${ABS_CONFIG_FILE}

# ------------------------------------
# restarta1p00120060725000030000.nest1
# 123456789012345678901234567890123456
# restarta1pNNN20060725000030000
# ------------------------------------
# Pull stuff out of the filename that we need - use bash string
# operations for this to avoid `echo $something | cut -c1-2`.
RESTART=`expr match "$INPUT_FILE" '.*\(restarta1p.*\)'`
RESTART_DTG=${RESTART:13:10}
RESTART_TAUH=${RESTART:23:3}
RESTART_TAUM=${RESTART:26:2}
RESTART_TAUS=${RESTART:28:2}

cat > ./convert.nml <<EOF
&convert
 ktaust  = $RESTART_TAUH,$RESTART_TAUM,$RESTART_TAUS,
 ktauf   = $RESTART_TAUH,$RESTART_TAUM,$RESTART_TAUS,
 cdtg    = '$RESTART_DTG',
 ndxnam  = $XCPUS,
 ndynam  = $YCPUS,
 npr0nam = $IOCPUS/
EOF

# Do the conversion to a DART restart file
echo "Converting $RESTART..."
${DART_WORK}/trans_coamps_to_dart &> ${DART_BASE}/convert.out
if [ $? != 0 ]; then
    echo "Problem with the translation program!"
    exit 1
fi
echo "  Renaming..."
mv ./dart_vector ${OUTPUT_FILE}

# Clean up the symbolic links
cd ${CURDIR}
rm -rf ${CONVERT_TEMP}

# Now we should have done everything we need in order to properly run
# everything!!!! 
#rm *.bak
echo "Finished!"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

