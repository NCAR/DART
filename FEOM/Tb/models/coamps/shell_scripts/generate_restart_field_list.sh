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
# Script to automatically generate the list of fields in a COAMPS restart file
# for use with the dynamic DART state vector definition, to avoid having to 
# retype everything or adjust quickly if the routine changes.
# 
# Command-line arguments are the path to the COAMPS data and the path
# to the parser script
######
 
validate_input()
{
    if [[ ! -e ${1} ]]; then
        echo "Could not find ${1}!"
        exit 1
    fi
}
COAMPS_DIR=$1
PARSE_DIR=$2

multi_io_restart_routine=${COAMPS_DIR}/atmos/libsrc/amlib/write_restart.F
single_io_restart_routine=${COAMPS_DIR}/atmos/libsrc/amlib/write_restart_io.F

PARSER=${PARSE_DIR}/parse_field_list.pl

validate_input ${multi_io_restart_routine}
validate_input ${single_io_restart_routine}
validate_input ${PARSER}

grep -A1 MPI_FILE_WRITE ${multi_io_restart_routine} | grep "adom(inest)" | grep '^[^c].*' > field_list_multi_io
grep -A1 rest_ ${single_io_restart_routine} | grep "adom(nn)" | grep '^[^c].*' > field_list_single_io

# Multiple processor
${PARSER} 'anca' 'cond' 'aalhs' 'w0avg' MULTIIO field_list_multi_io > fields_multi_io
${PARSER} 'anca' 'cond' 'aalhs' 'wbs' SINGLEIO field_list_single_io > fields_single_io

cat fields_*_io > fields_master_list
rm field_list_*_io
rm fields_*_io

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

