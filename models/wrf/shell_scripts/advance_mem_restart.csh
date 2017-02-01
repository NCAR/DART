#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# this script is expected to be called with the date to advance to,
# the ensemble member number, and the file which contains the global
# parameters for this run. 

# it converts the file into the right format, and sets up to call the
# advance_model.csh script (same one as used with 'filter') to run 
# the model to advance the data.  the resulting state vector is left
# in the file 'assim_model_state_ud.NNNN'.

set datea     = ${1}
set emember   = ${2}
set paramfile = ${3}

# read global parameters
source $paramfile

set start_time = `date +%s`

cd ${RUN_DIR}

set gdate = (`echo $datea $ASSIM_INT_HOURS -g | ${RUN_DIR}/advance_time`)

#  copy files to appropriate location
echo $start_time >! ${RUN_DIR}/start_member_${emember}
if ( -d ${RUN_DIR}/advance_temp${emember} )  ${REMOVE} ${RUN_DIR}/advance_temp${emember}
mkdir ${RUN_DIR}/advance_temp${emember}  ;  cd ${RUN_DIR}/advance_temp${emember}
set icnum = `printf "%04d" $emember`
${LINK} ${RUN_DIR}/filter_ic_new.${icnum} restart_file_input

#  run restart_file_tool to create assim_model_state_ic.* format file
set d = 1
while ( $d <= $domains )
   set dchar = `printf "%02d" $d`
   ${LINK} ${RUN_DIR}/wrfinput_d${dchar} .
   @ d++
end

cat >! script.sed << EOF
/new_advance_days/c\
new_advance_days             = ${gdate[1]},
/new_advance_secs/c\
new_advance_secs             = ${gdate[2]},
/num_domains/c\
num_domains = ${domains},
EOF
sed -f script.sed ${TEMPLATE_DIR}/input.nml.template >! input.nml
${RUN_DIR}/restart_file_tool
${MOVE} restart_file_output ${RUN_DIR}/assim_model_state_ic.${icnum}
cd $RUN_DIR  ;  ${REMOVE} ${RUN_DIR}/advance_temp${emember}

# construct the filter_control file that the advance_model.csh script
# is expecting to read.
echo $emember                      >! ${RUN_DIR}/filter_control${icnum}
echo assim_model_state_ic.${icnum} >> ${RUN_DIR}/filter_control${icnum}
echo assim_model_state_ud.${icnum} >> ${RUN_DIR}/filter_control${icnum}

#  integrate the model forward in time
${RUN_DIR}/advance_model.csh ${emember} 1 filter_control${icnum}
${REMOVE} ${RUN_DIR}/filter_control${icnum}

# print out how many seconds it took this script to run.
set end_time   = `date  +%s`
@ length_time  = $end_time - $start_time
echo "advance_mem_restart.csh duration in seconds = $length_time"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

