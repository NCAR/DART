#!/bin/csh

# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$

#--------------------------
# Script to advance one ensemble member one filter "time step"
# when the model advance is executed as a separate process.
# Called by filter_server.csh.
# Calls run-pc.csh, the CAM execution script.
# Calls 3 translation routines to translate time and model state.
# Runs on one of the compute nodes allotted to filter_server.csh.
#--------------------------

set CENTRALDIR = $1
set element = $2
set temp_dir = $3

if ($#argv == 4) set machine_file = $4
# 'machine_file' IS NOT USED by run-pc.csh by single-threaded executions; see alternate call below

echo "starting advance_model.csh for ens member $element at "`date` > cam_out_temp

if ( -d $temp_dir ) then
   cd   $temp_dir
   \rm -f *
   echo "advance_model; cd to $temp_dir" >> cam_out_temp
else
   echo "FATAL ERROR advance_model.csh ... temp_dir( ${temp_dir} ) does not exist."
   echo "FATAL ERROR advance_model.csh ... temp_dir( ${temp_dir} ) does not exist." >> cam_out_temp
   exit 99
endif

# Get information about this experiment from file "casemodel", created by the
# main controlling script (job.csh)
set caseinfo = `cat $CENTRALDIR/casemodel`
# set $case = the case and 
# $model = the directory name (in the central CAM directory) where CAM executable will be found. 
set case = $caseinfo[1]
set model = $caseinfo[2]
# set locations of the CAM and CLM input files
set cam_init = $caseinfo[3]
set clm_init = $caseinfo[4]

# output diagnostic information to the same file as the CAM list-directed output
echo "case $case model $model"    >> cam_out_temp
echo "cam init is $cam_init"      >> cam_out_temp
echo "clm init is $clm_init"      >> cam_out_temp

# get model state initial conditions for this ensemble member
ln -s $CENTRALDIR/assim_model_state_ic$element temp_ic
# get filter namelists for use by 
cp $CENTRALDIR/input.nml input.nml

echo "ls $temp_dir for element $element" >> cam_out_temp
echo junk > element$element
ls -lRt >> cam_out_temp
 
# Need a base CAM initial file into which to copy state vector from filter.
# c[al]minput_$element also carry along CAM/CLM fields which are not updated
#      by the filter (not part of the filter model state).
# First look for c[al]minput.nc resulting from the previous advance of this ensemble
#      member from within the same day/obs_seq.out time span (in CENTRALDIR)
# Failing that, look for the results of the last advance of this ensemble member
#      of the previous obs_seq.out (i.e. in CENTRALDIR/exp_name/day/CAM)
# Failing that (when starting an experiment which has no spun up set of members)
#      get a copy of a single CAM initial file (usually from somewhere independent
#      of this experiment, i.e. /scratch/.../New_state/T42_GWD/CAM/caminput_0.nc)

if (-e     ${CENTRALDIR}/caminput_$element.nc) then
   cp -p   ${CENTRALDIR}/caminput_$element.nc caminput.nc
   echo "CENTRALDIR caminput comes from ${CENTRALDIR}/caminput_num" >> cam_out_temp
else if (-e ${cam_init}$element.nc) then
   cp -p    ${cam_init}$element.nc caminput.nc
   echo "cam_init caminput comes from ${cam_init}num" >> cam_out_temp
else
   cp -p   ${cam_init}0.nc caminput.nc
   echo "DEFAULT caminput comes from ${cam_init}0.nc" >> cam_out_temp
endif

if (-e    ${CENTRALDIR}/clminput_$element.nc) then
    cp -p ${CENTRALDIR}/clminput_$element.nc clminput.nc
else if (-e ${clm_init}$element.nc) then
    cp -p   ${clm_init}$element.nc clminput.nc
else
    cp -p   ${clm_init}0.nc clminput.nc
endif

# create 'times' file for CAM from DART times in assim_model_state_ic#
# This info is passed to CAM through the creation of its namelist
if (-e temp_ic && -e $CENTRALDIR/trans_time) then
   echo 'advance_model; executing trans_time '`date` >> cam_out_temp
   $CENTRALDIR/trans_time
   ls -lt >> cam_out_temp
   cp -p times $CENTRALDIR
else
   echo 'either no ic file or trans_time available for trans_time'
   exit 1
endif

# Create an initial CAM.nc file from the DART state vector
# Times are handled separately in trans_time
$CENTRALDIR/trans_sv_pv
ls -ltR >> cam_out_temp

# advance cam 
if ($#argv == 4) then
    ${model:h}/run-pc.csh ${case}-$element $model $CENTRALDIR $machine_file >>& cam_out_temp
else
    ${model:h}/run-pc.csh ${case}-$element $model $CENTRALDIR >>& cam_out_temp
endif

grep 'END OF MODEL RUN' cam_out_temp > /dev/null
if ($status == 0) then
   # Extract the new state vector information from the new caminput.nc and
   # put it in temp_ud (time followed by state)
   $CENTRALDIR/trans_pv_sv

   # Move updated state vector and new CAM/CLM initial files back to experiment
   # directory for use by filter and the next advance.
   mv temp_ud     $CENTRALDIR/assim_model_state_ud$element
   mv clminput.nc $CENTRALDIR/clminput_$element.nc
   mv caminput.nc $CENTRALDIR/caminput_$element.nc
   mv namelist    $CENTRALDIR

   echo "finished advance_model.csh for ens member $element at "`date` >> cam_out_temp
   mv cam_out_temp $CENTRALDIR/cam_out_temp$element
   cd $CENTRALDIR
# orig;   rm -rf $temp_dir
# debug lightning
   rm $temp_dir/*
else
   echo "CAM stopped abnormally; filter_server will abort, see $temp_dir too" >> cam_out_temp
# orig;   mv cam_out_temp $CENTRALDIR/cam_out_temp$element
# debug lightning
   cp cam_out_temp $CENTRALDIR/cam_out_temp$element
   mkdir ../dead_$element
   mv * ../dead_$element
   mv ../dead_$element/bld .
# debug end 
endif
