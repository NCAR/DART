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
# Runs on one of the compute nodes alloted to filter_server.csh.
# 
# latest revision;  Kevin Raeder 11/11/04
#--------------------------

set PBS_O_WORKDIR = $1
set element = $2
set temp_dir = $3
if ($#argv == 4) set machine_file = $4
# THIS IS NOT USED by run-pc.csh on anchorage; see altered call below

# Get information about this experiment from file "casemodel", created by the
# main controlling script (job.csh)
set caseinfo = `cat $PBS_O_WORKDIR/casemodel`
# set $case = the case and 
# $model = the directory name (in the Central directory) where CAM will be found. 
set case = $caseinfo[1]
set model = $caseinfo[2]
# set locations of the CAM and CLM input files
set cam_init = $caseinfo[3]
set clm_init = $caseinfo[4]

# output diagnostic information to the same file as the CAM list-directed output
echo $case $model > cam_out_temp
echo $cam_init  >> cam_out_temp
echo $clm_init >> cam_out_temp

mkdir $temp_dir
cd $temp_dir

# get model state initial conditions for this ensemble member
cp $PBS_O_WORKDIR/assim_model_state_ic$element temp_ic
# get filter namelists for use by 
cp $PBS_O_WORKDIR/input.nml input.nml

echo ls $temp_dir for element $element
echo junk > element$element
ls -lRt >> cam_out_temp
 
# Need a base CAM initial file into which to copy state vector from filter.
# c[al]minput_$element also carry along CAM/CLM fields which are not updated
#      by the filter (not part of the filter model state).
# First look for c[al]minput.nc resulting from the previous advance of this ensemble
#      member from within the same day/obs_seq.out time span (in PBS_O_WORKDIR)
# Failing that, look for the results of the last advance of this ensemble member
#      of the previous obs_seq.out (i.e. in PBS_O_WORKDIR/exp_name/day/CAM)
# Failing that (when starting an experiment which has no spun up set of members)
#      get a copy of a single CAM initial file (usually from somewhere independent
#      of this experiment, i.e. /scratch/.../New_state/T42_GWD/CAM/caminput_0.nc)

if (-e ${PBS_O_WORKDIR}/caminput_$element.nc) then
    cp ${PBS_O_WORKDIR}/caminput_$element.nc caminput.nc
    echo caminput comes from ${PBS_O_WORKDIR}/caminput_num >> cam_out_temp
else if (-e ${cam_init}$element.nc) then
         cp ${cam_init}$element.nc caminput.nc
         echo caminput comes from ${cam_init}num >> cam_out_temp
else
   cp ${cam_init}0.nc caminput.nc
   echo caminput comes from ${cam_init}0.nc >> cam_out_temp
endif

if (-e ${PBS_O_WORKDIR}/clminput_$element.nc) then
    cp ${PBS_O_WORKDIR}/clminput_$element.nc clminput.nc
else if (-e ${clm_init}$element.nc) then
         cp ${clm_init}$element.nc clminput.nc
else
   cp ${clm_init}0.nc clminput.nc
endif

# create 'times' file for CAM from DART times in assim_model_state_ic#
# This info is passed to CAM through the creation of its namelist
if (-e temp_ic && -e $PBS_O_WORKDIR/trans_time) then
   echo 'advance_model; executing trans_time' >> cam_out_temp
   $PBS_O_WORKDIR/trans_time
   ls -lt 
   cp times $PBS_O_WORKDIR
else
   echo 'either no ic file or trans_time available for trans_time'
   exit 1
endif

# Create an initial CAM.nc file from the DART state vector
# Times are handled separately in trans_time
$PBS_O_WORKDIR/trans_sv_pv
ls -ltR 

# advance cam 
if ($#argv == 4) then
   $PBS_O_WORKDIR/$model/models/atm/cam/bld/run-pc.csh $case$element \
      $model $PBS_O_WORKDIR $machine_file > cam_out_temp
else
   $PBS_O_WORKDIR/$model/models/atm/cam/bld/run-pc.csh $case$element \
      $model $PBS_O_WORKDIR > cam_out_temp
endif
mv cam_out_temp $PBS_O_WORKDIR/cam_out_temp$element

# Extract the new state vector information from the new caminput.nc and
# put it in temp_ud (time followed by state)
$PBS_O_WORKDIR/trans_pv_sv

# Move updated state vector and new CAM/CLM initial files back to experiment
# directory for use by filter and the next advance.
mv temp_ud $PBS_O_WORKDIR/assim_model_state_ud$element
mv clminput.nc $PBS_O_WORKDIR/clminput_$element.nc
mv caminput.nc $PBS_O_WORKDIR/caminput_$element.nc

mv namelist $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
grep 'END OF MODEL RUN' cam_out_temp$element > /dev/null
if ($status == 0) rm -rf $temp_dir
