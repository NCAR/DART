#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$
#
# set echo verbose

set PBS_O_WORKDIR = $1
set element = $2

# set $case = the case and $model = the model we're running
set caseinfo = `cat $PBS_O_WORKDIR/casemodel`
set case = $caseinfo[1]
set model = $caseinfo[2]
echo $case $model

# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

# kdr this script is executed entirely in tempdir#, pulling files from elsewhere
#     except for the build in run-pc.csh, which exectutes in $blddir,
#     if I haven't created the executable cam in setup_advance_model.csh,
#     or init_adv....

# Ready for advance; make a copy of base cam and clm .nc files
# Should really keep a separate one for each ensemble member but not yet
# kdr; I am keeping them separate; in each tempdir#.

mkdir /scratch/local/dartcam
cd /scratch/local/dartcam
mv $PBS_O_WORKDIR/assim_model_state_ic$element temp_ic
echo ls /scratch/local/dartcam for element $element
echo junk > element$element
ls -lRt 

# Need a base nc file into which to copy modifications from filter
cp $PBS_O_WORKDIR/caminput.nc .
cp $PBS_O_WORKDIR/clminput.nc .

# Copy the initial condition file to the temp directory
# Need to strip out the current time (second line) and 
# leave advance to time (first line); Should all be 
# automated with more generalized model time handling
   head -1 temp_ic > temp2
   tail +3 temp_ic >> temp2
   mv temp2 temp_ic

#  is some of this caminput/CAM_FILE copying in order to preserve time
#  information; not use the updated time that comes back from run-pc.csh
#  in caminput? 


# Create an initial CAM.nc file from the DART state vector
   $PBS_O_WORKDIR/trans_sv_pv
   echo after trans_sv_pv
   ls -ltR 

# advance cam n hours; need generality
# this run-pc is resolution independent, and path relative (4/30/03)

   $PBS_O_WORKDIR/$model/models/atm/cam/bld/run-pc.csh $case $model $PBS_O_WORKDIR \
      > cam_out_temp

# Time not currently being advanced by model; need the target time
# as first line of updated state vector file
   head -1 temp_ic > temp_ud
# Generate the updated DART state vector
   $PBS_O_WORKDIR/trans_pv_sv
   echo after trans_pv_sv
   ls -lRt 

# For now, time is not being handled in model; just put state on
   tail +2 temp_ic >> temp_ud

mv temp_ud $PBS_O_WORKDIR/assim_model_state_ud$element; \
cd /scratch/local
rm -rf dartcam
