#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
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
cp $PBS_O_WORKDIR/assim_model_state_ic$element temp_ic
# Need a base nc file into which to copy modifications from filter
cp $PBS_O_WORKDIR/caminput.nc .
# CLM files used without modification
if (-e $PBS_O_WORKDIR/clminput_$element.nc) then
  cp $PBS_O_WORKDIR/clminput_$element.nc clminput.nc
else
  cp $PBS_O_WORKDIR/clminput.nc .
endif
cp $PBS_O_WORKDIR/input.nml input.nml

echo ls /scratch/local/dartcam for element $element
echo junk > element$element
ls -lRt 

# Copy the initial condition file to the temp directory
# Need to strip out the current time (second line) and 
# leave advance to time (first line); Should all be 
# automated with more generalized model time handling
set filetype = `file temp_ic`
echo temp_ic has file type $filetype[2]
if ($filetype[2] == ASCII) then
   head -1 temp_ic > temp2
   tail +3 temp_ic >> temp2
   mv temp2 temp_ic
else if ($filetype[2] == data) then
# binary version
   dd if=temp_ic of=time1.bin bs=4 count=4
   dd if=temp_ic of=data1.bin bs=4 skip=8
   cat time1.bin data1.bin >! temp_ic
else
   echo filetype of initial temp_ic not recognized
   stop
endif

#  is some of this caminput/CAM_FILE copying in order to preserve time
#  information; not use the updated time that comes back from run-pc.csh
#  in caminput? 


# Create an initial CAM.nc file from the DART state vector

$PBS_O_WORKDIR/trans_sv_pv

# test model advance
ls -ltR 

# advance cam n hours; need generality
# this run-pc is resolution independent, and path relative (4/30/03)

#if ($element == 1) then
#   $PBS_O_WORKDIR/$model/models/atm/cam/bld/run-pc.nonamelist.csh $case $model $PBS_O_WORKDIR \
#      > cam_out_temp
#else
   $PBS_O_WORKDIR/$model/models/atm/cam/bld/run-pc.csh $case $model $PBS_O_WORKDIR \
      > cam_out_temp
#endif

echo $element >> $PBS_O_WORKDIR/dump
tail cam_out_temp >> $PBS_O_WORKDIR/dump
cp cam_out_temp $PBS_O_WORKDIR/cam_out_temp$element

# Time not currently being advanced by model; need the target time (first of 2 time lines
# in _ic# files) as first line of updated state vector file
if ($filetype[2] == ASCII) then
   head -1 temp_ic > temp_ud
else if ($filetype[2] == data) then
   dd if=temp_ic of=time2.bin bs=4 count=4
else
   echo filetype of temp_ic not recognized
   stop
endif

# Generate the updated DART state vector and put it in temp_ic (time followed by state)
   $PBS_O_WORKDIR/trans_pv_sv

# For now, time is not being handled in model; just put state on
if ($filetype[2] == ASCII) then
   tail +2 temp_ic >> temp_ud
else if ($filetype[2] == data) then
   dd if=temp_ic of=data2.bin bs=4 skip=4
   cat time2.bin data2.bin >! temp_ud
endif
ls -lRt 

mv temp_ud $PBS_O_WORKDIR/assim_model_state_ud$element; \
mv clminput.nc $PBS_O_WORKDIR/clminput_$element.nc
cd /scratch/local
rm -rf dartcam
