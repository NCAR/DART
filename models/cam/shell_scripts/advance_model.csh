#!/bin/csh

# set echo verbose

set PBS_O_WORKDIR = $1
set element = $2
set temp_dir = $3

# set $case = the case and $model = the model we're running
set caseinfo = `cat $PBS_O_WORKDIR/casemodel`
set case = $caseinfo[1]
set model = $caseinfo[2]
echo $case $model

# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

# kdr this script is executed entirely in temp_dir, pulling files from elsewhere
#     except for the build in run-pc.csh, which executes in $blddir,
#     if I haven't created the executable cam in setup_advance_model.csh,
#     or init_adv....

mkdir $temp_dir
cd $temp_dir
cp $PBS_O_WORKDIR/assim_model_state_ic$element temp_ic
cp $PBS_O_WORKDIR/input.nml input.nml
echo ls $temp_dir for element $element
echo junk > element$element
ls -lRt 

# Need a base nc file into which to copy modifications from filter
cp $PBS_O_WORKDIR/caminput.nc .
# cp $PBS_O_WORKDIR/clminput.nc .
if (-e $PBS_O_WORKDIR/clminput_$element.nc) then
  cp $PBS_O_WORKDIR/clminput_$element.nc clminput.nc
else
  cp $PBS_O_WORKDIR/clminput.nc .
endif
cp $PBS_O_WORKDIR/input.nml input.nml


# Copy the initial condition file to the temp directory
# Need to strip out the current time (second line) and 
# leave advance to time (first line); Should all be 
# automated with more generalized model time handling


# Do we do this in the script instead of trans_sv_pv because of the form
# of temp_ic; the 2 times at the beginning are hard to read?
# Try moving this stuff into trans_sv_pv


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
   exit 1
endif

# Create an initial CAM.nc file from the DART state vector
   $PBS_O_WORKDIR/trans_sv_pv
   ls -ltR 

# advance cam 
   $PBS_O_WORKDIR/$model/models/atm/cam/bld/run-pc.csh $case $model $PBS_O_WORKDIR \
      > cam_out_temp
echo $element >> $PBS_O_WORKDIR/dump
tail cam_out_temp >> $PBS_O_WORKDIR/dump
cp cam_out_temp $PBS_O_WORKDIR/cam_out_temp$element

# Time not currently being advanced by model; need the target time 
# (first of 2 time lines in _ic# ) as first line of updated state vector file

# Do we do this in the script instead of trans_sv_pv because of the form
# of temp_ic; the 2 times at the beginning are hard to read?
# Try moving this stuff into trans_sv_pv

if ($filetype[2] == ASCII) then
   head -1 temp_ic > temp_ud
else if ($filetype[2] == data) then
   dd if=temp_ic of=time2.bin bs=4 count=4
else
   echo filetype of temp_ic not recognized
   exit 1
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
cd $PBS_O_WORKDIR
rm -rf $temp_dir
