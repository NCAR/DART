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

mkdir $temp_dir
cd $temp_dir
cp $PBS_O_WORKDIR/assim_model_state_ic$element temp_ic

cp $PBS_O_WORKDIR/input.nml input.nml
echo ls $temp_dir for element $element
echo junk > element$element
ls -lRt 

# Need a base nc file into which to copy modifications from filter
# c[al]minput_$element carries along CAM/CLM model state which is not updated
#      by the filter

if (-e $PBS_O_WORKDIR/caminput_$element.nc) then
   cp $PBS_O_WORKDIR/caminput_$element.nc caminput.nc
else
   cp $PBS_O_WORKDIR/caminput.nc .
endif
if (-e $PBS_O_WORKDIR/clminput_$element.nc) then
   cp $PBS_O_WORKDIR/clminput_$element.nc clminput.nc
else
   cp $PBS_O_WORKDIR/clminput.nc .
endif
cp $PBS_O_WORKDIR/input.nml input.nml

# create 'times' file from DART dates in assim_model_state_ic1
if (-e temp_ic && -e $PBS_O_WORKDIR/trans_time) then
   echo 'advance_model; executing trans_time'
   $PBS_O_WORKDIR/trans_time
   ls -lt 
   cp times $PBS_O_WORKDIR
else
   echo 'no ic file or trans_time available for trans_time'
   exit 1
endif

# Create an initial CAM.nc file from the DART state vector
# Times are handled separately in trans_time
$PBS_O_WORKDIR/trans_sv_pv
ls -ltR 

# advance cam 
$PBS_O_WORKDIR/$model/models/atm/cam/bld/run-pc.csh $case $model $PBS_O_WORKDIR \
   > cam_out_temp
echo $element >> $PBS_O_WORKDIR/dump
tail cam_out_temp >> $PBS_O_WORKDIR/dump
cp cam_out_temp $PBS_O_WORKDIR/cam_out_temp$element

# Generate the updated DART state vector and put it in temp_ic (time followed by state)
   $PBS_O_WORKDIR/trans_pv_sv

mv temp_ud $PBS_O_WORKDIR/assim_model_state_ud$element; \
mv clminput.nc $PBS_O_WORKDIR/clminput_$element.nc
mv caminput.nc $PBS_O_WORKDIR/caminput_$element.nc
mv namelist $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
rm -rf $temp_dir
