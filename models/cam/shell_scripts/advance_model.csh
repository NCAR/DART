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

# create 'times' file from DART dates in assim_model_state_ic1
if (-e temp_ic && -e $PBS_O_WORKDIR/trans_time) then
   echo 'advance_model; exectuting trans_time'
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

# Time not currently being advanced by model; need the target time 
# (first of 2 time lines in _ic# ) as first line of updated state vector file

# Generate the updated DART state vector and put it in temp_ic (time followed by state)
   $PBS_O_WORKDIR/trans_pv_sv

mv temp_ud $PBS_O_WORKDIR/assim_model_state_ud$element; \
mv clminput.nc $PBS_O_WORKDIR/clminput_$element.nc
cd $PBS_O_WORKDIR
rm -rf $temp_dir
