#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

# This script copies the necessary files into the temporary directory
# for a model run. It assumes that there is ${PBS_O_WORKDIR}/WRF directory
# where boundary conditions and namelist files reside.

set PBS_O_WORKDIR = $1
set element = $2
set temp_dir = $3

# Shell script to run the WRF model from DART input.

set verbose

rm -rf $temp_dir
mkdir  $temp_dir
cd     $temp_dir

# Copy the initial condition file to the temp directory

cp ${PBS_O_WORKDIR}/wrfinput .
mv ${PBS_O_WORKDIR}/assim_model_state_ic$element dart_wrf_vector
ln -s ${PBS_O_WORKDIR}/input.nml .

#set time = `cat ${PBS_O_WORKDIR}/async_may_go`
#set filetype = `file dart_wrf_vector`
#if ($filetype[2] == ASCII) then
#   set time = `head -1 dart_wrf_vector`
#else if ($filetype[2] == data) then
#   set time = `dd if=dart_wrf_vector bs=4 count=4`
#else
#   echo filetype of initial dart_wrf_vector not recognized
#   stop
#endif

set time = `tail -1 ${PBS_O_WORKDIR}/filter_control`

set secs = $time[1]
set days = $time[2]

# Copy the boundary condition file to the temp directory.
cp ${PBS_O_WORKDIR}/WRF/wrfbdy_${days}_${secs}_$element wrfbdy_d01

# Copy WRF input namelist to the temp directory.
ln -s  ${PBS_O_WORKDIR}/WRF/namelist.input_${days}_${secs}_$element namelist.input
ln -s  ${PBS_O_WORKDIR}/RRTM_DATA .
ln -s  ${PBS_O_WORKDIR}/LANDUSE.TBL .

# Convert DART to wrfinput

echo ".true." | ${PBS_O_WORKDIR}/dart_tf_wrf >& out.dart_to_wrf

#${PBS_O_WORKDIR}/trans_time >& out.trans_time

mv wrfinput wrfinput_d01

# Update boundary conditions

${PBS_O_WORKDIR}/update_wrf_bc >& out.update_wrf_bc

${PBS_O_WORKDIR}/wrf.exe >>& out_wrf_integration
mv wrfout_d01_000000 wrfinput

mv dart_wrf_vector dart_wrf_vector.input

# create new input to DART (taken from "wrfinput")
echo ".false." | ${PBS_O_WORKDIR}/dart_tf_wrf >& out.wrf_to_dart

mv dart_wrf_vector $PBS_O_WORKDIR/assim_model_state_ud$element

cd $PBS_O_WORKDIR
#rm -rf $temp_dir

exit
