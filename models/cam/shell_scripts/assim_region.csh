#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# Standard script for use in assimilation applications
# where regions are assimilated by separate executables.

# This script copies the necessary files into the temporary directory
# and then executes the fortran program assim_region.

set PBS_O_WORKDIR = $1
set element = $2
set temp_dir = $3

rm -rf $temp_dir
mkdir  $temp_dir
cd     $temp_dir

# Copy the executable, initial condition file, and other inputs to the temp directory

cp ${PBS_O_WORKDIR}/assim_region .
cp ${PBS_O_WORKDIR}/filter_assim_region__in$element filter_assim_region_in
cp ${PBS_O_WORKDIR}/input.nml .
cp ${PBS_O_WORKDIR}/filter_assim_obs_seq .
cp $PBS_O_WORKDIR/caminput.nc .
cp $PBS_O_WORKDIR/clminput.nc .

./assim_region > cam_reg_temp
# writes out filter_assim_region_out
# for filter...assim_tools_mod/filter_assim() to read in with
# ensemble member #s tacked onto the end

echo element $element >> cam_reg_temp
ls -lt >> cam_reg_temp

cat cam_reg_temp >> $PBS_O_WORKDIR/cam_reg_temp$element

mv filter_assim_region_out $PBS_O_WORKDIR/filter_assim_region_out$element

cd $PBS_O_WORKDIR
rm -rf $temp_dir
