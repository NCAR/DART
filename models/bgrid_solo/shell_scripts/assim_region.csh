#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# Standard script for use in assimilation applications
# where regions are assimilated by separate executables.
# Can be used with most low-order models, the bgrid model 
# and most large models.

# This script copies the necessary files into the temporary directory
# and then executes the fortran program assim_region.

set PBS_O_WORKDIR = $1
set element = $2
set temp_dir = $3

# Standard script for use in assimilation applications
# where regions are to be assimilated by separate executables.

# Create a clean temporary directory and go there
rm -rf $temp_dir
mkdir  $temp_dir
cd     $temp_dir

# Copy the initial condition file to the temp directory
cp ${PBS_O_WORKDIR}/filter_assim_region__in$element filter_assim_region_in

# Copy the DART namelist to the temp directory
cp ${PBS_O_WORKDIR}/input.nml .

# Copy the assim_region executable to the temporary directory
cp ${PBS_O_WORKDIR}/assim_region .

# Copy the observation sequence file to the temporary directory
cp ${PBS_O_WORKDIR}/filter_assim_obs_seq .

# The original version of the bgrid model required the following.
# Hawaii and above versions without MPI do not.
#cp ${PBS_O_WORKDIR}/diag_table .
#mkdir RESTART

# Assimilate the region, saving standard out
./assim_region > assim_region_out_temp

# Append the output from the assimilation to the file in the working directory
cat assim_region_out_temp >> $PBS_O_WORKDIR/assim_region_out_temp$element

# Move the updated region state to the working directory
mv filter_assim_region_out $PBS_O_WORKDIR/filter_assim_region_out$element

# Change back to working directory and get rid of temporary directory
cd $PBS_O_WORKDIR
rm -rf $temp_dir
