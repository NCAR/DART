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
cp ${PBS_O_WORKDIR}/assim_region .

# Copy the observation sequence file to the temporary directory
cp ${PBS_O_WORKDIR}/filter_assim_obs_seq .

echo ls $temp_dir for element $element
echo junk > element$element
ls -lRt 

cp ${PBS_O_WORKDIR}/input.nml input.nml
cp ${PBS_O_WORKDIR}/rose.nml rose.nml

   ./assim_region > rose_out_temp
   echo after executing assim_region

echo element $element >> rose_reg_temp
ls -lt >> rose_reg_temp
cat rose_reg_temp >> $PBS_O_WORKDIR/rose_reg_temp$element
mv filter_assim_region_out $PBS_O_WORKDIR/filter_assim_region_out$element

# Change back to working directory and get rid of temporary directory
cd $PBS_O_WORKDIR
rm -rf $temp_dir
