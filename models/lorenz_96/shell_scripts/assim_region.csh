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

# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

rm -rf $temp_dir
mkdir  $temp_dir
cd     $temp_dir

# Copy the initial condition file to the temp directory

cp ${PBS_O_WORKDIR}/filter_assim_region__in$element filter_assim_region_in
cp ${PBS_O_WORKDIR}/input.nml .
cp ${PBS_O_WORKDIR}/assim_region .
cp ${PBS_O_WORKDIR}/filter_assim_obs_seq .

./assim_region > lorenz_96_out_temp

cat lorenz_96_out_temp >> $PBS_O_WORKDIR/lorenz_96_out_temp$element

mv filter_assim_region_out $PBS_O_WORKDIR/filter_assim_region_out$element

cd $PBS_O_WORKDIR
rm -rf $temp_dir
