#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

# This script copies the necessary files into the temporary directory
# and then executes the fortran program integrate_model.

set PBS_O_WORKDIR = $1
set element = $2
set temp_dir = $3

# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

rm -rf $temp_dir
mkdir  $temp_dir
cd     $temp_dir

# Copy the initial condition file to the temp directory

cp ${PBS_O_WORKDIR}/assim_model_state_ic$element temp_ic
cp ${PBS_O_WORKDIR}/input.nml .
cp ${PBS_O_WORKDIR}/integrate_model .

./integrate_model > lorenz_96_out_temp

cat lorenz_96_out_temp >> $PBS_O_WORKDIR/lorenz_96_out_temp$element

mv temp_ud $PBS_O_WORKDIR/assim_model_state_ud$element

cd $PBS_O_WORKDIR
rm -rf $temp_dir
