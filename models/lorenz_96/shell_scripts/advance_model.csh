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
# where the model advance is executed as a separate process.
# Can be used with most low-order models and the bgrid model which
# can be advanced using the integrate_model executable.

# This script copies the necessary files into the temporary directory
# and then executes the fortran program integrate_model.

set PBS_O_WORKDIR = $1
set element = $2
set temp_dir = $3

# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

# Create a clean temporary directory and go there
rm -rf $temp_dir
mkdir  $temp_dir
cd     $temp_dir

# Copy the initial condition file to the temp directory
cp ${PBS_O_WORKDIR}/assim_model_state_ic$element temp_ic

# Copy the DART namelist to the temp directory
cp ${PBS_O_WORKDIR}/input.nml .

# Copy the integrate_model executable to the temporary directory
cp ${PBS_O_WORKDIR}/integrate_model .

# Advance the model, saving standard out
./integrate_model > integrate_model_out_temp

# Append the output from the advance to the file in the working directory
cat integrate_model_out_temp >> $PBS_O_WORKDIR/integrate_model_out_temp$element

# Move the updated state vector to the working directory
mv temp_ud $PBS_O_WORKDIR/assim_model_state_ud$element

# Change back to working directory and get rid of temporary directory
cd $PBS_O_WORKDIR
rm -rf $temp_dir
