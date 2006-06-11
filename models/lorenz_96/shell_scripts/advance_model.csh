#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section 
# University Corporation for Atmospheric Research
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

set WORKDIR = $1
set element = $2
set temp_dir = $3

set REMOVE = 'rm -rf'
set   COPY = 'cp -p'
set   MOVE = 'mv'

# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

# Create a clean temporary directory and go there
${REMOVE} $temp_dir
mkdir -p  $temp_dir
cd        $temp_dir

# Copy the initial condition file to the temp directory
${COPY} ${WORKDIR}/assim_model_state_ic$element temp_ic

# Copy the DART namelist to the temp directory
${COPY} ${WORKDIR}/input.nml .

# Copy the integrate_model executable to the temporary directory
${COPY} ${WORKDIR}/integrate_model .

# Advance the model, saving standard out
./integrate_model > integrate_model_out_temp

# Append the output from the advance to the file in the working directory
cat integrate_model_out_temp >> $WORKDIR/integrate_model_out_temp$element

# Move the updated state vector to the working directory
${MOVE} temp_ud $WORKDIR/assim_model_state_ud$element

# Change back to working directory and get rid of temporary directory
cd $WORKDIR
${REMOVE} temp_dir
