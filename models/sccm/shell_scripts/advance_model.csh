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
# Script to advance the sccm model. Requires converting dart format
# files into sccm input format, running sccm, and converting files back
# to dart format. The input dart files have a target time pre-pended.
# The files going back to dart do not have a target time. 
# The programs trans_dart_to_sccm and trans_sccm_to_dart handle the
# translation.

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
cp ${PBS_O_WORKDIR}/assim_model_state_ic$element dart_file_in

# Copy the DART namelist to the temp directory
cp ${PBS_O_WORKDIR}/input.nml .

# Copy the trans executables to the temporary directory
cp ${PBS_O_WORKDIR}/trans_dart_to_sccm .
cp ${PBS_O_WORKDIR}/trans_sccm_to_dart .

# Yuqiong PLEASE COPY YOUR MODEL PROGRAM, SCRIPTS and ANY NEEDED
# FILES TO THIS TEMPORARY DIRETORY HERE

# Translate the dart file to an sccm input file (dart_data.dat)
./trans_dart_to_sccm

# Yuqiong: Advance the model (advanced state overwrites dart_data.dat)
# ADVANCE MODEL COMMAND GOES HERE

# Append the output from the advance to the file in the working directory
# NEED TO DO THIS TO KEEP RECORD OF RUN. Yuqiong, keep your output somewher
# if you need it.

# Convert the sccm file to a dart file (ends up in dart_file_out)
./trans_sccm_to_dart

# Move the updated state vector to the working directory
mv dart_file_out $PBS_O_WORKDIR/assim_model_state_ud$element

# Change back to working directory and get rid of temporary directory
cd $PBS_O_WORKDIR
rm -rf $temp_dir
