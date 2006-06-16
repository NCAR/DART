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
# Script to advance the sccm model. Requires converting dart format
# files into sccm input format, running sccm, and converting files back
# to dart format. The input dart files have a target time pre-pended.
# The files going back to dart do not have a target time. 
# The programs trans_dart_to_sccm and trans_sccm_to_dart handle the
# translation.

# This script copies the necessary files into the temporary directory
# and then executes the fortran program integrate_model.

set      myname = $0
set  CENTRALDIR = $1
set     element = $2
set    temp_dir = $3

set REMOVE = 'rm -rf'
set   COPY = 'cp -p'
set   MOVE = 'mv -f'

# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

echo "starting ${myname} for ens member $element at "`date`
echo "CENTRALDIR is ${CENTRALDIR}"
echo "temp_dir is ${temp_dir}"

# Create a clean temporary directory and go there
${REMOVE} ${temp_dir}
mkdir -p  ${temp_dir}
cd        ${temp_dir}

# Copy the initial condition file to the temp directory
${COPY} ${CENTRALDIR}/assim_model_state_ic$element dart_file_in

# Copy the DART namelist to the temp directory
${COPY} ${CENTRALDIR}/input.nml .

# Copy the trans executables to the temporary directory
${COPY} ${CENTRALDIR}/trans_dart_to_sccm .
${COPY} ${CENTRALDIR}/trans_sccm_to_dart .

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
${MOVE} dart_file_out ${CENTRALDIR}/assim_model_state_ud$element

# Change back to working directory and get rid of temporary directory
cd ${CENTRALDIR}
#${REMOVE} ${temp_dir}
