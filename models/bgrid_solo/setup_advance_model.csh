#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
# 
# Script to guide in getting required files for running
# a b-grid model DART assimilation starting in the PWD.
# Might want to make this more modification friendly.

# Set directory for finding source code
set SOURCE_DIR = ~jla/fms_havana/work


# Copy required files to current directory to run
cp $SOURCE_DIR/filter .
cp $SOURCE_DIR/perfect_model_obs .
cp $SOURCE_DIR/integrate_model .

# Copy the obs_seq.in file
cp ~jla/fms_havana/work/obs_seq.in .

# Copy b-grid required files
cp -r ~jla/fms_havana/work/RESTART .
cp ~jla/fms_havana/work/diag_table .
cp ~jla/fms_havana/work/input.nml .

# If needed, copy some filter_ics and perfect_model_ics
# Make sure destinations are consistent with namelist
cp /data/jla/exp38/filter_ics_base filter_ics
cp /data/jla/exp38/perfect_ics_base  perfect_ics

# Copy the DART scripts needed to run 
cp ~jla/DART/shell_scripts/async_filter.csh .
cp ~jla/DART/shell_scripts/long_run.csh .
cp ~jla/DART/shell_scripts/async_long_run.csh .

# Copy all the other required scripts for bgrid model
cp ~jla/DART/models/bgrid_solo/advance_model.csh .
cp ~jla/DART/models/bgrid_solo/init_advance_model.csh .

