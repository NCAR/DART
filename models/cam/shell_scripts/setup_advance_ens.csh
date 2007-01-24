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
# || marks changes for || (and batch) setup.

# Script to guide in getting required files for running
# a cam model DART assimilation starting in the PWD.

# Assumes execution from .../cam_test

# this is written for anchorage [dart]

# Should I have a call to run-pc.csh in here to do the 0th iteration;
# create the first caminput.nc and clminput.nc, using original files in /fs/cgd
# Then cam would exist for all the "elements", and I wouldn't have to 
# manually create c[al]minput.nc and namelists before running the filter.
#   YES, but later
# Or would it work better in init_advance_model.csh?
#   NO

# set $case = the case and $model = the model we're running
set case = T5
set model = cam2.0.1
# output info for advance_model.csh to use
echo $case $model > casemodel
# set caseinfo = `cat casemodel`

# Set directory for finding source code
set SOURCE_DIR =  ~raeder/CamPathRel/Source
# set SOURCE_DIR =  ~jla/cam_test/work

# Set directory for scripts, model dependent and independent
set SCRIPT_DIR = ~raeder/CamPathRel/Scripts
# set SCRIPT_DIR = ~raeder/CamPathRel

# Set directory for cam and case input
set DATA_DIR = ~raeder/CamPathRel/Caminput/$case

# Copy required files to current directory to run
cp $SOURCE_DIR/filter .
cp $SOURCE_DIR/perfect_model_obs .
cp $SOURCE_DIR/trans_pv_sv .
cp $SOURCE_DIR/trans_sv_pv .
cp $SOURCE_DIR/trans_time .
# cp /home/jla/execs/cam_filter ./filter
# cp /home/jla/execs/cam_perfect_model_obs ./perfect_model_obs
# cp /home/jla/execs/trans_pv_sv .
# cp /home/jla/execs/trans_sv_pv .

# Copy the obs_seq.in file and namelist
cp $SOURCE_DIR/obs_seq.in .
cp $SOURCE_DIR/input.nml .

# If needed, copy some filter_ics and perfect_model_ics
# Make sure destinations are consistent with namelist
cp $DATA_DIR/cam_filter_ics .
cp $DATA_DIR/cam_perfect_ics .
# cp $SOURCE_DIR/cam_perfect_restart.base  cam_perfect_ics
# I deleted a line here; similar but for filter

# Copy the DART scripts needed to run 
# [cp ~jla/DART/shell_scripts/async_filter.csh .]
# [cp ~jla/DART/shell_scripts/long_run.csh .]
# [cp ~jla/DART/shell_scripts/async_long_run.csh .]
# kdr; on anchorage these would have come in 'cp -Rf ~jla/cam_test'
# || cp $SCRIPT_DIR/async_filter.csh .    
cp $SCRIPT_DIR/sync_submit.csh .    
cp $SCRIPT_DIR/async_long_run.csh .  
 
# Copy all the other required scripts for cam model
cp $SCRIPT_DIR/init_advance_model.csh .
cp $SCRIPT_DIR/advance_model.csh .
# || This is the batch script submitted by qsub in sync_submit.csh
cp $SCRIPT_DIR/advance_ens.csh .

# Copy an initial condition for cam .nc files in cam working directory
# unnecessary copy?; just reference the permanent storage area where needed?
# No, we want all this changable stuff in setup_advance_, not in advance_model
   cp $DATA_DIR/caminput.nc caminput.nc
   cp $DATA_DIR/clminput.nc clminput.nc
#   cp ~jla/cam_test/cam/$case/caminput1.nc caminput.nc
#   cp ~jla/cam_test/cam/$case/clminput1.nc clminput.nc

# copy cam directory of source code and scripts into cwd
cp -r $SOURCE_DIR/$model .
# cp -r ~jla/cam_test/cam2.0.1/ .

# copy case specific cam script input 
cp           $DATA_DIR/config_cache_defaults.xml \
   $model/models/atm/cam/bld/config_cache_defaults.xml
cp $DATA_DIR/namelist   namelistin


