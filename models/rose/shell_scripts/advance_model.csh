#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.htm
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$

# set echo verbose

set PBS_O_WORKDIR = $1
set element = $2
set temp_dir = $3

# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

mkdir $temp_dir
cd $temp_dir

cp $PBS_O_WORKDIR/assim_model_state_ic$element temp_ic
echo ls $temp_dir for element $element
echo junk > element$element
ls -lRt 

# Need a base "binary" restart file into which to copy modifications from filter
if (-e $PBS_O_WORKDIR/rose_restart_$element.dat) then
   cp $PBS_O_WORKDIR/rose_restart_$element.dat rose_restart.dat
   set ROSE_old_restart = .false.
   set ROSE_nstart = 1
else
   cp $PBS_O_WORKDIR/NMC_SOC.day151_2002.dat rose_restart.dat
   set ROSE_old_restart = .true.
   set ROSE_nstart = 0
endif

cp $PBS_O_WORKDIR/input.nml input.nml
cp $PBS_O_WORKDIR/rose.nml rose.nml

# Get target_time from temp_ic
if (-e temp_ic && -e $PBS_O_WORKDIR/trans_time) then
   $PBS_O_WORKDIR/trans_time
   echo after trans_time
   ls -lRt
set ROSE_target_time = `cat times`
else
   echo 'no ic file or trans_time available for trans_time'
   exit 1
endif

# Create rose namelist
   cp rose.nml rose.nml_default
   set NML = namelist.in
   echo $ROSE_old_restart 
   echo $ROSE_old_restart   > $NML
   echo $ROSE_nstart      
   echo $ROSE_nstart       >> $NML
   echo $ROSE_target_time  
   echo $ROSE_target_time  >> $NML
   $PBS_O_WORKDIR/nmlbld_rose < $NML
   echo after nmlbld_rose
   ls -ltR

# Create an initial rose_restart file from the DART state vector
   $PBS_O_WORKDIR/trans_sv_pv
   echo after trans_sv_pv
   ls -ltR 

# advance ROSE "target_time" hours
   $PBS_O_WORKDIR/rose > rose_out_temp
   echo after executing rose

# Generate the updated DART state vector
# Needs to change "old_restart" switch in rose.nml
   cp rose.nml rose.nml_default
   set ROSE_old_restart = .false.
   set ROSE_nstart = 1
   echo $ROSE_old_restart  
   echo $ROSE_old_restart   > $NML
   echo $ROSE_nstart       
   echo $ROSE_nstart       >> $NML
   echo $ROSE_target_time  
   echo $ROSE_target_time  >> $NML
   $PBS_O_WORKDIR/nmlbld_rose < $NML
   echo after nmlbld_rose
   ls -ltR

   $PBS_O_WORKDIR/trans_pv_sv
   echo after trans_pv_sv
   ls -lRt 

mv temp_ud $PBS_O_WORKDIR/assim_model_state_ud$element; \
mv rose_restart.dat $PBS_O_WORKDIR/rose_restart_$element.dat
mv rose.nml $PBS_O_WORKDIR

cd $PBS_O_WORKDIR
rm -rf $temp_dir
