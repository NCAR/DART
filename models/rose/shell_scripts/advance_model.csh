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
${COPY} ${CENTRALDIR}/assim_model_state_ic$element temp_ic
echo "ls $temp_dir for element $element"
echo junk >! element$element
ls -lRt 

# Need a base "binary" restart file into which to copy modifications from filter
if (  -e   ${CENTRALDIR}/rose_restart_$element.dat) then
   ${COPY} ${CENTRALDIR}/rose_restart_$element.dat rose_restart.dat
   echo "Using newish rose_restart.dat file."
else
   echo "Using default rose_restart.dat file."
   ${COPY} ${CENTRALDIR}/NMC_SOC.day151_2002.dat rose_restart.dat
endif

${COPY} ${CENTRALDIR}/input.nml input.nml
${COPY} ${CENTRALDIR}/rose.nml rose.nml

# Get target_time from temp_ic
if (-e temp_ic && -e ${CENTRALDIR}/trans_time) then
   ${CENTRALDIR}/trans_time
   echo "after trans_time"
   ls -lRt
set ROSE_target_time = `cat times`
else
   echo 'no ic file or trans_time available for trans_time'
   exit 1
endif

# Create rose namelist
   ${COPY} rose.nml rose.nml_default
   set NML = namelist.in
   echo $ROSE_target_time  
   echo $ROSE_target_time >! $NML
   ${CENTRALDIR}/nmlbld_rose < $NML
   echo "after nmlbld_rose"
   ls -ltR

# Create an initial rose_restart file from the DART state vector
   echo "Running trans_sv_pv at "`date`
   ${CENTRALDIR}/trans_sv_pv
   echo "after trans_sv_pv"
   ls -ltR 

# advance ROSE "target_time" hours
   echo "Running rose at "`date`
   ${CENTRALDIR}/rose >! rose_out_temp
   echo "Finished rose at "`date`

# Generate the updated DART state vector
   ${CENTRALDIR}/trans_pv_sv
   echo "after trans_pv_sv"
   ls -lRt 

${COPY} temp_ud          ${CENTRALDIR}/assim_model_state_ud$element
${COPY} rose_restart.dat ${CENTRALDIR}/rose_restart_$element.dat
${COPY} rose.nml         ${CENTRALDIR}

# Change back to working directory and get rid of temporary directory
cd ${CENTRALDIR}
#${REMOVE} ${temp_dir}
