#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
# 
# This script advances the MITgcm annulus model.  It copies the necessary
# files into the relevant directory, converts the dart output files to
# an MITgcm pickup file, runs the model, converts the MITgcm pickup
# file back into a dart file, and sends the file back to the dart directory.

# The script advance_ens.csh passes a working directory, an ensemble
# member number, and a temp directory.  First thing to do is place those
# passed values into variables

set  CENTRALDIR = $1
set     element = $2
set    temp_dir = $3

# Lots of logging
#set verbose

# The MITgcm configuration has directories associated with the ensemble
# member number, thus single digit numbers must be prefixed with a 0.
if ($element == 1) set ensdir = 01
if ($element == 2) set ensdir = 02
if ($element == 3) set ensdir = 03
if ($element == 4) set ensdir = 04
if ($element == 5) set ensdir = 05
if ($element == 6) set ensdir = 06
if ($element == 7) set ensdir = 07
if ($element == 8) set ensdir = 08
if ($element == 9) set ensdir = 09
if ($element >  9) set ensdir = $element

# Set the path to where the MITgcm works its magic, and go to that directory
set MITGCM = ${CENTRALDIR}/MITgcm/verification/osse/da/${ensdir}/assimilate
cd $MITGCM

# Move the dart initial condition file to the MITgcm directory
mv ${CENTRALDIR}/assim_model_state_ic$element dart_vector

# Link in the input.nml namelist for use by translation programs
ln -s ${CENTRALDIR}/input.nml .

# Do a byteswap on the old restart files so can be read by trans code
${CENTRALDIR}/byteswap pickup.in pickup.in.s -w 8
${CENTRALDIR}/byteswap pickup_nh.in pickup_nh.in.s -w 8

# Convert the dart initial condition file to a MITgcm pickup file
${CENTRALDIR}/trans_dart_to_MITgcm

# Do a byteswap on the new restart files so MITgcm can read them
${CENTRALDIR}/byteswap pickup.out.s pickup.out -w 8
${CENTRALDIR}/byteswap pickup_nh.out.s pickup_nh.out -w 8

# Run the MITgcm
sleep 0.05
echo "Running the MITgcm on ensemble member "
echo $element
./mitgcmuv >& logfile
sleep 0.05

# Do a byteswap on the new restart files so can be read by trans code
${CENTRALDIR}/byteswap pickup.in pickup.in.s -w 8
${CENTRALDIR}/byteswap pickup_nh.in pickup_nh.in.s -w 8

# Remove the old dart_vector file
\rm -f dart_vector

# Convert the resulting MITgcm pickup file to a dart file
${CENTRALDIR}/trans_MITgcm_to_dart

# Move the dart file back to the working directory
sleep 0.05
mv dart_vector ${CENTRALDIR}/assim_model_state_ud$element
sleep 0.05

exit 0


