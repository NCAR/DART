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
# Shell script to work with asynchronous filter integration
# This script needs to be piped to the filter program with the
# filter namelist async variable set to .true.

#--------------------------------------------------------------------
# If this is first of recursive calls need to get rid of async_may_go
# Technically, this could lock, but seems incredibly unlikely
if ($?First) then
   setenv First no
   rm -f async_may_go

# Copy an initial condition for cam .nc files in cam working directory
   cp /home/jla/cam_test/cam/T5iter/caminput1.nc /home/jla/cam_test/cam/T5iter/caminput.nc
   cp /home/jla/cam_test/cam/T5iter/clminput1.nc /home/jla/cam_test/cam/T5iter/clminput.nc

endif
#--------------------------------------------------------------------


while(1 == 1)
   ls async_may_go > .async_garb
   if($status == 0) break
   echo waiting_for_async_may_go_file
   sleep 5
end
echo found_async_may_go_file

# First line of filter_control should have number of model states to be integrated
set num = `head -1 filter_control`


# Ready for advance; make a copy of base cam and clm .nc files
# Should really keep a separate one for each ensemble member but not yet
cp /home/jla/cam_test/cam/T5iter/caminput.nc CAM_FILE.nc
# Also need a base land file; NOTE REALLY NEED TO KEEP SEPARATE FOR EACH ENS
cp /home/jla/cam_test/cam/T5iter/clminput.nc CLM_FILE.nc


# Create a directory for each member to run in for namelists
set element = 1
while($element <= $num)
# Make a temporary directory for this element's run
   mkdir tempdir$element
# Copy the initial condition file to the temp directory
# Need to strip out the current time (second line) and 
# leave advance to time (first line); Should all be 
# automated with more generalized model time handling
   head -1 assim_model_state_ic$element > tempdir$element/temp_ic
   tail +3 assim_model_state_ic$element >> tempdir$element/temp_ic
   rm -f assim_model_state_ic$element
# Need template for cam model_mod to get resolution information
   cp T5H0-12icl.cam2.i.0001-09-01-43200.nc tempdir$element
# Need a base nc file into which to copy modifications from filter
   cp CAM_FILE.nc tempdir$element
# Also need a base land file; NOTE REALLY NEED TO KEEP SEPARATE FOR EACH ENS
   cp CLM_FILE.nc tempdir$element
   cd tempdir$element
# Create an initial CAM.nc file from the DART state vector
   ../trans_sv_pv
# Copy the CAM and CLM initial files to cam2.0.1 hierarchy for run
   cp CAM_FILE.nc /home/jla/cam_test/cam/T5iter/caminput.nc
   cp CLM_FILE.nc /home/jla/cam_test/cam/T5iter/clminput.nc
# Run the 12 hour advance for the cam; need generality
   /home/jla/cam_test/cam2.0.1/models/atm/cam/bld/run-pc.csh > cam_out_temp
echo DONE_RUNNING_CAM
# Need to remove accumulating h0, h1, and .r. cam files
   rm -f /home/jla/cam_test/cam/T5iter/T5iter*.r.*
   rm -f /home/jla/cam_test/cam/T5iter/*.h0.*.nc
   rm -f /home/jla/cam_test/cam/T5iter/*.h1.*.nc
# Move back up a level and do next model; might want to do these sequentially
   cd ..



#------------------ Block removed for sequential application ----
# Needed for memory conservation on anchorage root node
#   @ element++



#end

# Wait for all processes to finish up
#wait

# All model runs complete, move the updated file up
#set element = 1
#while($element <= $num)
#------------------------------------------------------------


echo JUST_PAST_THE_WAIT_BLOCK

   cd tempdir$element
# Get the updated cam init file
   cp /home/jla/cam_test/cam/T5iter/caminput.nc CAM_FILE.nc

# Time not currently being advanced by model; need the target time
# as first line of updated state vector file
   head -1 temp_ic > ../assim_model_state_ud$element
# Generate the updated DART state vector
   ../trans_pv_sv
   cd ..

# For now, time is not being handled in model; just put state on
   tail +2 tempdir$element/temp_ic >> assim_model_state_ud$element

# Clear out the working directory
   rm -rf tempdir$element
   @ element++
end

# Remove the semaphore file
rm -f async_may_go

# Cleaned up; let the filter know it can proceed
echo All_done:Please_proceed

# Doing recursive call
csh ./cam_async.csh
