#!/bin/csh
# Shell script to work with asynchronous filter integration
# This script needs to be piped to the filter program with the
# filter namelist async variable set to .true.

# If this is first of recursive calls need to get rid of async_may_go
# Technically, this could lock, but seems incredibly unlikely
if ($?First) then
   setenv First no
   rm -f async_may_go
endif

while(1 == 1)
   ls async_may_go > .async_garb
   if($status == 0) break
   echo waiting_for_async_may_go_file
   sleep 5
end
echo found_async_may_go_file

# First line of filter_control should have number of model states to be integrated
set num = `head -1 filter_control`

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
   cp T5H0-12icl.cam2.i.0001-09-01-43200.nc tempdir$element
   cp T5H0-12icl.cam2.i.0001-09-01-43200.nc tempdir$element/CAM_FILE.nc
   cd tempdir$element
   ../trans_sv_pv
   cd ..
   @ element++
end

# Wait for all processes to finish up
#wait

# All model runs complete, move the updated file up
set element = 1
while($element <= $num)
#   mv tempdir$element/temp_ud assim_model_state_ud$element
   cd tempdir$element
   ../trans_pv_sv
   cd ..
   mv tempdir$element/temp_ic assim_model_state_ud$element
   rm -rf tempdir$element
   @ element++
end

# Remove the semaphore file
rm -f async_may_go

# Cleaned up; let the filter know it can proceed
echo All_done:Please_proceed

# Doing recursive call
csh ./cam_async.csh
