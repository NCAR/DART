#!/bin/csh
# Shell script to work with asynchronous filter integration
# This script needs to be piped to the filter program with the
# filter namelist async variable set to .true.

# If this is first of recursive calls need to get rid of async_may_go
# Technically, this could lock, but seems incredibly unlikely
if ($?First) then
# Do nothing if this is not the first time
else
   setenv First no
###   rm -f async_may_go
# Clean up any assim_model_ic and ud files and temp directories
   rm -f assim_model_ic*
   rm -f assim_model_ud*
   rm -rf tempdir*
# Call the model's initialization script to allow it to set up if needed
   csh ./init_advance_model.csh
endif

while(1 == 1)
   rm -f .async_garb
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
# Copy the appropriate ics file to the temp directory
   mv assim_model_state_ic$element tempdir$element/temp_ic
# Change to temp directory and run integrate
   cd tempdir$element
   csh ../advance_model.csh 
###   csh advance_model.csh &
# Pop back up to main directory
   cd ..
   @ element++
end

# Wait for all processes to finish up
wait

# All model runs complete, move the updated file up
set element = 1
while($element <= $num)
   mv tempdir$element/temp_ud assim_model_state_ud$element
   rm -rf tempdir$element
   @ element++
end

# Remove the semaphore file
rm -f async_may_go

# Cleaned up; let the filter know it can proceed
echo All_done:Please_proceed

# Doing recursive call
csh ./async_filter.csh
