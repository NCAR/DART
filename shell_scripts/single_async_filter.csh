#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# Shell script to work with asynchronous filter integration
# This script needs to be piped to the filter program with the
# filter namelist async variable set to .true.
# Make sure that file filter_control is cleared out by a higher level
# script.

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
echo number_of_lines_expected_is_$num


# Would be nice to check to make sure total filter_control
# file length is consistent with this number of model sets
#wc -l filter_control


# Create a directory for each member to run in for namelists
set element = 1
while($element <= $num)
   echo working_on_elemnet_$element
# Make a temporary directory for this element's run
   mkdir tempdir$element
# Copy the initial condition file to the temp directory
   cp -r RESTART tempdir$element
   mv assim_model_state_ic$element tempdir$element/temp_ic
   cp diag_table tempdir$element
   cp input.nml tempdir$element
# Change to temp directory and run integrate
   cd tempdir$element
   ../integrate_model
   cd ..
   @ element++
end


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
async_filter.csh
