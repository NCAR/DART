#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
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

set time = `cat async_may_go`
set secs = $time[1]
set days = $time[2]

# First line of filter_control should have number of model states to be integrated
set num = `head -1 filter_control`

# Create a directory for each member to run in for namelists
set element = 1
while($element <= $num)

# Make a temporary directory for this element's run
   mkdir tempdir$element
   cp wrfinput tempdir${element}/wrfinput

# Copy the initial condition file to the temp directory.
   mv assim_model_state_ic$element tempdir${element}/dart_wrf_vector

# Copy the boundary condition file to the temp directory.
   cp /ocotillo1/caya/GEN_TRUTH/wrfbdy_${days}_${secs}_81 tempdir${element}/wrfbdy_d01

# Copy WRF input namelist to the temp directory.
   cp /ocotillo1/caya/GEN_TRUTH/namelist.input_${days}_${secs}_81 tempdir${element}/namelist.input

# Change to temp directory and run integrate
   cd tempdir$element
   ../integrate_wrf >> &out.integration .$element &
   cd ..
   @ element++
end

# Wait for all processes to finish up
wait

# All model runs complete, move the updated file up
set element = 1
while($element <= $num)
   mv tempdir${element}/wrfinput_d01 /ocotillo1/caya/GEN_TRUTH/wrfinput_${days}_${secs}
   mv tempdir${element}/temp_ud assim_model_state_ud$element
   cat tempdir${element}/out.integration >> out.integration$element
   cat tempdir${element}/out_wrf_integration >> out_wrf_integration$element
   cat tempdir${element}/out.dart_to_wrf >> out.dart_to_wrf$element
   cat tempdir${element}/out.wrf_to_dart >> out.wrf_to_dart$element
   rm -rf tempdir$element
   @ element++
end

# Remove the semaphore file
rm -f async_may_go

# Cleaned up; let the filter know it can proceed
echo All_done:Please_proceed

# Doing recursive call
csh ./async_perfect_obs.csh
