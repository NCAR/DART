#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# Shell script to work with asynchronous filter integration
# This script needs to be piped to the filter program with the
# filter namelist async variable set to .true.

set startnode = 1
set endnode = 10

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
set inode = $startnode
set cyclenode = 1
while($element <= $num)

# Make a temporary directory for this element's run
   (rsh -n node$inode "rm -rf /var/tmp/tempdir$element" )
   (rsh -n node$inode "mkdir /var/tmp/tempdir$element" )
   rcp wrfinput node${inode}:/var/tmp/tempdir${element}/wrfinput

# Copy the initial condition file to the temp directory.
   rcp assim_model_state_ic$element node${inode}:/var/tmp/tempdir${element}/dart_wrf_vector
   rm assim_model_state_ic$element

# Copy the boundary condition file to the temp directory.
   rcp /ocotillo1/caya/GEN_INIT_ENS/wrfbdy_${days}_${secs}_$element node${inode}:/var/tmp/tempdir${element}/wrfbdy_d01

# Copy WRF input namelist to the temp directory.
   rcp /ocotillo1/caya/GEN_INIT_ENS/namelist.input_${days}_${secs}_1 node${inode}:/var/tmp/tempdir${element}/namelist.input
   rcp RRTM_DATA node${inode}:/var/tmp/tempdir${element}/
   rcp LANDUSE.TBL node${inode}:/var/tmp/tempdir${element}/

# Copy the required programs
   rcp integrate_wrf2 node${inode}:/var/tmp/tempdir${element}/
   rcp dart_tf_wrf node${inode}:/var/tmp/tempdir${element}/
   rcp wrf.exe node${inode}:/var/tmp/tempdir${element}/

# Run integrate
   set workdir = /var/tmp/tempdir${element}
   (rsh -n node$inode "cd $workdir; integrate_wrf2 >>& out.integration .$element" ) &

   @ element++
   if ($inode == $endnode) then
      set inode = $startnode
      if ($cyclenode == 1) then
         @ cyclenode ++
      else
         set cyclenode = 1
         wait
      endif
   else
      @ inode ++
   endif

   sleep 2

end

# Wait for all processes to finish up
wait

#sleep 10

# All model runs complete, move the updated file up
set element = 1
set inode = $startnode
while($element <= $num)
   rcp node${inode}:/var/tmp/tempdir${element}/temp_ud assim_model_state_ud$element
#   cat tempdir${element}/out.integration >> out.integration$element
#   cat tempdir${element}/out_wrf_integration >> out_wrf_integration$element
#   cat tempdir${element}/out.dart_to_wrf >> out.dart_to_wrf$element
#   cat tempdir${element}/out.wrf_to_dart >> out.wrf_to_dart$element
#   cat tempdir${element}/out.update_wrf_bc >> out.update_wrf_bc$element
   (rsh -n node$inode "rm -rf /var/tmp/tempdir$element" )
   @ element++
   if ($inode == $endnode) then
      set inode = $startnode
   else
      @ inode ++
   endif
end

# Remove the semaphore file
rm -f async_may_go

# Cleaned up; let the filter know it can proceed
echo All_done:Please_proceed

# Doing recursive call
csh ./async_filter_wrf_remote2.csh
