#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# Shell script to do repeated segment integrations
#
# $Id$

set output_dir = out_

# Have an overall outer loop
set i = 3
while($i <= 4)

# Run perfect model obs to get truth series
   csh cam_async.csh | ./cam_perfect_model_obs

# Run filter
   csh cam_async.csh | ./cam_filter

# Move the netcdf files to an output directory
   mkdir $output_dir$i
   mv True_State.nc Prior_Diag.nc Posterior_Diag.nc $output_dir$i

# Copy the perfect model and filter restarts to start files
   cp filter_restart filter_ics
   cp cam_perfect_restart cam_perfect_ics

# Move along to next iteration
   @ i++

end
