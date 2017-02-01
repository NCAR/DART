#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# Shell script to do repeated segment integrations
# BUT, don't know how to do MATLAB from script yet

set output_dir = out_

# Have an overall outer loop
set i = 1
while($i <= 2)

# Run perfect model obs to get truth series
   ./perfect_model_obs

# Run filter
   ./filter

# Not yet; Run diagnostics to generate some standard figures

# Move the netcdf files to an output directory
   mkdir $output_dir$i
   mv True_State.nc Prior_Diag.nc Posterior_Diag.nc $output_dir$i

# Copy the perfect model and filter restarts to start files
   cp filter_restart filter_ics
   cp perfect_restart perfect_ics

# Move along to next iteration
   @ i++

end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

