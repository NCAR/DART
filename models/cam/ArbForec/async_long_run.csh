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
# Shell script to do repeated segment integrations
# BUT, don't know how to do MATLAB from script yet

set output_dir = out_

# Have an overall outer loop
set i = 6
while($i <= 6)

# Run perfect model obs to get truth series
# ! sync_submit.csh runs advance_ens.csh, which hardcodes # of processors
#   to use.  But perfect_model_obs uses only one, so others will be wasted
#   while True_State is generated.
   csh ./sync_submit.csh | ./perfect_model_obs

# Run filter
   csh ./sync_submit.csh | ./filter

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
