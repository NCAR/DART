#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# README: you must denote all observations as 'evaluate only' and turn OFF any/all inflation.
# To get identical values (to the existing ascii restarts), I had to change the 
# trunk/models/brid_solo/model_mod.f90 to create the diagnostic variables as nf90_double values 
# (they were nf90_real).  
#
# 1) takes the first timestep from True_State.nc and creates a single netCDF restart file.
# 2) takes the first timestep from Prior_Diag.nc and creates a (multi-member) single netCDF restart file.
#    There are no posterior inflation values in this case.
#
# only the non-obvious NCO tasks are explained here.
# ncks -O -d time,0 -d copy,82 Prior_Diag.nc bob.nc         rip out the right copy/copies
# ncwa -O -a copy,0            bob.nc prior_inflation.nc    removes 'copy' dimension
# ncks -O -x -v copy,inputnml  prior_inflation.nc bob.nc    removes copy, inputnml variables
#
# Since version control systems cannot store incremental differences, we actually convert
# the netCDF to it's ASCII representation (a .cdl file). It's tolerable for smallish files.
# The prior inflation standard deviation had to be hand-edited from 0 to a more meaningful 0.6
# The prior inflation value was already a useful 1.0
#
# This script will generate perfect_input.nc and filter_input.nc
# To create the matching .cdl files is easy:
#   ncdump perfect_input.nc > perfect_input.cdl
#   ncdump  filter_input.nc >  filter_input.cdl
# and to create a netCDF file from a .cdl file is also easy:
#   ncgen -o perfect_input.nc perfect_input.cdl
#   ncgen -o  filter_input.nc  filter_input.cdl
#
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# bgrid_solo - True_State.nc -> perfect_input.nc

ncks -O -d time,0 -d copy,0 True_State.nc perfect_input.nc
ncrename -O -d copy,member -v copy,member perfect_input.nc temp.$$.nc
ncks -O -x -v inputnml                    temp.$$.nc perfect_input.nc

# remove extraneous history, modifies a couple
ncatted  -O -h \
         -a assim_model_source,global,d,, \
         -a assim_model_revision,global,d,, \
         -a assim_model_revdate,global,d,, \
         -a netcdf_version,global,d,, \
         -a creation_date,global,d,, \
         -a model_source,global,d,, \
         -a model_revision,global,d,, \
         -a model_revdate,global,d,, \
         -a NCO,global,d,, \
         -a long_name,time,m,c,'valid time of the model state' \
         -a title,global,m,c,'a spun-up model state' \
         -a history,global,m,c,'same values as in perfect_ics r1025 (circa Dec 2004)' \
            perfect_input.nc

#-------------------------------------------------------------------------------
# bgrid_solo - Prior_Diag.nc -> filter_input.nc
#
# Like I said, there is no posterior inflation in this test case.
# The prior inflation mean was 1.0

# Generate the prior inflation mean netCDF file

ncks -O -d time,0 -d copy,82 Prior_Diag.nc temp.$$.nc
ncwa -O -a copy,0            temp.$$.nc prior_inf_mean.$$.nc
ncks -O -x -v copy,inputnml  prior_inf_mean.$$.nc temp.$$.nc
ncrename -O \
         -v ps,ps_priorinf_mean \
         -v t,t_priorinf_mean \
         -v u,u_priorinf_mean \
         -v v,v_priorinf_mean \
         temp.$$.nc prior_inf_mean.$$.nc

ncatted -O -h \
        -a  units,,d,, \
        -a  cell_methods,,d,, \
        -a  units_long_name,,d,, \
        -a  long_name,ps_priorinf_mean,m,c,'prior inflation value for surface pressure' \
        -a  long_name,t_priorinf_mean,m,c,'prior inflation value for temperature' \
        -a  long_name,u_priorinf_mean,m,c,'prior inflation value for u wind component' \
        -a  long_name,v_priorinf_mean,m,c,'prior inflation value for v wind component' \
             prior_inf_mean.$$.nc

# Generate the prior inflation sd netCDF file

ncks -O -d time,0 -d copy,83 Prior_Diag.nc temp.$$.nc
ncwa -O -a copy,0            temp.$$.nc prior_inf_sd.$$.nc
ncks -O -x -v copy,inputnml  prior_inf_sd.$$.nc temp.$$.nc

ncrename -O \
         -v ps,ps_priorinf_sd \
         -v t,t_priorinf_sd \
         -v u,u_priorinf_sd \
         -v v,v_priorinf_sd \
         temp.$$.nc prior_inf_sd.$$.nc

ncatted -O -h \
        -a  units,,d,, \
        -a  cell_methods,,d,, \
        -a  units_long_name,,d,, \
        -a  long_name,ps_priorinf_sd,m,c,'prior inflation standard deviation for surface pressure' \
        -a  long_name,t_priorinf_sd,m,c,'prior inflation standard deviation for temperature' \
        -a  long_name,u_priorinf_sd,m,c,'prior inflation standard deviation for u wind component' \
        -a  long_name,v_priorinf_sd,m,c,'prior inflation standard deviation for v wind component' \
             prior_inf_sd.$$.nc

# Generate the input un-inflated ensemble member netCDF file

ncks     -O -d time,0 -d copy,2,81 Prior_Diag.nc members.$$.nc
ncrename -O -d copy,member -v copy,member members.$$.nc temp.$$.nc
ncks     -O -x -v inputnml temp.$$.nc members.$$.nc

# Generate a scalar for the advance_to_time

ncks -O -d time,0 -d copy,0 -v time Prior_Diag.nc adv_to_time.$$.nc
ncrename -O -d time,advance_to_time -v time,advance_to_time adv_to_time.$$.nc temp.$$.nc
ncwa -O -a advance_to_time,0 temp.$$.nc adv_to_time.$$.nc
ncatted -O -h \
        -a  cell_methods,,d,, \
        -a  long_name,advance_to_time,m,c,'desired time at end of the next model advance' adv_to_time.$$.nc

# combine everything into one netCDF file

ncks -A -v ps_priorinf_mean,t_priorinf_mean,u_priorinf_mean,v_priorinf_mean \
           prior_inf_mean.$$.nc members.$$.nc
ncks -A -v ps_priorinf_sd,t_priorinf_sd,u_priorinf_sd,v_priorinf_sd \
           prior_inf_sd.$$.nc members.$$.nc
ncks -A -v advance_to_time \
           adv_to_time.$$.nc members.$$.nc
ncatted  -O -h \
         -a  long_name,time,m,c,'valid time of the model state' \
         -a  title,global,m,c,'an ensemble of spun-up model states' \
         members.$$.nc

# remove extraneous history
ncatted -O -h -a assim_model_source,global,d,, \
              -a assim_model_revision,global,d,, \
              -a assim_model_revdate,global,d,, \
              -a netcdf_version,global,d,, \
              -a creation_date,global,d,, \
              -a model_source,global,d,, \
              -a model_revision,global,d,, \
              -a model_revdate,global,d,, \
              -a NCO,global,d,, \
              -a history_of_appended_files,global,d,, \
              -a history,global,m,c,'same values as in filter_ics r1025 (circa Dec 2004)' \
                 members.$$.nc

mv members.$$.nc filter_input.nc

\rm temp.$$.nc adv_to_time.$$.nc prior_inf_mean.$$.nc prior_inf_sd.$$.nc

exit 0

