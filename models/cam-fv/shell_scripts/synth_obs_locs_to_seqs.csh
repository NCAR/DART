#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# https://www2.cisl.ucar.edu/software/dart/download

# ==============================================================================
# Call create_obs_sequence to convert a text file of obs locations and types
# into a set_def.out file.  
# Call create_fixed_network_seq to create separate files for each time period. 
# This one makes multiple files/day.  There's a single time in each file.
# 2 and 4 are the most common numbers of files/day. 

# It requires only a few seconds to run a single date with ~100,000 observations, 
# in which case it can be run interactively 
# and the batch job directives can be ignored.

# (Copy to and) Submit from a case_dir, where there's 
# + a text file of locations and types for input to create_obs_sequence,
#   (The text file CAN NOT request identity obs unless 
#    suitable cam{input,_phis}.nc files are provided.)
# + in input.nml which was used to build the create_* programs.
#   (minimize printed output by listing only the relevant obs types 
#    in assimilate_these_obs_types, or utilities_nml:write_nml = 'none').
# Caminput.nc and cam_phis.nc will be linked from the executable directory
# where create_* programs were created.
# The resolution is irrelevant because no identity obs will be requested.
# The time period must be within a calendar month.
# Edit the values below to change the dates and intervals.

# ==============================================================================
#PBS  -N synth_obs_locs_to_seqs.csh
#PBS  -A P86850054
#PBS  -q casper
#PBS  -l select=1:ncpus=1:mpiprocs=1
#PBS  -l walltime=00:20:00
#PBS  -o synth_obs_locs_to_seqs_2017-12.out
#PBS  -k eod
#PBS  -j oe 

#----------------------------------------------------------------------------
#BSUB  -J synth_obs_locs_to_seqs_2017-12
#BSUB  -n 1 
#BSUB  -R "span[ptile=1]"
#BSUB  -q shared_node_queue_for_this_setup_script
#BSUB  -P your_account_there
#BSUB  -W 2:00
#BSUB  -u you@email.org
#BSUB  -N  
#BSUB  -a poe 
#BSUB  -o synth_obs_locs_to_seqs_2017-12.out
#BSUB  -e synth_obs_locs_to_seqs_2017-12.out
# ==============================================================================
# standard commands:

set VERBOSE = '-v'
set    MOVE = '/usr/bin/mv'
set    COPY = '/usr/bin/cp --preserve=timestamps'
set    LINK = '/usr/bin/ln -s'
set    LIST = '/usr/bin/ls'
set  REMOVE = '/usr/bin/rm'

set nonomatch       # suppress "rm" warnings if wildcard does not match anything

# ==============================================================================
# User set parameters

@ year      = 2017
@ month     = 12
@ day       = 25
@ day_last  = 31 
@ hour_incr = 6

# case_dir      = Directory from which this script is submitted, 
#                 where job diagnostics will go.
# obs_text_file = Text file of obs locations which will be fed to create_obs_sequence.
# obs_out_dir   = Directory where obs_seq.in files will be created.
# exec_dir      = Directory where create_* programs were built. It is the location 
#                 of the default input files needed by cam-fv's static_init_model.

set case_dir      = `pwd`
set obs_text_file = ${case_dir}/even_create_input_all
set obs_out_dir   = /glade/p/cisl/dares/Observations/Synthetic/UVTRadiosonde_3456
set exec_dir      = /glade/u/home/raeder/DART/reanalysis_git/models/cam-fv/work

# This section could be modified to include running even_sphere.m as part of a single job
# to generate synthetic observations evenly distributed on the sphere.
#     \rm -rf matlab_input.m
# 
#     cat >> matlab_input.m << EndOfInput
# 
#        nprofiles   = 30;
#        levels      = [1000  850  500  300  200  100];
#        T_error_var = [1.44 0.64 0.64 0.81 1.44 0.64];
#        W_error_var = [1.96 2.25 4.41 9.00 7.29 4.41];
#        even_sphere(nprofiles, 'levels', levels, ...
#                   'T_error_var', T_error_var, 'W_error_var', W_error_var)
#        fname = sprintf('even_sphere_%d_profiles',nprofiles);
#        orient landscape
#        print(fname,'-dpdf')
# 
#     EndOfInput
# 
#     matlab -nosplash -nodesktop -r "try; cd $PWD; matlab_input; catch; end; exit";

# ==============================================================================

if (! -d $obs_out_dir) mkdir $obs_out_dir
cd $obs_out_dir

# Link to files needed by CAM's static_init_model subroutine,
# which is used by create_fixed_network_seq.
if (-f input.nml) $MOVE input.nml input.nml.$$
$COPY ${case_dir}/input.nml .

if (-f ${case_dir}/caminput.nc) then
   $COPY ${case_dir}/caminput.nc .
else
   echo "ERROR: failed to find ${case_dir}/caminput.nc"
   exit 10
endif

if (-f ${case_dir}/cam_phis.nc) then
   $COPY ${case_dir}/cam_phis.nc .
else
   echo "ERROR: failed to find ${case_dir}/cam_phis.nc"
   exit 20
endif

# Transform the raw obs location data into a set_def.out file
if (! -f set_def.out) then
   ${exec_dir}/create_obs_sequence < $obs_text_file > /dev/null
   if ($status != 0) then
      echo "create_obs_sequence failed"
      exit 30
   endif
endif

# Tranform the set_def.out into a series of obs_seq.in files
# whose names have the CESM date string format in them,
# for use by cam-fv/shell_scripts/assimilate.csh.
while($day <= $day_last)
   @ hour = 0
   while ($hour < 24)
      # cam-fv's perfect_model.csh wants to see this form.
      set fstring = `printf obs_seq%04d%02d%02d%02d $year $month $day $hour`
      if (-f $fstring) then
         echo "$fstring exists: move it and try again"
         exit 40
      endif

      # Create the file which will be input to create_fixed_network_seq
      echo "set_def.out"               >! create_fixed_input
      echo 1                           >> create_fixed_input
      echo 1                           >> create_fixed_input
      echo $year $month $day $hour 0 0 >> create_fixed_input
      echo '0 0'                       >> create_fixed_input
      echo $fstring                    >> create_fixed_input
   
      # Create an obs_seq.in file
      ${exec_dir}/create_fixed_network_seq < create_fixed_input

      @ hour += $hour_incr
   end

   @ day++
end

rm create_fixed_input

exit 0

