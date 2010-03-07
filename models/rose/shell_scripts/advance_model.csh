#!/bin/csh
#
# DART software - Copyright © 2004 - 2010 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# This script has 4 logical 'blocks':
# 1) creates a clean, temporary directory in which to run a model instance 
#    and copies the necessary files into the temporary directory
# 2) converts the DART output to input expected by the model
# 3) runs the model
# 4) converts the model output to input expected by DART
#
# The error code from the script reflects which block it failed.
#
# Arguments are the 
# 1) process number of caller, 
# 2) the number of state copies belonging to that process, and 
# 3) the name of the filter_control_file for that process

set process = $1
set num_states = $2
set control_file = $3

#-------------------------------------------------------------------------
# Block 1: populate a run-time directory with the bits needed to run rose.
#-------------------------------------------------------------------------

# Get unique name for temporary working directory for this process's stuff
set temp_dir = 'advance_temp'${process}

# Create a clean temporary directory and go there
\rm -rf  $temp_dir
mkdir -p $temp_dir
cd       $temp_dir

# Get the data files needed to run rose. One directory up is 'CENTRALDIR'

cp ../input.nml .
cp ../rose.nml rose.nml_default

# Loop through each state
set state_copy = 1
set ensemble_member_line = 1
set input_file_line = 2
set output_file_line = 3

while($state_copy <= $num_states)

   set ensemble_member = `head -$ensemble_member_line ../$control_file | tail -1`
   set input_file      = `head -$input_file_line      ../$control_file | tail -1`
   set output_file     = `head -$output_file_line     ../$control_file | tail -1`

   #----------------------------------------------------------------------
   # Block 2: Convert the DART output file to form needed by model.
   # We are going to take a POP netCDF restart file and simply overwrite the
   # appropriate variables. The DART output file also has the 'advance_to'
   # time - which must be communicated to the model ...
   #----------------------------------------------------------------------

   # The EXPECTED input DART 'initial conditions' file name is 'temp_ic'
   # The dart_to_pop_nml:advance_time_present = .TRUE. must be set

   ln -sfv ../$input_file temp_ic || exit 2
   cp -p   ../rose_restart.nc  .  || exit 2

   ../dart_to_model || exit 2

   # Convey the new 'advance_to' time to rose via the namelist
   # trans_time creates a teeny file called 'times' that contains
   # the 'advance_to' time from DART in the rose format
   # The program nmlbld_rose takes the rose template namelist and
   # inserts the proper time and ensemble member bits.

   ../trans_time

   echo `cat times`        >! namelist.in
   echo $ensemble_member   >> namelist.in
#  echo `cat a_tunes`      >> namelist.in
#  echo `cat p_tunes`      >> namelist.in

   ../nmlbld_rose  < namelist.in
   echo "advance_model: after nmlbld_rose"

   ls -lrt   

   #----------------------------------------------------------------------
   # Block 3: Run the model
   #----------------------------------------------------------------------

   ../rose |& tee rose_out_$ensemble_member

   ls -lrt

   #----------------------------------------------------------------------
   # Block 4: Convert the model output to form needed by DART
   #----------------------------------------------------------------------

   ../model_to_dart

   mv temp_ud ../$output_file || exit 4

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3
end

# Change back to original directory 
cd ..

# After you are assured this script works as expected, you can actually 
# remove the temporary directory. For now ... leave this commented OUT.
#\rm -rf $temp_dir

# Remove the filter_control file to signal completion
\rm -rf $control_file

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

