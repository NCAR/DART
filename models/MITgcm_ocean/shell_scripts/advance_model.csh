#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$

# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.
# Can be used with most low-order models and the bgrid model which
# can be advanced using the integrate_model executable.

# This script copies the necessary files into the temporary directory
# and then executes the fortran program integrate_model.

# Arguments are the 
# 1) process number of caller, 
# 2) the number of state copies belonging to that process, and 
# 3) the name of the filter_control_file for that process
set process = $1
set num_states = $2
set control_file = $3

# Get unique name for temporary working directory for this process's stuff
set temp_dir = 'advance_temp'${process}

# Create a clean temporary directory and go there
\rm -rf  $temp_dir
mkdir -p $temp_dir
cd       $temp_dir

# Get input.nml ... is it used by trans_sv_pv ...?
cp ../input.nml .

# Loop through each state
set state_copy = 1
set ensemble_member_line = 1
set input_file_line = 2
set output_file_line = 3
while($state_copy <= $num_states)
   
   set ensemble_member = `head -$ensemble_member_line ../$control_file | tail -1`
   set input_file      = `head -$input_file_line      ../$control_file | tail -1`
   set output_file     = `head -$output_file_line     ../$control_file | tail -1`
   
   # Get the ics file for this state_copy
   # or ... run trans_sv_pv ... muck about with MIT namelists for model control
   mv ../$input_file assim_model_state_ic

   ../trans_sv_pv

   move the new MIT namelist output ... 

   # Advance the model saving standard out
   # integrate_model is hardcoded to expect input in temp_ic and it creates
   # temp_ud as output.
   ./integrate_model >! integrate_model_out_temp

   echo "some time index relating to the expected [S.xxxxxxxxxx.data]" | ../trans_pv_sv

   # Append the output from the advance to the file in the working directory
   #cat integrate_model_out_temp >> ../integrate_model_out_temp$process

   # Move the updated state vector back up
   # (temp_ud was created by integrate_model.)
   mv assim_model_state_ud ../$output_file

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3
end

# Change back to original directory and get rid of temporary directory
cd ..
\rm -rf $temp_dir

# Remove the filter_control file to signal completion
# Is there a need for any sleeps to avoid trouble on completing moves here?
\rm -rf $control_file

