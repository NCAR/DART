#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.
# Can be used as-is with most low-order models and the bgrid model which
# can be advanced using the integrate_model executable.
# 
# Arguments are (created by 'filter' or 'perfect_model_obs' and include):
# 1) the process number of caller,
# 2) the number of ensemble members/state copies belonging to that process, and 
# 3) the name of the control_file for that process.
# 
# If this script finishes and the 'control_file' still exists, it is
# an ERROR CONDITION and means one or more of the ensemble members did
# not advance properly. Despite our best attempts to trap on this
# condition, some MPI installations simply hang, some properly terminate.
#
# This script loops over all the entries in the control_file to advance 
# any/all of the ensemble members.  The number of trips through the 
# loop is the second argument to this script. The control_file contains 
# the information about which ensemble members are to be advanced by THIS TASK.
# Sometimes it may be just one ensemble member, sometimes all of them.
# Read DART/doc/html/filter_async_modes.html and the mpi_intro.html
# for an overview.
#
# This script has 4 logical 'blocks':
# 1) creates a clean, temporary directory in which to run a model instance
#    and copies the necessary files into the temporary directory
# 2) copies/converts the DART state vector to something the model can ingest
# 3) runs the model
# 4) copies/converts the model output to input expected by DART

set      process = $1
set   num_states = $2
set control_file = $3

#----------------------------------------------------------------------
# Block 1: copy necessary input files/executables/files common
#          to all model advances to a clean, temporary directory.
#          These will be used by ALL of the ensemble
#          members being advanced by this script.
#----------------------------------------------------------------------

# Create a unique temporary working directory for this process's stuff
# The run-time directory for the entire experiment is called CENTRALDIR;
# we need to provide a safe haven for each TASK ... in 'temp_dir'.

set temp_dir = 'advance_temp'${process}

# Create a clean temporary directory and go there
\rm -rf  $temp_dir  || exit 1
mkdir -p $temp_dir  || exit 1
cd       $temp_dir  || exit 1

# Get the program and input.nml
ln -s ../integrate_model  .  || exit 1
cp    ../input.nml        .  || exit 1

# Loop through each state
set state_copy = 1
set ensemble_member_line = 1
set      input_file_line = 2
set     output_file_line = 3

while($state_copy <= $num_states)
   
   set ensemble_member = `head -$ensemble_member_line ../$control_file | tail -1`
   set input_file      = `head -$input_file_line      ../$control_file | tail -1`
   set output_file     = `head -$output_file_line     ../$control_file | tail -1`
   
   #-------------------------------------------------------------------
   # Block 2: copy/convert the DART state vector to something the 
   #          model can ingest. In this case, just copy.
   #          In general, there is more to it. 
   #
   #          * copy/link ensemble-member-specific files
   #          * convey the advance-to-time to the model
   #          * convert the DART state vector to model format 
   #-------------------------------------------------------------------

   mv ../$input_file temp_ic || exit 2

   #-------------------------------------------------------------------
   # Block 3: advance the model
   #          In this case, we are saving the run-time messages to
   #          a LOCAL file, which makes debugging easier.
   #          integrate_model is hardcoded to expect input in temp_ic 
   #          and it creates temp_ud as output. 
   #          Your model will likely be different.
   #-------------------------------------------------------------------

   ./integrate_model >! integrate_model_out_temp || exit 3

   # (OPTIONAL) Append the run-time messages to the file in the CENTRALDIR
   #cat integrate_model_out_temp >> ../integrate_model_out_temp${ensemble_member}

   #-------------------------------------------------------------------
   # Block 4: Move the updated state vector back to CENTRALDIR
   #          (temp_ud was created by integrate_model and is in the 
   #          right format already.) In general, you must convert your 
   #          model output to a DART ics file with the proper name.
   #-------------------------------------------------------------------

   mv temp_ud ../$output_file || exit 4

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3
end

# Change back to original directory and get rid of temporary directory.
# If all goes well, there should be no need to keep this directory.
# If you are debugging, you may want to keep this directory. 

cd ..
\rm -rf $temp_dir

# MANDATORY - Remove the control_file to signal completion. If it still
# exists in CENTRALDIR after all the ensemble members have been advanced,
# it means one or more of the advances failed and is an ERROR CONDITION.

\rm -rf $control_file

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

