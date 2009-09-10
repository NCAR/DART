#!/bin/tcsh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2009, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$
#
# This script has 4 logical 'blocks':
# 1) creates a clean, temporary directory in which to run a model instance 
#    and copies the necessary files into the temporary directory
# 2) converts the DART output to input expected by the ocean model
# 3) runs the ocean model
# 4) converts the ocean model output to input expected by DART
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
# Block 1: populate a run-time directory with the bits needed to 
# run the ocean model.
#-------------------------------------------------------------------------

# Get unique name for temporary working directory for this process's stuff
set temp_dir = 'advance_temp'${process}
echo "temp_dir is $temp_dir"

# Create a clean temporary directory and go there
\rm -rf  $temp_dir
mkdir -p $temp_dir
cd       $temp_dir

# Get the 'changing' namelist files from CENTRALDIR
# Only the namelists in CENTRALDIR have the updated information about
# the state of the model partway through an assimilation experiment.
foreach FILE ( ../pop_in.part1 \
               ../pop_in.part2 \
               ../input.nml )
   cp -pv $FILE . || exit 1
end

# copy the files used by 
foreach FILE ( ../horiz_grid.gx3v5.* \
               ../topography.gx3v5.* \
               ../vert_grid.gx3v5    \
               ../*_contents )
   ln -sf $FILE . || exit 1
end

echo 'listing now that the table has been set ...'
ls -l

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
   # Block 2: Convert the DART output file to form needed by ocean model.
   # We are going to take a POP netCDF restart file and simply overwrite the
   # appropriate variables. The DART output file also has the 'advance_to'
   # time - which must be communicated to the model ...
   #----------------------------------------------------------------------

   # The EXPECTED input DART 'initial conditions' file name is 'dart.ic'
   # The dart_to_pop_nml:advance_time_present = .TRUE. must be set

   ln -sfv ../$input_file dart.ic || exit 2

   # CENTRALDIR will always contain a pointer file containing the name
   # of the most recent POP restart file for this ensemble member.
   # Locally, the POP restart file name is always pop.r.nc

   if ( -e ../rpointer.ocn.${ensemble_member}.restart ) then
      # dereference the pointer file
      set RESTARTFILE = `head -1 ../rpointer.ocn.${ensemble_member}.restart`
      cp -pv ../${RESTARTFILE} pop.r.nc || exit 2
   
   else
      echo "Pointer file for ensemble member $ensemble_member is missing."
      echo "Looking for "`pwd`" ../rpointer.ocn.${ensemble_member}.restart"
      echo "Exiting ... (pointer file not found in CENTRALDIR)"
      exit 2
   endif

   # create a pop_in to satisfy dart_pop_mod:initialize_module()

   cat pop_in.part1 pop_in.part2 >! pop_in

   ../dart_to_pop || exit 2

   # Convey the new POP 'advance_to' time to POP via the namelist
   cat pop_in.DART pop_in.part2 >! pop_in

   # POP needs a pointer file containing the restart filename
   echo "pop.r.nc"       >! rpointer.ocn.restart
   echo "RESTART_FMT=nc" >> rpointer.ocn.restart

   #----------------------------------------------------------------------
   # Block 3: Run the ocean model
   # The CCSM version has a pointer file that contains the name of the
   # last restart. The LANL version has no such mechanism, but the 
   # filename can be predicted from the pop_in namelist information.
   #----------------------------------------------------------------------

   mpirun.lsf ../pop || exit 3
   
   #----------------------------------------------------------------------
   # Block 4: Convert the ocean model output to form needed by DART
   #----------------------------------------------------------------------

   ls -lrt

   # POP makes a new restart file and updates the pointer file
   set RESTARTFILE = `head -1 rpointer.ocn.restart`
   echo "POP member $ensemble_member made restart file $RESTARTFILE"
   ln -svf ${RESTARTFILE} pop.r.nc || exit 2
   
   # pop_to_dart reads the restart file after the model advance and writes
   # out an updated DART 'initial conditions' file. This initial conditions
   # file contains a header with the valid time of the ensuing model state.
   # The POP restart files contain the valid time of the model state.

   ../pop_to_dart || exit 4

   # The (new,updated) DART restart file name is called 'dart.ud'
   # Move the updated files back to 'centraldir'
   mv -v dart.ud ../$output_file || exit 4
   mv -v rpointer.ocn.restart ../rpointer.ocn.${ensemble_member}.restart || exit 4
   mv -v ${RESTARTFILE} ../${RESTARTFILE} || exit 4

   # bookkeeping

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3
end

# must communicate the time_manager_nml:stop_count 
# cp -pv pop_in.DART ../pop_in

# Change back to original directory and get rid of temporary directory
cd ..
# \rm -rf $temp_dir

# Remove the filter_control file to signal completion
# Is there a need for any sleeps to avoid trouble on completing moves here?
\rm -rf $control_file

