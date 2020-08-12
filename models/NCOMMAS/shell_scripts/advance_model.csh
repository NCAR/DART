#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# This script has 4 logical 'blocks':
# 1) creates a clean, temporary directory in which to run a model instance 
#    and copies the necessary files into the temporary directory
# 2) converts the DART output to input expected by NCOMMAS
# 3) runs NCOMMAS
# 4) converts NCOMMAS output to input expected by DART
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
# run the NCOMMAS model.
#-------------------------------------------------------------------------

# Get unique name for temporary working directory for this process's stuff
set temp_dir = 'advance_temp'${process}
echo "temp_dir is $temp_dir"

# FIXME: have the startup script echo the name of the experiment into a
# known, fixed filename.   for now, assume it is going to be called 'experiment.txt'.
set exp=`cat experiment.txt`

# get the name of the experiment out of the attributes namelist

# If directory doesn't exist, create it.  If it already exists,
# assume what's in it is good.
if ( ! -d $temp_dir) mkdir -p $temp_dir
cd       $temp_dir

# Get the 'changing' namelist files from CENTRALDIR
# Only the namelists in CENTRALDIR have the updated information about
# the state of the model partway through an assimilation experiment.
# NCOMMAS namelist, anything else
foreach FILE ( ../${exp}.3d.run \    
               ../input.nml )
   cp -pv $FILE . || exit 1
end

# Ensure that the input.nml has the required value for
# dart_to_ncommas_nml:advance_time_present for this context.

echo '1'                      >! ex_commands
echo '/dart_to_ncommas_nml'   >> ex_commands
echo '/advance_time_present'  >> ex_commands
echo ':s/\.false\./\.true\./' >> ex_commands
echo ':wq'                    >> ex_commands

( ex input.nml < ex_commands ) >& /dev/null
\rm -f ex_commands

# Link the files that are common to all model advances;
# you cannot put any files that need to be updated in this list.
# (Put them above in the copy list).
# ?empty for now
#foreach FILE ( ../*_contents )
#   if ( ! -e $FILE) ln -sfv $FILE . || exit 1
#end

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
   # Block 2: Convert the DART output file to form needed by NCOMMAS.
   # We are going to take a NCOMMAS netCDF restart file and simply overwrite the
   # appropriate variables. The DART output file also has the 'advance_to'
   # time - which must be communicated to the model ...
   #----------------------------------------------------------------------

   # The EXPECTED input DART 'initial conditions' file name is 'dart_restart'
   # The dart_to_ncommas_nml:advance_time_present = .TRUE. must be set

   # why is this a link and not a move?
   ln -sfv ../$input_file dart_restart || exit 2

   # CENTRALDIR should contain the restart files for all the ensemble members.

   # if this is a problem with a symbolic link, make it a hard link or
   # change it to a move and then move it back when done.
   set RESTARTFILE = `printf %s%03d%s ${exp}.r. $ensemble_member .nc`
   if ( -e ../$RESTARTFILE) then
      ln -sfv ../${RESTARTFILE} ncommas_restart.nc || exit 2
   else
      echo "ERROR: Restart file for ensemble member $ensemble_member is missing."
      echo "Cannot find ../$RESTARTFILE"
      echo "Exiting."
      exit 2
   endif

   # Use the contents of the dart state vector to update the local history file.

   ../dart_to_ncommas || exit 2

   # now give it the right filename which includes the ensemble member number
   mv ncommas_restart.nc ${RESTARTFILE}

   # assume the converter has put two times on one line, the current time
   # and the advance-to time.  These are in seconds from time0.  The name of
   # the file is just 'times'.
   set start_stop = `cat times`


   #----------------------------------------------------------------------
   # Block 3: Run the NCOMMAS model
   #----------------------------------------------------------------------
   # the value of MPI could be inherited - null for now
   set MPI = ''

   if ( "$MPI" == '' ) then
      ${MPI} ../commas ${exp}.3d.run $ensemble_member $start_stop || exit 3
   else
      # ${MPI} ../commas_mpi ${exp}.3d.run $ensemble_member $start_stop || exit 3
   endif

   # get the return status somehow - either from $status, or grepping
   # for a success message in a log file.
   grep "Successful completion of ncommas run" log.*
   set ncommasstatus = $status
   if ( $ncommasstatus != 0 ) then
      echo "ERROR - NCOMMAS ensemble member $ensemble_member did not complete successfully" 
      echo "ERROR - NCOMMAS ensemble member $ensemble_member did not complete successfully" 
      exit 3 
   endif
   
   #----------------------------------------------------------------------
   # Block 4: Convert the model output to form needed by DART
   #----------------------------------------------------------------------

   ls -lrt

   # NCOMMAS updates the restart file.   extract the state vector data
   # and put it into the name requested by dart.

   ln -s $RESTARTFILE restart.nc 
  
   # ncommas_to_dart reads the restart file after the model advance and writes
   # out an updated DART 'initial conditions' file. This initial conditions
   # file contains a header with the valid time of the ensuing model state.
   # The NCOMMAS restart files contain the valid time of the model state.

   ../ncommas_to_dart || exit 4


   # The (new,updated) DART restart file name is called 'dart_ics'
   # Move the updated files back to 'centraldir'
   mv -v dart_ics ../$output_file || exit 4
   rm restart.nc

   # the restart name was a link, so the central version should
   # have been updated.  if not, move it up there.
   #mv -v ${RESTARTFILE} ../${RESTARTFILE} || exit 4

   # bookkeeping

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3
end

# must communicate the time_manager_nml:stop_count 
# cp -pv ncommas_in.DART ../ncommas_in

# Change back to original directory and get rid of temporary directory
cd ..
# \rm -rf $temp_dir

# Remove the filter_control file to signal completion
# Is there a need for any sleeps to avoid trouble on completing moves here?
\rm -rf $control_file

exit 0


