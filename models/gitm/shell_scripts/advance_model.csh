#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
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
#
# Alex: if you read clean.sh, this is exactly the same thing, just adapted into
# DART framework and written in CSH instead of BASH (so beware of @ and setenv and `)

set      process = $1
set   num_states = $2
set control_file = $3

#----------------------------------------------------------------------
# set up commands to supersede any user invocations
#----------------------------------------------------------------------

unalias cd
unalias ls
set  REMOVE = '/bin/rm -rf'
set    COPY = '/bin/cp --preserve=timestamps'
set    MOVE = '/bin/mv -f'
set    LINK = '/bin/ln -sf'

#----------------------------------------------------------------------
# Block 1: copy necessary input files/executables/files common
#          to all model advances to a clean, temporary directory.
#          These will be used by ALL of the ensemble
#          members being advanced by this script.
#----------------------------------------------------------------------

# Loop through each state
set state_copy = 1
set ensemble_member_line = 1
set      input_file_line = 2
set     output_file_line = 3

while($state_copy <= $num_states)

   set ensemble_member = `head -$ensemble_member_line $control_file | tail -1`
   set input_file      = `head -$input_file_line      $control_file | tail -1`
   set output_file     = `head -$output_file_line     $control_file | tail -1`

   # Create a unique temporary working directory for this process's stuff
   # The run-time directory for the entire experiment is called CENTRALDIR;
   # we need to provide a safe haven for each TASK ... in 'temp_dir'.

   set temp_dir = 'advance_temp_e'${ensemble_member}
   cd $temp_dir  || exit 3

   # making the directory is done in script called clean in work dir.
   # CD is done later too

   # Get the program and input.nml DONE IN start_GITM.sh

   #-------------------------------------------------------------------
   # Block 2: copy/convert the DART state vector to something the
   #          model can ingest.
   #
   #          * copy/link ensemble-member-specific files
   #          * convey the advance-to-time to the model
   #          * convert the DART state vector to model format
   #-------------------------------------------------------------------

   ${MOVE} ../$input_file dart_restart || exit 2

   echo 'advance_model.csh is currently in' `pwd`

   ${COPY} restartOUT/header.rst restartOUT.out/.

   ./dart_to_gitm || exit 4

   echo 'advance_model.csh: dart_to_gitm succeeded'

   # remove restart files so that if GITM fails while running, we know for sure
   ${REMOVE} UA/restartOUT/*.rst

   #-------------------------------------------------------------------
   # Block 3: advance GITM
   #-------------------------------------------------------------------

   #create a full UAM.in (the one with start and end times)
   cat DART_GITM_time_control.txt UAM.in >! UAM.in.temporary
   ${MOVE} UAM.in.temporary UAM.in

   echo ' '    >> UAM.in
   echo '#END' >> UAM.in

   echo 'advance_model.csh: UAM.in was successfully created'

   if ( $?PBS_NODEFILE ) then

      # run each ensemble member in the background on separate 
      # processor set specified in hf
      # 'nop' is changed by pbs_file - TJH (fixme?)
      @ nop = 10
      @ b = $nop * $ensemble_member
      head -n $b $PBS_NODEFILE | tail -n $nop >! hf
      echo "member ${ensemble_member} using processors:"
      cat hf

      (mpiexec -hostfile hf -n $nop ./GITM.exe >& ens_log.txt) &

   else

      # Under LSF, each ensemble member gets the full processor set. 
      # No point in running in the background
      mpirun.lsf ./GITM.exe || exit 5
      
   endif

   cd ..

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3

end

sleep 1

# Wait for all background processes to finish.
wait

# Loop through each state again and collect the advanced files when ready

set state_copy = 1
set ensemble_member_line = 1
set      input_file_line = 2
set     output_file_line = 3

while($state_copy <= $num_states)

   set ensemble_member = `head -$ensemble_member_line $control_file | tail -1`
   set input_file      = `head -$input_file_line      $control_file | tail -1`
   set output_file     = `head -$output_file_line     $control_file | tail -1`

   set temp_dir = 'advance_temp_e'${ensemble_member}

   cd $temp_dir || exit 6

   # Check to make sure GITM is done.
   test -e GITM.DONE
   set s = $?
   while ( $s == 1 ) #while it doesn't exist
      echo "G.D doesn't exist in " $temp_dir
      sleep 1
      test -e GITM.DONE
      set s=$?
   end
   ${REMOVE} GITM.DONE

   pwd
   echo 'advance_model.csh: Listing files because sometimes the restart files were different size.'
   echo '                   Check that they have the same size, as the blocks are equal!'
   ls -l UA/restartOUT

   ${REMOVE} UAM.in                   # remove the used full UAM.in
   ${COPY}   UAM.in.truncated UAM.in  # prepare a UAM.in for the next assimilation window
   ./gitm_to_dart || exit 7           # convert GITM restart files to DART file

   echo 'advance_model.csh: ./gitm_to_dart successful'

   # move the updated DART state into place and give it the expected name
   # (something like filter_0003)
   ${MOVE} dart_ics ../$output_file || exit 8

   cd ..

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3

end

# MANDATORY - Remove the control_file to signal completion. If it still
# exists in CENTRALDIR after all the ensemble members have been advanced,
# it means one or more of the advances failed and is an ERROR CONDITION.

${REMOVE} $control_file

exit 0


