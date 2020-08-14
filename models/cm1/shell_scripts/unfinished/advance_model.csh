#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#-----------------------------------------------------------------------------
#
# Script for use in assimilation applications
# where the model advance is executed as a separate process.
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
# 1) determine how many ensemble members (num_states) this process 
#    will need to advance.
# 2) convey the forecast length to the model
# 3) run the model (make the forecast)
# 4) determine if there are more ensemble members to advance

set      process = $1
set   num_states = $2
set control_file = $3

echo "advance_model: process $process"
echo "advance_model: number of states to advance $num_states"
echo "advance_model: control file name is $control_file"

set CENTRALDIR = `pwd`

# Loop through each model state

set ensemble_member_line = 1
set      input_file_line = 2
set     output_file_line = 3
set    current_time_line = 4
set advance_to_time_line = 5

set  state_copy = 1

while($state_copy <= $num_states)

   #-------------------------------------------------------------------
   # Block 1: Parse the information from the control file
   #
   # NOTE: the input_file and output_file are not used by this script
   #       because of the implementation of the input_filelist.txt constuct.
   #       They are here merely to maintain consistency with old versions
   #       of advance_model.csh scripts. Could easily be removed.
   #-------------------------------------------------------------------

   set ensemble_member = `head -n $ensemble_member_line  $control_file | tail -n 1`
   set input_file      = `head -n $input_file_line       $control_file | tail -n 1`
   set output_file     = `head -n $output_file_line      $control_file | tail -n 1`
   set current_time    = `head -n $current_time_line     $control_file | tail -n 1`
   set advance_to_time = `head -n $advance_to_time_line  $control_file | tail -n 1`

   set ensemble_member = `printf %03d ${ensemble_member}`

   # The run-time directory for the entire experiment is called CENTRALDIR;
   # we need to provide a unique directory for each model advance.

   set temp_dir = `printf dir_model_%03d ${ensemble_member}`

   cd $temp_dir || exit 1

   #-------------------------------------------------------------------
   # Block 2:  Convey the forecast length to the model.
   #
   # Create two time strings that we can difference to get the number
   # of seconds to advance which must be conveyed to CM1. We need to
   # parse the times from the control_file (whose format is determined
   # by DART in consideration of all the models) to a time useful for CM1
   #-------------------------------------------------------------------

   # The dart function 'advance_time' needs an input.nml

   ln -svf  ../input.nml  .  || exit 1

   set bob=`echo $current_time`
   set year   = $bob[2]
   set month  = $bob[3]
   set day    = $bob[4]
   set hour   = $bob[5]
   set minute = $bob[6]
   set second = $bob[7]
   set cur_datestr = `printf %04d%02d%02d%02d%02d%02d $year $month $day $hour $minute $second`
   set current_d=(`echo ${cur_datestr} 0 -g | ../advance_time`)

   #echo "We think the current  date string is: ${cur_datestr}"
   #echo "Current  days and seconds: $current_d"
   #echo "Started with: ${current_time}"

   set bob=`echo $advance_to_time`
   set year   = $bob[2]
   set month  = $bob[3]
   set day    = $bob[4]
   set hour   = $bob[5]
   set minute = $bob[6]
   set second = $bob[7]
   set fcst_datestr = `printf %04d%02d%02d%02d%02d%02d $year $month $day $hour $minute $second`
   set forecast_d =(`echo ${fcst_datestr} 0 -g | ../advance_time`)

   #echo "We think the forecast date string is: ${fcst_datestr}"
   #echo "Forecast days and seconds: $forecast_d"
   #echo "Started with: ${advance_to_time}"

   @ totaldays = ( $forecast_d[1] - $current_d[1] )
   @ totalsecs = ( $forecast_d[2] - $current_d[2] )

   @ run_time = $totaldays * 86400 + $totalsecs

   #echo "DEBUG:Run time: $run_time"

   # Now that we know the forecast length, must convey that to CM1 by
   # replacing the bogus string with the real forecast length.
   # Both run_time and rstfrq (restart frequency) are changed.
   # FIXME: check the impact of changing the restart frequency
   # has on whole thing by varying the forecast lengths.
   #
   # NOTE: irst == 1, always. Indicates using an existing state.
   # By forcing rstnum == 1, we can always link the most current
   # file to be cm1out_rst_000001.nc

   rm -f namelist.input
   sed -e "s/CM1_FORECAST_LENGTH/${run_time}/" namelist.input.template > namelist.input

   grep CM1_FORECAST_LENGTH namelist.input
   if ($status == 0) then
      echo "The CM1 namelist file 'namelist.input' did not get the new run_time."
      echo "Aborting."
      exit 2
   endif

   #-------------------------------------------------------------------
   # Block 3: Advance the model.
   #
   # Saving the run-time messages from each file to a unique log file.
   # This is intended to make debugging easier.
   #
   # CM1 always expects a restart file called "cm1out_rst_000001.nc".
   #
   # DART is going to save the most current state with a date/time tag 
   #    and must also link that most current state to the static 
   #    "cm1out_rst_000001.nc" name.
   #-------------------------------------------------------------------

   rm -rf cm1log.${fcst_datestr}.txt

   ../cm1.exe |& tee cm1log.${fcst_datestr}.txt
   grep 'Program terminated normally' cm1log.${fcst_datestr}.txt
   if ($status != 0) then
      echo "ERROR: cm1 model advance failed."
      echo "ERROR: check $temp_dir/cm1log.${fcst_datestr}.txt"
      exit 3
   else
      echo "CM1 should now be at $fcst_datestr"
   endif

   # Trying to grab the most recent restart file
   set outfile = `ls -1 cm1out_rst_??????.nc | tail -n 1`
   mv -v ${outfile} cm1out.$ensemble_member.$fcst_datestr.nc

   # Update filename referenced by CENTRALDIR/input_filelist.txt
   # to point to the most recent restart file.
   ln -svf cm1out.$ensemble_member.$fcst_datestr.nc cm1out_rst_000001.nc

   #-------------------------------------------------------------------
   # Block 4:
   # Update the location in the control file with the information
   # for the next ensemble member
   #-------------------------------------------------------------------

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 5
   @ input_file_line = $input_file_line + 5
   @ output_file_line = $output_file_line + 5
   @ current_time_line = $current_time_line + 5
   @ advance_to_time_line = $advance_to_time_line + 5

   # Change back to original directory
   cd $CENTRALDIR
end

# MANDATORY - Remove the control_file to signal completion. If it still
# exists in CENTRALDIR after all the ensemble members have been advanced,
# it means one or more of the advances failed and is an ERROR CONDITION.

\rm -rf $control_file

exit 0


