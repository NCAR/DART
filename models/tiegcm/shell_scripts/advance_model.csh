#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
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

echo "Starting advance_model for process $process at "`date`

#-------------------------------------------------------------------------
# Block 1: populate a run-time directory with the bits needed to run tiegcm.
#-------------------------------------------------------------------------

# Get unique name for temporary working directory for this process's stuff
set temp_dir = `printf "advance_temp%04d" $process`

# Create a clean temporary directory and go there
\rm -rf  $temp_dir
mkdir -p $temp_dir
cd       $temp_dir

# The following is a list of the REQUIRED input.nml namelist settings
#
# &model_nml         tiegcm_restart_file_name   = 'tiegcm_restart_p.nc',
#                    tiegcm_secondary_file_name = 'tiegcm_s.nc',
#                    tiegcm_namelist_file_name  = 'tiegcm.nml'
# &dart_to_model_nml file_in                    = 'dart_restart',
#                    file_namelist_out          = 'namelist_update',
# &model_to_dart     file_out                   = 'dart_ics'
#
# Ensure that the input.nml has the required value for
# dart_to_model_nml:advance_time_present for this context.

sed -e "/advance_time_present /c\ advance_time_present = .TRUE." \
       ../input.nml >! input.nml || exit 1

# Check to see if you are running async==4 ... an mpirun.lsf situation
set MYSTRING = `grep -A 42 filter_nml input.nml | grep async`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set ASYNC = $MYSTRING[2]

if ($ASYNC == 4) then
   set RUN_CMD = mpirun.lsf
else
   set RUN_CMD = ''
endif

# Loop through each state
set state_copy = 1
set ensemble_member_line = 1
set input_file_line = 2
set output_file_line = 3

while($state_copy <= $num_states)

   set ensemble_member = `head -n $ensemble_member_line ../$control_file | tail -n 1`
   set input_file      = `head -n $input_file_line      ../$control_file | tail -n 1`
   set output_file     = `head -n $output_file_line     ../$control_file | tail -n 1`

   # make sure we have a clean logfile for this entire advance
   set logfile = `printf "log_advance.%04d.txt"     $ensemble_member`

   echo "control_file is ../$control_file"             >! $logfile
   echo "working on ensemble_member $ensemble_member"  >> $logfile
   echo "input_file  is $input_file"                   >> $logfile
   echo "output_file is $output_file"                  >> $logfile
   echo "Starting dart_to_model at "`date`             >> $logfile

   #----------------------------------------------------------------------
   # Block 2: Convert the DART output file to form needed by model.
   # Overwrite the appropriate variables of a TIEGCM netCDF restart file.
   # The DART output file (namelist_update) has the 'advance_to' time 
   # which must be communicated to the model ... through the tiegcm namelist
   #----------------------------------------------------------------------

   set tiesecond   = `printf "tiegcm_s.nc.%04d"         $ensemble_member`
   set tierestart  = `printf "tiegcm_restart_p.nc.%04d" $ensemble_member`
   set tieinp      = `printf "tiegcm.nml.%04d"          $ensemble_member`

   cp -pv ../$input_file dart_restart        >>&   $logfile || exit 2
   cp -pv ../$tiesecond  tiegcm_s.nc         >>&   $logfile || exit 2
   cp -pv ../$tierestart tiegcm_restart_p.nc >>&   $logfile || exit 2
   cp -pv ../$tieinp     tiegcm.nml          >>&   $logfile || exit 2

   ../dart_to_model >>& $logfile || exit 2  # dart_to_model generates namelist_update

   # update tiegcm namelist variables by grabbing the values from namelist_update  
   # and then overwriting whatever is in the tiegcm.nml
   # There is a danger that you match multiple things ... F107 and F107A,
   # SOURCE_START and START, for example ... so try to grep whitespace too ...

   set start_year   = `grep " START_YEAR "   namelist_update`
   set start_day    = `grep " START_DAY "    namelist_update`
   set source_start = `grep " SOURCE_START " namelist_update`
   set start        = `grep " START "        namelist_update`
   set secstart     = `grep " SECSTART "     namelist_update`
   set stop         = `grep " STOP "         namelist_update`
   set secstop      = `grep " SECSTOP "      namelist_update`
   set hist         = `grep " HIST "         namelist_update`
   set sechist      = `grep " SECHIST "      namelist_update`
   set f107         = `grep " F107 "         namelist_update`

   # FIXME TOMOKO ... do we want F107 and F107A to be identical
   #
   # the way to think about the following sed syntax is this:
   # / SearchStringWithWhiteSpaceToMakeUnique  /c\ the_new_contents_of_the_line 
   
   sed -e "/ START_YEAR /c\ ${start_year}" \
       -e "/ START_DAY /c\ ${start_day}" \
       -e "/ SOURCE_START /c\ ${source_start}" \
       -e "/ START /c\ ${start}" \
       -e "/ STOP /c\ ${stop}" \
       -e "/ HIST /c\ ${hist}" \
       -e "/ SECSTART /c\ ${secstart}" \
       -e "/ SECSTOP /c\ ${secstop}" \
       -e "/ SECHIST /c\ ${sechist}" \
       -e "/ F107 /c\ ${f107}" \
       tiegcm.nml >! tiegcm.nml.updated

   if ( -e tiegcm.nml.updated ) then
      echo "tiegcm.nml updated with new start/stop time for ensemble member $ensemble_member"
      mv -v tiegcm.nml tiegcm.nml.original
      mv -v tiegcm.nml.updated tiegcm.nml
   else
      echo "ERROR tiegcm.nml did not update correctly for ensemble member $ensemble_member."
      exit 2
   endif

   #----------------------------------------------------------------------
   # Block 3: Run the model
   #----------------------------------------------------------------------

   echo "ensemble member $ensemble_member : before tiegcm" >>  $logfile
   ncdump -v mtime tiegcm_restart_p.nc                     >>& $logfile
   echo "Starting tiegcm at "`date`                        >>  $logfile

   ${RUN_CMD} ../tiegcm < tiegcm.nml >>& $logfile

   grep -q "NORMAL EXIT" $logfile
   set tiegcm_status = $status

   if ($tiegcm_status != 0) then
      echo "ERROR: tiegcm model advance failed."
      echo "ERROR: check $temp_dir/$logfile"
      exit 3
   endif

   echo "ensemble member $ensemble_member : after tiegcm" >>  $logfile
   ncdump -v mtime tiegcm_restart_p.nc                    >>& $logfile

   #----------------------------------------------------------------------
   # Block 4: Convert the model output to form needed by DART
   # AT this point, the model has updated the information in tiegcm_restart_p.nc
   # We need to get that information back into the DART state vector.
   #
   # The updated information needs to be moved into CENTRALDIR in
   # preparation for the next cycle.
   #----------------------------------------------------------------------

   echo "Starting model_to_dart at "`date`  >>  $logfile
   ../model_to_dart                         >>& $logfile || exit 4

   mv -v dart_ics            ../$output_file || exit 4
   mv -v tiegcm_s.nc         ../$tiesecond   || exit 4
   mv -v tiegcm_restart_p.nc ../$tierestart  || exit 4
   mv -v tiegcm.nml          ../$tieinp      || exit 4

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line      = $input_file_line + 3
   @ output_file_line     = $output_file_line + 3
   echo "Finished model_to_dart at "`date`  >> $logfile
end

# Change back to original directory 
cd ..

# After you are assured this script works as expected, you can actually 
# remove the temporary directory. For now ... leave this commented OUT.
#\rm -rf $temp_dir

# Remove the filter_control file to signal completion
\rm -fv $control_file

echo "Finished advance_model for process $process at "`date`

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

