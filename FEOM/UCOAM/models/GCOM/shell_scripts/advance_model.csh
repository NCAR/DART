#!/bin/tcsh
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

# echo "My process ID is $process"
# echo "I am responsible for advancing $num_states states."
# echo "I am going to read my filenames from $control_file"

setenv   MOVE 'mv -v'
setenv   COPY 'cp -v --preserve=timestamps'
setenv   LINK 'ln -vs'
setenv REMOVE 'rm -fr'

#-------------------------------------------------------------------------
# Block 1: populate a run-time directory with the bits needed to
# run the ocean model.
#-------------------------------------------------------------------------

# Create a unique name for temporary working directory for this process.
# Create a clean temporary directory and go there.

set temp_dir = `printf "temp_advance_%04d" $process`

${REMOVE}  $temp_dir
mkdir -p   $temp_dir
cd         $temp_dir

# Copy all the static data from the previous directory "CENTRALDIR"
# into this directory where the model advance takes place.

${LINK} ../gcom_geometry.nc      .  || exit 1
${LINK} ../Grid.dat              .  || exit 1
${LINK} ../Gridll.dat              .  || exit 1
${LINK} ../ProbSize.dat          .  || exit 1

# Ensure that the input.nml has the required value for
# dart_to_gcom_nml:advance_time_present for this context.
# the way to think about the following sed syntax is this:
# / SearchStringWithWhiteSpaceToMakeUnique  /c\ the_new_contents_of_the_line

sed -e "/ advance_time_present /c\ advance_time_present = .true." \
       ../input.nml >! input.nml

if (-z input.nml) then
   echo "the advance_time_present sed failed ..."
   exit 1
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

   echo "control_file is ../${control_file}"           >! $logfile
   echo "working on ensemble_member $ensemble_member"  >> $logfile
   echo "input_file  is ${input_file}"                 >> $logfile
   echo "output_file is ${output_file}"                >> $logfile
   echo "starting advance_model.csh at "`date`         >> $logfile

   #----------------------------------------------------------------------
   # Block 2: Convert the DART output file to form needed by ocean model.
   # Copy the 'most current' gcom_restart file for this ensemble member
   # from CENTRALDIR to the local directory and simply overwrite the
   # appropriate variables with the state from the DART assimilation.
   # dart_to_gcom also creates dart_gcom_timeinfo.txt, which contains
   # the current time of the model as well as the desired stopping time.
   # It has this information as several character strings that are then
   # extracted and used to modify namelists or construct filenames.
   #----------------------------------------------------------------------

   # The EXPECTED DART dart_to_gcom_input_file is 'dart_restart'
   # The dart_to_gcom_nml:advance_time_present = .TRUE. must be set
   # CENTRALDIR will always contain the gcom_restart_nnnn.nc file.
   # gcom_restart_nnnn.nc can be considered a 'pointer' file.
   # The latest,greatest GCOM state will always be gcom_restart_nnnn.nc

   set GCOM_POINTER = `printf gcom_restart_%04d.nc ${ensemble_member}`

   echo "GCOM RESTARTFILE is [${GCOM_POINTER}]"    >> $logfile || exit 2

   if ( -e dart_restart    ) ${REMOVE} dart_restart
   if ( -e gcom_restart.nc ) ${REMOVE} gcom_restart.nc

   ${LINK} ../${input_file}     dart_restart      >> $logfile || exit 2
   ${COPY} ../${GCOM_POINTER}   gcom_restart.nc   >> $logfile || exit 2

   ../dart_to_gcom                     >> $logfile || exit 2

   ls -lR                              >> $logfile
   echo "finished dart_to_gcom "`date` >> $logfile

   # Convey the new start/stop times to gcom via param.dat

   set newstartdate = `grep Start_Time      dart_gcom_timeinfo.txt`
   set  sec2advance = `grep Stop_Time_sec   dart_gcom_timeinfo.txt`
   set   currenttag = `grep currenttimetag  dart_gcom_timeinfo.txt | sed -e "s/[a-z,= ]//g"`
   set  forecasttag = `grep forecasttimetag dart_gcom_timeinfo.txt | sed -e "s/[a-z,= ]//g"`

   echo "gcomstartdate is ${newstartdate}" >> $logfile
   echo "sec2advance   is ${sec2advance}"  >> $logfile

   sed -e "/ Start_Time /c\ ${newstartdate}" \
       -e "/ Stop_Time_sec /c\ ${sec2advance}" \
          ../param.dat >! param.dat || exit 2

   if ( -e param.dat ) then
      echo "param.dat updated with of ${sec2advance} for member $ensemble_member" >> $logfile
   else
      echo "ERROR param.dat did not update for ensemble member $ensemble_member" >> $logfile
      echo "ERROR Stop Time in sec is ${sec2advance}" >> $logfile
      exit 1
   endif

   echo "before running gcom "`date` >> $logfile

   #----------------------------------------------------------------------
   # Block 3: Run the ocean model ... produces gcom_output.nc
   #          If the model advanced successfully, rename the output file
   #          so it contains the ensemble member and current model time.
   #          Must push this output back to CENTRALDIR and update the link
   #          to gcom_restart_nnnn.nc so that the next time gcom_restart.nc
   #          is copied to TEMPDIR, we get the most current state.
   #----------------------------------------------------------------------

   ../gcom.serial.exe >> $logfile || exit 3

   echo "after running gcom "`date`       >> $logfile
   echo "Should now be at ${forecasttag}" >> $logfile

   grep "UCOAM Finished successfully." $logfile
   set gcomstatus = $status
   if ( $gcomstatus != 0 ) then
      echo "ERROR - gcom ensemble member $ensemble_member did not complete successfully"
      echo "ERROR - gcom ensemble member $ensemble_member did not complete successfully"
      echo "ERROR - gcom ensemble member $ensemble_member did not complete successfully" >> $logfile
      echo "ERROR - gcom ensemble member $ensemble_member did not complete successfully" >> $logfile
      exit 3
   endif

   if ( -e gcom_output.nc ) then
      # gcom_to_dart expects a filename of gcom_restart.nc 
      ${MOVE} gcom_output.nc gcom_restart.nc || exit 3
   else
      echo "ERROR - gcom ensemble member $ensemble_member did not create gcom_output.nc"
      echo "ERROR - gcom ensemble member $ensemble_member did not create gcom_output.nc"
      echo "ERROR - gcom ensemble member $ensemble_member did not create gcom_output.nc" >> $logfile
      echo "ERROR - gcom ensemble member $ensemble_member did not create gcom_output.nc" >> $logfile
      exit 3
   endif

   #----------------------------------------------------------------------
   # Block 4: Convert the ocean model output to form needed by DART
   #----------------------------------------------------------------------

   # gcom_to_dart reads the GCOM file after the model advance and writes
   # out an updated DART 'initial conditions' file named 'dart_ics'.
   # Move the updated files back to 'centraldir' and update the pointer file
   # to point to the newest GCOM restart file for the next advance.

   ../gcom_to_dart                    >> $logfile || exit 4
   ${MOVE} dart_ics ../${output_file} >> $logfile || exit 4

   # This just renames the file to contain the instance number and valid time
   set NEWNAME = `printf gcom_restart_%04d.${forecasttag}.nc ${ensemble_member}`

   ${MOVE} gcom_restart.nc ../${NEWNAME} >> $logfile || exit 4

   ${REMOVE} ../${GCOM_POINTER}
   (cd ..; ${LINK} ${NEWNAME} ${GCOM_POINTER}; cd -) >> $logfile || exit 4

   # bookkeeping

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3
end

# Change back to original directory and get rid of temporary directory
cd ..
# \rm -rf ${temp_dir}

# Remove the filter_control file to signal completion
\rm -rfv ${control_file}

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

