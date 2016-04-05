#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Script to advance one ensemble member one filter "time step"
# when the model advance is executed as a separate process.
# Called by the filter executable (for async=2 or 4)
# Calls run-cam.csh, the CAM execution script.
# Calls translation program to translate time and model state.
# Runs on one of the compute nodes allotted to the filter executable.
#
# Arguments are 
# arg#1  the process number of caller
# arg#2  the number of state copies belonging to that process
# arg#3  the name of the filter_control_file for that process

set process = $1
set num_states = $2
set control_file = $3

# Get unique name for temporary working directory for this process's stuff
set temp_dir = 'advance_temp'${process}

# args to previous version of this script
set     myname = $0
set CENTRALDIR = `pwd`

# Create a clean temporary directory and go there
\rm -rf  $temp_dir
mkdir -p $temp_dir
cd       $temp_dir

# People have the craziest aliases. These prevent the obsessive-compulsive
# from causing themselves no end of angst.
if ( ! $?REMOVE ) then
  set REMOVE = 'rm -rf'
endif
if ( ! $?COPY ) then
  set COPY = 'cp -fp'
endif
if ( ! $?MOVE ) then
  set MOVE = 'mv -f'
endif
if ( ! $?LINK ) then
  set LINK = 'ln -fs'
endif

echo "advance_model.csh args = $1 $2 $3"                    >  cam_out_temp
echo "CENTRALDIR is ${CENTRALDIR}"                          >> cam_out_temp
echo "temp_dir is $temp_dir"                                >> cam_out_temp

# Get information about this experiment from file "casemodel", 
# created by the main controlling script (job.csh)

# set $case = the case and 
# $model = the directory name (in the CAM source tree) 
# where CAM executable will be found. 
# set locations of the CAM and CLM input files
set case        = `head -1 ${CENTRALDIR}/casemodel | tail -1`
set model       = `head -2 ${CENTRALDIR}/casemodel | tail -1`
set cam_init    = `head -3 ${CENTRALDIR}/casemodel | tail -1`
set clm_init    = `head -4 ${CENTRALDIR}/casemodel | tail -1`
set list        = `head -5 ${CENTRALDIR}/casemodel | tail -1`
set ice_init    = $list[1]

# output diagnostic information to the same file as the CAM list-directed output
echo "case $case model $model"    >> cam_out_temp
echo "cam init is $cam_init"      >> cam_out_temp
echo "clm init is $clm_init"      >> cam_out_temp
echo "ice init is $ice_init"      >> cam_out_temp

# Loop through each ensemble this task is responsible for advancing.
set ensemble_number_line = 1
set input_file_line      = 2
set output_file_line     = 3
set state_copy = 1
while($state_copy <= $num_states)

   # loop through the control file, extracting lines in groups of 3.
   set ensemble_number = `head -$ensemble_number_line ../$control_file | tail -1`
   set input_file      = `head -$input_file_line      ../$control_file | tail -1`
   set output_file     = `head -$output_file_line     ../$control_file | tail -1`

   # the previous script used element instead of ensemble_number.  make them
   # the same for now.
   set element = $ensemble_number
   echo "starting ${myname} for ens member $element at "`date` >> cam_out_temp

   # get model state initial conditions for this ensemble member
   ${LINK} ${CENTRALDIR}/$input_file dart_restart

   # get filter namelists for use by cam
   ${COPY} ${CENTRALDIR}/input.nml input.nml

   # this just creates a file that helps you figure out which member is
   # being advanced in this directory. FYI only, you don't need it.
   echo $element >! element
   cp element element$element

   echo "ls $temp_dir for element $element" >> cam_out_temp
   ls -lRt                                  >> cam_out_temp
    
   # Need a base CAM initial file into which to copy state vector from filter.
   # c[al]minput_$element also carry along CAM/CLM fields which are not updated
   #  by the filter (not part of the filter model state).
   # Look for c[al]minput.nc resulting from the previous advance of this ensemble
   #  member from within the same day/obs_seq.out time span (in CENTRALDIR)
   # When starting an experiment which has no spun up set of members, you
   #  should have already copied the single CAM initial file (e.g. campinput_0.nc)
   #  from the CAM build directory into the CENTRALDIR and made N copies of it.
   
   if (-e     ${CENTRALDIR}/caminput_${element}.nc) then
      ${COPY} ${CENTRALDIR}/caminput_${element}.nc caminput.nc
      echo "caminput comes from ${CENTRALDIR}/caminput_${element}.nc" >> cam_out_temp
   else
      echo ERROR - need $CENTRALDIR/caminput_${element}.nc to exist
      echo ERROR - need $CENTRALDIR/caminput_${element}.nc to exist >> cam_out_temp
      exit -${element}
   endif
   
   if ( -e     ${CENTRALDIR}/clminput_${element}.nc) then
       ${COPY} ${CENTRALDIR}/clminput_${element}.nc clminput.nc
   else
      echo ERROR - need $CENTRALDIR/clminput_${element}.nc to exist
      echo ERROR - need $CENTRALDIR/clminput_${element}.nc to exist >> cam_out_temp
      exit -${element}
   endif
   
   if ( -e     ${CENTRALDIR}/iceinput_${element}.nc) then
       ${COPY} ${CENTRALDIR}/iceinput_${element}.nc iceinput.nc
   else
      echo ERROR - need $CENTRALDIR/iceinput_${element}.nc to exist
      echo ERROR - need $CENTRALDIR/iceinput_${element}.nc to exist >> cam_out_temp
      exit -${element}
      # or if no ice restart file available; start it with ice_in = 'default' 
      # via existence of iceinput in run-cam.csh
   endif
   
   # topography infomation
   ${LINK} ${CENTRALDIR}/cam_phis.nc .

   # translate DART state vector into a CAM caminput.nc file, and create an
   # ascii 'times' file, which will be used to set the namelist for cam to tell
   # it how far to advance the model.
   if (-e dart_restart && -e ${CENTRALDIR}/dart_to_cam) then
      echo ' '                                           >> cam_out_temp
      echo 'advance_model: executing dart_to_cam '`date` >> cam_out_temp
      ${CENTRALDIR}/dart_to_cam                          >> cam_out_temp
      ls -lt                                             >> cam_out_temp
      ${COPY} times ${CENTRALDIR}
   else
      echo "ERROR: either dart_restart file for $element or dart_to_cam not available" >> cam_out_temp
      exit -${element}
   endif
   
   # advance cam 
   echo executing: ${CENTRALDIR}/run-cam.csh ${case}-$element $model ${CENTRALDIR}  >> cam_out_temp
   ${CENTRALDIR}/run-cam.csh ${case}-$element $model ${CENTRALDIR}  >>& cam_out_temp
   
   grep 'END OF MODEL RUN' cam_out_temp > /dev/null
   if ($status == 0) then
      # Extract the new state vector information from the new caminput.nc and
      # put it in '$output_file' (time followed by state)
      echo ' '                           >> cam_out_temp
      echo 'Executing cam_to_dart'       >> cam_out_temp
      ${CENTRALDIR}/cam_to_dart          >> cam_out_temp
   
      # Move updated state vector and new CAM/CLM initial files back to experiment
      # directory for use by filter and the next advance.

      ${MOVE} dart_ics        ${CENTRALDIR}/$output_file
      ${MOVE} namelist        ${CENTRALDIR}
      ${MOVE} caminput.nc     ${CENTRALDIR}/caminput_${element}.nc
      ${MOVE} clminput.nc     ${CENTRALDIR}/clminput_${element}.nc
      ${MOVE} iceinput.nc     ${CENTRALDIR}/iceinput_${element}.nc
   
      echo "finished ${myname} for ens member $element at "`date` >> cam_out_temp
      ${COPY} cam_out_temp ${CENTRALDIR}/H${hour}/cam_out_temp$element
      ${MOVE} cam_out_temp ${CENTRALDIR}/cam_out_temp$element
   else
      echo "WARNING - CAM $element stopped abnormally"
      echo "WARNING - CAM $element stopped abnormally" >> cam_out_temp
      echo "=========================================" >> cam_out_temp
      exit -${element}
   endif
end

   # if this process needs to advance more than one model, read the next set of
   # filenames and ensemble number at the top of this loop.

   @ state_copy++
   @ ensemble_number_line = $ensemble_number_line + 3
   @ input_file_line      = $input_file_line + 3
   @ output_file_line     = $output_file_line + 3
end

cd ${CENTRALDIR}

${REMOVE} $temp_dir/*

# Remove the filter_control file to signal completion
\rm -rf $control_file

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

