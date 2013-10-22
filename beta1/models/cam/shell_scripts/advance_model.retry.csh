#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
#----------------------------------------------------------------------
# advance_model.csh
#
# Script to advance one ensemble member one filter "time step"
# when the model advance is executed as a separate process.
# Called by the filter executable (for async=2 or 4)
# Calls run-cam.csh, the CAM execution script.
# Calls 3 translation routines to translate time and model state.
# Runs on one of the compute nodes allotted to the filter executable.
#
# Arguments are 
# arg#1  the process number of caller
# arg#2  the number of state copies belonging to that process
# arg#3  the name of the filter_control_file for that process

#----------------------------------------------------------------------

set process = $1
set num_states = $2
set control_file = $3

set retry_max = 2

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
   touch cam_out_temp
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
   #      by the filter (not part of the filter model state).
   # First look for c[al]minput.nc resulting from the previous advance of this ensemble
   #      member from within the same day/obs_seq.out time span (in CENTRALDIR)
   # Failing that, look for the results of the last advance of this ensemble member
   #      of the previous obs_seq.out (i.e. in CENTRALDIR/exp_name/day/CAM)
   # Failing that (when starting an experiment which has no spun up set of members)
   #      get a copy of a single CAM initial file (usually from somewhere independent
   #      of this experiment, i.e. /scratch/.../New_state/T42_GWD/CAM/caminput_0.nc)
   
   if (-e     ${CENTRALDIR}/caminput_${element}.nc) then
      ${COPY} ${CENTRALDIR}/caminput_${element}.nc caminput.nc
      echo "CENTRALDIR caminput comes from ${CENTRALDIR}/caminput_${element}.nc" >> cam_out_temp
   else if (-e ${cam_init}${element}.nc) then
      ${COPY}  ${cam_init}${element}.nc caminput.nc
      echo "cam_init caminput comes from ${cam_init}${element}.nc" >> cam_out_temp
   else
      ${COPY}  ${cam_init}0.nc caminput.nc
      echo "DEFAULT caminput comes from ${cam_init}0.nc" >> cam_out_temp
   endif
   
   if ( -e     ${CENTRALDIR}/clminput_${element}.nc) then
       ${COPY} ${CENTRALDIR}/clminput_${element}.nc clminput.nc
   else if (-e ${clm_init}${element}.nc) then
       ${COPY} ${clm_init}${element}.nc clminput.nc
   else
       ${COPY} ${clm_init}0.nc clminput.nc
   endif
   
# Pond restart files have unchangable names like
# FV_2deg-noleap-O2-Dev20-5-10.cice.r.volpn.2003-06-01-21600
# Get iceinput and the related restart files for meltpond and aero
# Early versions of CAM3.6 (<3.56.59)
#   if ( -e     ${CENTRALDIR}/iceinput_${element}.tar) then
#       tar -x -f ${CENTRALDIR}/iceinput_${element}.tar
#   else if (-e ${ice_init}${element}.tar) then
#       tar -x -f ${ice_init}${element}.tar 
   if ( -e     ${CENTRALDIR}/iceinput_${element}.nc) then
       ${COPY} ${CENTRALDIR}/iceinput_${element}.nc iceinput.nc
   else if (-e ${ice_init}${element}.nc) then
       ${COPY} ${ice_init}${element}.nc iceinput.nc
   else
       # no ice restart file available; start it with ice_in = 'default' via existence of iceinput
       # in run-cam.csh
   endif
   
   ${LINK} ${CENTRALDIR}/cam_phis.nc .

   # Create 'times' file for CAM from DART times in assim_model_state_ic#
   # This info is passed to CAM through the creation of its namelist(s)
   # Also extract state variables from assim_model_state_ic# into the caminput.nc file
   if (-e dart_restart && -e ${CENTRALDIR}/dart_to_cam) then
      echo ' '                           >> cam_out_temp
      echo 'Executing dart_to_cam'       >> cam_out_temp
      ${CENTRALDIR}/dart_to_cam          >> cam_out_temp
      ls -ltR                            >> cam_out_temp
      ${COPY} times ${CENTRALDIR}
   else
      echo "ERROR: either ic file $element or dart_to_cam not available for dart_to_cam" >> cam_out_temp
      exit 1
   endif
   
   # advance cam 
   #   echo executing: ${model:h}/run-cam.csh ${case}-$element $model ${CENTRALDIR} >> cam_out_temp
      set retry = 0
      while ($retry < $retry_max)
         echo executing: ${CENTRALDIR}/run-cam.csh ${case}-$element $model ${CENTRALDIR} \
              >> cam_out_temp
         ${CENTRALDIR}/run-cam.csh ${case}-$element $model ${CENTRALDIR}  >>& cam_out_temp
# kdr added 12/28/2010
         set retry_status = $status
         if ($retry_status == 0) then
            grep 'END OF MODEL RUN' cam_out_temp > /dev/null
            set retry_status = $status
         else 
            echo 'run-cam.csh DIED with status'  >> cam_out_temp
            exit $retry_status
         endif
#         grep 'END OF MODEL RUN' cam_out_temp > /dev/null
#         if ($status == 0) then
         if ($retry_status == 0) then
# kdr end
            set retry = $retry_max
      # Extract the new state vector information from the new caminput.nc and
      # put it in temp_ud (time followed by state)
            echo ' '                           >> cam_out_temp
            echo 'Executing cam_to_dart'       >> cam_out_temp
            ${CENTRALDIR}/cam_to_dart          >> cam_out_temp
   
      # Save CLM and CAM files for storage of analyses in CAM initial file format (analyses2initial)
            # get the forecast time, which is the time of this CLM initial file
            set seconds = (`head -1 times`)
            if ($seconds[2] == 0) then
               set hour = 24
            else
               @ hour = $seconds[2] / 3600
            endif
            if ($hour < 10) set hour = 0$hour
            # Directory for storing CAM initial files until they can be averaged during
            # archiving for analyses
            if (! -d ${CENTRALDIR}/H${hour}) mkdir ${CENTRALDIR}/H${hour}
      # Move updated state vector and new CAM/CLM initial files back to experiment
      # directory for use by filter and the next advance.
      # Store them in H## directories for later generation of  ensemble average CAM and CLM 
      # initial file after all members are done.
            ${MOVE} temp_ud         ${CENTRALDIR}/$output_file
            ${MOVE} namelist        ${CENTRALDIR}
            ${MOVE} caminput.nc     ${CENTRALDIR}/H${hour}/caminput_${element}.nc
            ${MOVE} clminput.nc     ${CENTRALDIR}/H${hour}/clminput_${element}.nc
            ${MOVE} iceinput.nc     ${CENTRALDIR}/H${hour}/iceinput_${element}.nc 
            # Earlier versions      ${MOVE} iceinput.tar ${CENTRALDIR}/H${hour}/iceinput_${element}.tar 

            # link the new initial files into the CENTRAL directory where filter will find them.
            ${LINK} ${CENTRALDIR}/H${hour}/caminput_${element}.nc \
                    ${CENTRALDIR}/caminput_${element}.nc
            ${LINK} ${CENTRALDIR}/H${hour}/clminput_${element}.nc \
                    ${CENTRALDIR}/clminput_${element}.nc
            ${LINK} ${CENTRALDIR}/H${hour}/iceinput_${element}.nc  \
                    ${CENTRALDIR}/iceinput_${element}.nc 
            # Earlier versions
            # ${LINK} ${CENTRALDIR}/H${hour}/iceinput_${element}.tar  \
            #         ${CENTRALDIR}/iceinput_${element}.tar 
   
            echo "finished ${myname} for ens member $element at "`date` >> cam_out_temp
               ${COPY} cam_out_temp ${CENTRALDIR}/H${hour}/cam_out_temp$element
            ${MOVE} cam_out_temp ${CENTRALDIR}/cam_out_temp$element
         else
            @ retry++
            if ($retry < $retry_max) then
# Add section to make CAM write out something every time step during this retry.
# Could be added to casemodel, but be careful of how run-cam.csh uses the number of
# lines in casemodel.
               echo "WARNING - CAM $element stopped abnormally; will be retried"
               echo "WARNING - CAM $element stopped abnormally; will be retried" >> cam_out_temp
               echo "===========================================================" >> cam_out_temp
            else
               set DEADDIR = ${temp_dir}_dead
               echo "WARNING - CAM $element stopped abnormally; see $DEADDIR"
               echo "WARNING - CAM $element stopped abnormally; see $DEADDIR" >> cam_out_temp
               ${COPY} cam_out_temp ${CENTRALDIR}/cam_out_temp${element}_died
               mkdir $DEADDIR
               ${MOVE} * $DEADDIR
               exit -${element}
            endif
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

