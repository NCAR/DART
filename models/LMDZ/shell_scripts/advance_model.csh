#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# Script to advance one ensemble member one filter "time step"
# when the model advance is executed as a separate process.
# Calls run_lmdz.csh, the LMDZ execution script.
#
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


set      process = $1    # the process number of caller
set   num_states = $2    # the number of state copies belonging to that process
set control_file = $3    # the name of the filter_control_file for that process


# Create a unique temporary working directory for this process's stuff
set temp_dir = 'advance_temp'${process}

# The run-time directory for the entire experiment is called CENTRALDIR;
# we need to provide a safe haven for each TASK ... in 'temp_dir'.
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

 echo "advance_model.csh args = $1 $2 $3"                    >  lmdz_out_temp
 echo "CENTRALDIR is ${CENTRALDIR}"                          >> lmdz_out_temp
 echo "temp_dir is $temp_dir"                                >> lmdz_out_temp

# Loop through each state
set ensemble_member_line = 1
set      input_file_line = 2
set     output_file_line = 3


set  state_copy = 1

while($state_copy <= $num_states)
   
   set ensemble_member = `head -$ensemble_member_line ../$control_file | tail -1`
   set input_file      = `head -$input_file_line      ../$control_file | tail -1`
   set output_file     = `head -$output_file_line     ../$control_file | tail -1`

 set element = $ensemble_member
 echo "starting ${myname} for ens member $element at "`date` >> lmdz_out_temp

   # get model state initial conditions for this ensemble member
    ${LINK} ${CENTRALDIR}/$input_file temp_ic
   # get filter namelists for use by lmdz
    ${COPY} ${CENTRALDIR}/input.nml input.nml

   # this just creates a file that helps to figure out which member is
   # being advanced in this directory. FYI only, you don't need it.
       echo $element >! element
       cp element element$element
   
       echo "ls $temp_dir for element $element" >> lmdz_out_temp
       ls -lRt                                  >> lmdz_out_temp

   # Need a base LMDZ initial file into which to copy state vector from filter.
   # Look for start.nc resulting from the previous advance of this ensemble
   # member from within the same day/obs_seq.out time span (in CENTRALDIR)
   # When starting an experiment which has no spun up set of members, you
   # should have already copied the single LMDZ initial file (e.g. start_0.nc)
   # from the LMDZ build directory into the CENTRALDIR and made N copies of it.


   if (-e  ${CENTRALDIR}/start_${element}.nc) then
      ${COPY} ${CENTRALDIR}/start_${element}.nc start.nc
      ${COPY} start.nc start.nc_prior 
      echo "start comes from ${CENTRALDIR}/start_${element}.nc" 
      echo "start comes from ${CENTRALDIR}/start_${element}.nc" >> lmdz_out_temp
   else
      echo ERROR - need $CENTRALDIR/start_${element}.nc to exist
      echo ERROR - need $CENTRALDIR/start_${element}.nc to exist >> lmdz_out_temp
      exit -${element}
   endif

   if (-e  ${CENTRALDIR}/startphy_${element}.nc) then
      ${COPY} ${CENTRALDIR}/startphy_${element}.nc startphy.nc

      echo "startphy comes from ${CENTRALDIR}/startphy_${element}.nc"
      echo "startphy comes from ${CENTRALDIR}/startphy_${element}.nc" >> lmdz_out_temp
   else
      echo ERROR - need $CENTRALDIR/startphy_${element}.nc to exist
      echo ERROR - need $CENTRALDIR/startphy_${element}.nc to exist >> lmdz_out_temp
      exit -${element}
   endif



   # translate DART state vector into a LMDZ start.nc file
    if (-e temp_ic && -e ${CENTRALDIR}/dart_to_lmdz) then
      echo ' '                                            >> lmdz_out_temp
      echo 'advance_model: executing dart_to_lmdz '`date` >> lmdz_out_temp
      echo  ${CENTRALDIR}
     ${CENTRALDIR}/dart_to_lmdz                 
     ${COPY} start.nc start.nc_posterior
      ls -lt                                              >> lmdz_out_temp
    else
      echo "ERROR: either temp_ic file for $element or dart_to_lmdz not available" >> lmdz_out_temp
      exit -${element}
    endif

   # advance LMDZ
   echo executing: ${CENTRALDIR}/run_lmdz.csh $element 
   echo executing: ${CENTRALDIR}/run_lmdz.csh $element  >> lmdz_out_temp
   ${CENTRALDIR}/run_lmdz.csh >& gcm.log
 
   
   grep 'GLOB '               gcm.log  
   grep 'Simulation finished' gcm.log  
   grep 'Everything is cool'  gcm.log 
   grep 'Everything is cool'  gcm.log > /dev/null 

  if ($status == 0) then
   # Extract the new state vector information from the new start.nc and
   # put it in '$output_file' (time followed by state)
    echo ' '                           >> lmdz_out_temp
    echo 'Executing lmdz_to_dart'      >> lmdz_out_temp
    ${CENTRALDIR}/lmdz_to_dart         >> lmdz_out_temp


   # Move updated state vector and new LMDZ initial files back to experiment
   # directory for use by filter and the next advance.
    
     ${MOVE} dart_ics       ${CENTRALDIR}/$output_file
     ${COPY} start.nc     ${CENTRALDIR}/start_${element}.nc
     ${COPY} startphy.nc  ${CENTRALDIR}/startphy_${element}.nc
     echo "finished ${myname} for ens member $element at "`date` >> lmdz_out_temp
     #${COPY} lmdz_out_temp ${CENTRALDIR}/H${hour}/lmdz_out_temp$element
     ${MOVE} lmdz_out_temp ${CENTRALDIR}/lmdz_out_temp$element

   else
     echo "WARNING - LMDZ $element stopped abnormally"
     echo "WARNING - LMDZ $element stopped abnormally" >> lmdz_out_temp
     echo "========================================="  >> lmdz_out_temp
     exit -${element}
   endif


   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3
end

cd ${CENTRALDIR}

#${REMOVE} $temp_dir/*


\rm -rf $control_file

exit 0


