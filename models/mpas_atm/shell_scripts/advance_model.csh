#!/bin/csh
#set echo
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Shell script to run the MPAS-A(tmostphere) model from DART input.
#
# This script is called by 'filter' or 'perfect_model_obs'
# after the analysis step is done at each cycle.
#
# This script performs the following:
# 1.  Creates a temporary directory to run an MPAS-A realization (see options)
# 2.  Copies or links the files necessary for the model run into the temporary directory
# 3.  Converts DART state vectors to the mpas input at the beginning
# 4.  Updates an MPAS namelist from a template with new dates
# 5.  Runs the MPAS-A model in a restart mode until target time is reached
# 7.  Checks for incomplete runs
# 8.  Converts mpas output to a DART binary file for the next analysis cycle
# 9.  Saves the mpas analysis file for each member if save_analysis = true
#10.  Saves the mpas forecast file for each member if save_forecast = true
#11.  Saves u on the edges if save_analysis = false and use_u_wind = false (for extended forecasts later)
#
# Note: 1. MPAS is run in a restart mode during the cycles, which means
#          both input and output will be restart files.
#          This also means that one cannot delete temp_dir since we need a
#          place to keep the restart file for each member during the cycle.
#       2. If save_analysis = false and save_forecast = false,
#          the mpas analysis file will be overwritten by 
#          the mpas forecast file to be used as a background for the next analysis cycle.
#       3. dart_to_model expects to have advance_time_present = .true. in input.nml
#          to generate 'mpas_time' for the current and target time info for the forecast run.
#       4. For the required data to run this script, check the section of 'dependencies'.
#
# Arguments for this script (created by 'filter' or 'perfect_model_obs') are:
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

set      process = $1
set   num_states = $2
set control_file = $3

#----------------------------------------------------------------------
# Block 1: copy necessary input files/executables/files common
#          to all model advances to a clean, temporary directory.
#          These will be used by ALL of the ensemble
#          members being advanced by this script.
#----------------------------------------------------------------------

# The run-time directory for the entire experiment is called CENTRALDIR;
set CENTRALDIR = `pwd`

# Do you want to save the analysis file?
set save_analysis = true

# Do you want to save the forecast at the target time?
set save_forecast = false

# Do you want to save normal speed (u) in the forecast at the target time?
set save_prior_uedge = true 

#
set  REMOVE = 'rm -rf'
set    COPY = 'cp -pf'
set    MOVE = 'mv -f'
set    LINK = 'ln -sf'
unalias cd
unalias ls

# if process 0 go ahead and check for dependencies here
if ( $process == 0 ) then

   foreach fn ( advance_time dart_to_model model_to_dart )
   if ( ! -x ${CENTRALDIR}/$fn ) then
     echo ABORT\: advance_model.csh could not find required executable dependency ${CENTRALDIR}/$fn
     exit 1
   endif
   end

   if ( ! -d ${CENTRALDIR}/MPAS_RUN ) then
      echo ABORT\: advance_model.csh could not find required data directory ${CENTRALDIR}/MPAS_RUN, 
      echo         which contains all the MPAS run-time input files
      exit 1
   endif

   if ( ! -x ${CENTRALDIR}/MPAS_RUN/atmosphere_model ) then
     echo ABORT\: advance_model.csh could not find required executable dependency 
     echo         ${CENTRALDIR}/MPAS_RUN/atmosphere_model    
     exit 1
   endif

   if ( ! -r ${CENTRALDIR}/input.nml ) then
     echo ABORT\: advance_model.csh could not find required readable dependency ${CENTRALDIR}/input.nml
     exit 1
   endif

   if ( ! -r ${CENTRALDIR}/namelist.atmosphere ) then
     echo ABORT\: advance_model.csh could not find required readable dependency ${CENTRALDIR}/namelist.atmosphere
     exit 1
   endif

endif # process 0 dependency checking


# Loop through each state
set state_copy = 1
set ensemble_member_line = 1
set      input_file_line = 2
set     output_file_line = 3

while($state_copy <= $num_states)
   
   set ensemble_member = `head -$ensemble_member_line ${CENTRALDIR}/$control_file | tail -1`
   set input_file      = `head -$input_file_line      ${CENTRALDIR}/$control_file | tail -1`
   set output_file     = `head -$output_file_line     ${CENTRALDIR}/$control_file | tail -1`
   
   # Create a new temp directory for each member unless requested to keep and it exists already.
   set temp_dir = 'advance_temp'${ensemble_member}

   if(! -d $temp_dir) mkdir -p $temp_dir                || exit 1
   cd $temp_dir  

   # Get the program and necessary files for the model
   ${LINK} ${CENTRALDIR}/atmosphere_model     .         || exit 2
   ${LINK} ${CENTRALDIR}/streams.atmosphere   .         || exit 3
   ${LINK} ${CENTRALDIR}/stream_list.atmosphere.* .     || exit 4
   ${LINK} ${CENTRALDIR}/*BL                  .         || exit 5
   ${LINK} ${CENTRALDIR}/*DATA                .         || exit 6
   ${LINK} ${CENTRALDIR}/advance_time         .         || exit 7

   # Get the namelists
   if( -e input.nml) ${REMOVE} input.nml
   ${COPY} ${CENTRALDIR}/input.nml      .               || exit 8
   ${COPY} ${CENTRALDIR}/namelist.atmosphere namelist.atmosphere.template  || exit 9

   # Get the grid info files - now for PIO
   set is_grid_info = `grep config_block_decomp_file_prefix namelist.atmosphere.template | wc -l`
   if( $is_grid_info != 1 ) then
       echo Cannot find grid info. Stop.
       exit
   endif
   set fs_grid = `grep config_block_decomp_file_prefix namelist.atmosphere.template | awk '{print $3}' | sed -e "s/'//g"`
   set  x_grid = `echo $fs_grid | cut -d . -f1`
   set  n_grid = `echo $fs_grid | cut -d . -f2`
   set  f_grid = $x_grid.$n_grid
   ${LINK} ${CENTRALDIR}/${fs_grid}* .

   # Surface update - FIXME: No idea how to grab the surface file name (defined in streams.atmosphere)
   set if_sfc_update = `grep config_sst_update namelist.atmosphere.template | awk '{print $3}'`
   if($if_sfc_update == .true.) then
      set fsfc = `ls -1 ${CENTRALDIR}/*sfc_update*.nc | tail -1`
      if(-e $fsfc) then
         ${LINK} $fsfc .
         ls -lL $fsfc						|| exit
      else
         echo $fsfc does not exist.  || exit
      endif
   endif

   # Get the in/out file names for converters and the model
   set f1 = `grep  dart_to_model_input_file input.nml | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g"`
   set f2 = `grep model_to_dart_output_file input.nml | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g"`
   set f3 = `grep   model_analysis_filename input.nml | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g"`
   set f4 = `grep  grid_definition_filename input.nml | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g"`

   if(! -e $f4) ${LINK} ${CENTRALDIR}/$f4 .

   #----------------------------------------------------------------------
   # Block 2: move/convert the DART state vector to the model netcdf file.
   #----------------------------------------------------------------------
   set ff = `basename $f3 .nc`
   set fn = ${ff}.e${ensemble_member}.nc

   if(! -e ${f3}) then
      if(! -e ${CENTRALDIR}/$fn) then
         echo ABORT\: ${CENTRALDIR}/$fn does not exist.
         exit
      endif
      echo ${COPY} ${CENTRALDIR}/$fn ${f3}      
           ${COPY} ${CENTRALDIR}/$fn ${f3}           || exit 2
   endif
   
   ${MOVE} ${CENTRALDIR}/$input_file $f1 || exit 2
   
   # Overwrite a template file (or prior) with the analysis from filter.
   # That is, f3 is updated by f1 here.
   ${CENTRALDIR}/dart_to_model >&! out.dart_to_model
    
   # The program dart_to_model has created an ascii file named mpas_time.
   # Time information is extracted from the file.
   set curr_utc = `head -1 mpas_time | tail -1`		|| exit
   set targ_utc = `head -2 mpas_time | tail -1`
   set intv_utc = `head -3 mpas_time | tail -1`
   set targ_pio = `echo ${targ_utc} | sed -e 's/:/\./g'`

   # Sanity check - if there is any Not-A-Number of Infinity values.
   set drange = `grep "min/max" out.dart_to_model | awk '{print $4,$5}'`
   set is_nan = `echo $drange | grep '[a-zA-Z]' | wc -l`
   if($is_nan == 1) then
      echo Error in dart_to_model: Something bad happened in the analysis. Stop.
      exit
   endif

   ${MOVE} out.dart_to_model out.dart_to_model.${curr_utc}.txt

   # Rename the restart file for PIO
   set f3new = ${ff}.${curr_utc}.nc
   set f3pio = `echo ${f3new} | sed -e 's/:/\./g'`

   # In a restart mode, the initial file name should be provided in a new form 
   # including the date info in PIO format. On the other hand, DART/input.nml
   # does not use the date info in the model_analysis_filename.
   # So we just make a link between the two different filenames here.
   ${LINK} $f3 ${f3pio}

   #----------------------------------------------------------------------
   # Block 3: advance the model
   #          Make sure the file name is consistent in the namelist.atmosphere.
   #          Mar-21-2013: To save the variables at pressure levels,
   #          we print out output files as well. 
   #----------------------------------------------------------------------
   cat >! script.sed << EOF
   /config_start_time/c\
   config_start_time = '$curr_utc'
   /config_stop_time/c\
   config_stop_time = '$targ_utc'
   /config_run_duration/c\
   config_run_duration = '$intv_utc'
EOF
   sed -f script.sed namelist.atmosphere.template >! namelist.atmosphere

   # clean out any old rsl files
   if ( -e log.0000.out ) ${REMOVE} log.*

   # mpi run on Yellowstone
   mpirun.lsf ./atmosphere_model        || exit 3
   sleep 60

   # Check the status
   ls -lrt > list.${curr_utc}.txt
  
   # Check if it is in a restart mode or not
   set if_restart = `grep config_do_restart   namelist.atmosphere | awk '{print $3}'`
   set if_cycling = `grep config_do_DAcycling namelist.atmosphere | awk '{print $3}'`
   if($if_restart == .false. && $if_cycling == .false.) then
      set fh = `basename $ff .init`
      set ff = $fh.restart
   endif

   # Model output at the target time
   set fout = ${ff}.`echo ${targ_utc} | sed -e 's/:/\./g'`.nc
   set date_utc = `ncdump -v xtime $fout | tail -2 | head -1 | cut -d";" -f1 | sed -e 's/"//g'`
   set targ_grg = `echo $date_utc 0 -g | advance_time`
   set targ_day = $targ_grg[1]
   set targ_sec = $targ_grg[2]

   # Check if the model was succefully completed.
   set fatal_err = `grep -i abort log.0000.err | cat | wc -l`
   if(($fatal_err != 0) || ($date_utc != $targ_utc)) then
      mv log.0000.err log.failed.${targ_utc}.txt
      echo $ensemble_member >>! ${CENTRALDIR}/blown_${targ_utc}.out
      echo "Model failure! Check file " ${CENTRALDIR}/blown_${targ_utc}.out
      exit 1
   endif

   #-------------------------------------------------------------------
   # Back up some fields before cleaning up.
   #-------------------------------------------------------------------
   set f3utc = `echo ${curr_utc} 0 | advance_time`
   set if_u_used = `grep use_u_for_wind input.nml | awk '{print $3}' | cut -d ',' -f1`
   if ( $if_u_used == .false. ) then
        set f_wind = analysis.uedge.${f3utc}.nc
        echo ncks -O -v xtime,u $f3 ${f_wind}
             ncks -O -v xtime,u $f3 ${f_wind}
   else
        set f_wind = analysis.uv.${f3utc}.nc
        ncks -O -v xtime,uReconstructZonal,uReconstructMeridional $f3 ${f_wind}
   endif
   ls -l ${f_wind}

   set f3out = mpas_anal.${f3utc}.nc
   ${MOVE} $f3 ${f3out}

   if( ! -e ${f_wind} ) then
       if(! -e ${f3out} ) then
	  echo Neither ${f_wind} or ${f3out} exists. We stop here.
          echo You cannot run extended forecasts from the mean EnKF analysis at this cycle later. 
          exit
       else 
          set date_anl = `ncdump -v xtime ${f3out} | tail -2 | head -1 | cut -d";" -f1 | sed -e 's/"//g'`
          set    tanal = `echo $date_anl 0 | advance_time`
          if($tanal != ${f3utc}) then
             echo ${f3out} has wrong time. $tanal should be ${f3utc}. Stop.
             exit
          endif  
       endif
   endif
   if( -e ${f3pio}) ${REMOVE} ${f3pio}

   #-------------------------------------------------------------------
   # Block 4: Convert your model output to a DART format ics file,
   #          then move it back to CENTRALDIR
   #          We also want to keep $f3 for the next cycle under this 
   #          temp directory
   #-------------------------------------------------------------------
   # Overwrite the analysis file with the forecast at target time (for the next cycle).
   set futc = `echo ${targ_utc} 0 | advance_time`
   set fsav = mpas_fcst.${futc}.nc
   if ( $save_prior_uedge == true ) then
        if ( $if_u_used == .false. ) then
        set f_wind = forecast.uedge.${futc}.nc
        echo ncks -O -v xtime,u $fout ${f_wind}
             ncks -O -v xtime,u $fout ${f_wind}
        else
        set f_wind = forecast.uv.${futc}.nc
        ncks -O -v xtime,uReconstructZonal,uReconstructMeridional $fout ${f_wind}
        endif
   endif
   if ( $save_forecast == true ) then
        echo ${COPY} $fout $fsav
             ${COPY} $fout $fsav	|| exit 5
   endif

   # Get ready to run with full state vectors in a restart mode from the second cycle.
   if($if_restart == .false. && $if_cycling == .false.) then
      ${COPY} ${CENTRALDIR}/input.nml.filter input.nml || exit
      set f3 = `grep model_analysis_filename input.nml | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g"`
   endif

   echo ${MOVE} $fout $f3
        ${MOVE} $fout $f3   	 		|| exit 6

   ${CENTRALDIR}/model_to_dart >&! out.model_to_dart.${date_utc}

   # Sanity check - if there is any Not-A-Number or Infinity values.
   set vrange = `grep "min/max" out.model_to_dart.${date_utc} | awk '{print $4,$5}'`
   set is_nan = `echo $vrange | grep '[a-zA-Z]' | wc -l`
   if($is_nan == 1) then
      echo Error in model_to_dart: Something bad happened in mpas forecast. Stop.
      exit
   endif

   echo ${MOVE} $f2 ${CENTRALDIR}/$output_file 
        ${MOVE} $f2 ${CENTRALDIR}/$output_file || exit 7

   if ( $save_analysis == false ) ${REMOVE} ${f3out}

   # Change back to original directory.
   cd $CENTRALDIR

   echo "Ensemble Member $ensemble_member completed"

   # and now repeat the entire process for any other ensemble member that
   # needs to be advanced by this task.
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

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

