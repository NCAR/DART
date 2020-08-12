#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#***********************************************************************************
# Tarkeshwar Singh
# July 2014
# tarkphysics87@gmail.com

# purpose :  prepare all initial conditions and submit filter run 
#
#***********************************************************************************
# Set alias
# #
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
#

# remove any pre existing temporary filles and directories
$REMOVE start_* startphy_*  start.nc assim_model_state_* lmdz_out_temp*
$REMOVE -rf advance_temp*
#--------------------------------------------------------------------------------------
source Control_File.csh
#---------------------------------Copy--------------------------------------------------
#Copy dart exe
$COPY $DART_LMDZ5/work/dart_to_lmdz .
$COPY $DART_LMDZ5/work/lmdz_to_dart .
$COPY $DART_LMDZ5/work/filter .
$COPY $DART_LMDZ5/work/fill_inflation_restart .
$COPY $DART_LMDZ5/work/trans_time .

#Copy shell scripts 
$COPY $DART_LMDZ5/shell_scripts/advance_model.csh .
$COPY $DART_LMDZ5/shell_scripts/run_lmdz.csh .

# Determine the number of ensemble members from input.nml,
 set ENSEMBLESTRING = `grep -A 42 filter_nml input.nml | grep ens_size`
 set ensemble_size  = `echo $ENSEMBLESTRING[3] | sed -e "s#,##"`
 set num_ens        = ${ensemble_size}
 echo "There are ${num_ens} ensemble members."
# Copy final_full.$num_ens for sampling error corrections
 
if ( $sampling_error_correction == 1) then
   $COPY $DART_LMDZ5/../../system_simulation/final_full_precomputed_tables/final_full.$num_ens .
endif
 

#--------copy start_#.nc,startphy_#.nc  & filter_ic_old.### from storage and rename -
set n = 1
while($n <= ${num_ens})

   $COPY $DART_ics/start_$n.nc .
   $COPY $DART_ics/startphy_$n.nc .

   set from = `printf "%s.%04d" $DART_ics/filter_ic_new $n`
   set to   = `printf "%s.%04d" filter_ic_old        $n`
   echo copying $from to $to
   $COPY $from $to
   @ n++
end

$COPY start_1.nc start.nc
#----------------- alter the data timestamp in a DART restart file filter_ic_old --
# fill required time in restart_file_tool_nml inside input.nml

#if ($change_timestamp == 1) then
#  ./restart_file_tool
#endif


#---------------------------------------------------------------------------------
#create inflation restart files so the values inf_initial_from_restart and 
#inf_sd_initial_from_restart items in the &filter_nml namelist can be .TRUE. 
#from the beginning. 

 echo $inf_initial $inf_sd_initial | ./fill_inflation_restart
 $COPY inflate_ics prior_inf_ic_old 


#                      ------  MAIN JOB -----------
#
#************************ LOOP OVER OBSERVATIONS ***********************************************
set day_count = 1

foreach file (`cat $obs_seq_list`)

   echo "=================================================================="
   echo "Assimilating Observation file $file"
   echo "=================================================================="
   set date = `echo $file | cut -c8-16`
   mkdir OUTPUT_$date

   #************************************************************************************************
   # Store all restart files at give day frequency : $restart_store_freq
   if ( $day_count % $restart_store_freq == 0 ) then
     echo 'Storing Filter restart files in OUTPUT_'$date' directory' 
     $COPY start_*.nc        OUTPUT_$date
     $COPY startphy_*.nc     OUTPUT_$date
     $COPY filter_ic_new.*   OUTPUT_$date
     $COPY prior_inf_ic_new  OUTPUT_$date 
   endif
   @ day_count++
 
  #-----------
  # Link obs_seq.out file  to be assimilated 
  $REMOVE obs_seq.out
  $LINK $OBS_PATH/$file obs_seq.out

  #-------run filter--------
  ulimit -s unlimited
  date > time_filter
   /data/opt/mpi/openmpi-1.6.3/bin/mpirun --hostfile $Host_File_Path/$host_file -np $num_proc filter >filter.log
  date >> time_filter

  grep 'filter: End of main filter assimilation loop'  filter.log
  grep 'filter: End of main filter assimilation loop'  filter.log > /dev/null

  #-------copy outputs file in OUTPUT_#### directory ----------------------------------------------
  if ($status == 0) then
    echo "Filter Finished for observation $file"
    $MOVE analysis.nc       OUTPUT_$date
    $MOVE preassim.nc       OUTPUT_$date
    $MOVE obs_seq.final     OUTPUT_$date
    $MOVE histhf*.nc*       OUTPUT_$date
    $MOVE histins*.nc*      OUTPUT_$date
    $MOVE filter.log        OUTPUT_$date
    $MOVE dart_log.out      OUTPUT_$date
  else
    echo "WARNING - FILTER  stopped abnormally $file"
    echo "=========================================" 
    exit -${file}
  endif

  #-------move prior_inf_ic_new to prior_inf_ic_old for next restart run --------------------------
  if (-e  prior_inf_ic_new ) then
    $COPY prior_inf_ic_new prior_inf_ic_old
   else
    echo ERROR - need prior_inf_ic_new to exist 
    exit -$file
  endif

  #---- move filter_ic_new.## to filter_ic_old.## for next restart run -----------------------------
  set n = 1
  while($n <= ${num_ens})
      set from = `printf "%s.%04d" filter_ic_new $n`
      set to   = `printf "%s.%04d" filter_ic_old $n`
      echo copying $from to $to
      $COPY $from $to
      @ n++
  end

  #************************************************************************************************
end  # loop $file

exit 0


