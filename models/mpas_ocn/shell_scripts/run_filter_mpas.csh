#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
####################################################################################
#
#  run_filter_mpas.csh
#
#  THIS IS A TOP-LEVEL DRIVER SCRIPT TO RUN FILTER AND MPAS MODEL IN AN MPI VERSION.
#
#  This is a sample script for a retrospective case study, and was tested on
#  NCAR IBM Supercomputer (bluefire) using a "bsub" command.
#  Any back-up option to a mass storage is not supported here. 
#  Instead, all the outputs will be locally stored. 
#  Check your disk space before running this script.
#
#  Required file to run this script:
#  1. namelist.input          (for mpas run)
#  2. advance_model.csh       (for mpas run)
#  3. mpas_init.nc            (for mpas run) - you can change the filename below (mpas_fname).
#  4. input.nml.template      (for filter run)  
#  5. filter.template.lsf     (for mpi filter run)
#  6. advance_model.template.lsf (for mpi mpas run)
#  You may need to check other namelist configuration in your input.nml.template
#  before running this script. See 'Note' below.
#
#  Input files to run filter:
#  1. FG_DIR/filter_ics - first guess in dart format for the initial cycle
#  2. OBS_DIR/YYYYMMDDHH/obs_seq.out - obs sequence file for each analysis cycle
#
#  Output from this script will be all stored in $RUN_DIR/${expname}/YYYYMMDDHH/.
#
#  Written by So-Young Ha (MMM/NCAR) Sep-30-2011
#
####################################################################################
# Experiment name and the cycle period
#--------------------------------------------------------------------------
set  expname = Test2			# experiment name
set date_ini = 2010-10-23_06:00:00	# initial cycle for this experiment
set date_beg = 2010-10-23_06:00:00	# start date to run this script
set date_end = 2010-10-23_07:00:00	# end date to run this script
set intv_day = 0			# cycling frequency - assimilation_period_days    in input.nml
set intv_sec = 3600			# cycling frequency - assimilation_period_seconds in input.nml
#--------------------------------------------------------------------------
# Directories
#--------------------------------------------------------------------------
set DART_DIR = /mmm/mmmtmp/syha/MPAS/DART/models/mpas_atm/work	# where dart executables exist
set MPAS_DIR = /mmm/mmmtmp/syha/MPAS             		# where all the aux files exist to run mpas
set  OBS_DIR = OBS_SEQ						# where all obs sequence files exist
set   FG_DIR = MPAS_FG 						# where filter_ics exists for initial cycle
set  RUN_DIR = `pwd`						# top-level working directory
set  OUT_DIR = $RUN_DIR/$expname				# storage for output
#--------------------------------------------------------------------------
# Assimilation parameters
# -----------------------
# Note: For your own complete filter design, you may want to check your input.nml
# for the parameters which are not set up in this section. At least it would be nice to
# double-check &filter_nml, &obs_kind_nml, &model_nml, &location_nml and &mpas_vars_nml,
# then rename it as input.nml.template before running this script.
#--------------------------------------------------------------------------
set nens            = 2		# ensemble size for ens_size in input.nml
set cutoff          = 0.05	# horizontal location - cutoff in input.nml
set horiz_dist_only = true	# horizontal localization only (true or false)
set adaptive_inf    = true	# adaptive_inflation - If true, this script only supports 
 				# spatially-varying state space prior inflation.
				# And you also need to edit inf_sd_initial, inf_damping,
				# inf_lower_bound, and inf_sd_lower_bound in &filter_nml.
set single_restart = false	# true if all copies read from and written to a single file in filter.
#--------------------------------------------------------------------------
# LSF setup on NCAR bluefire (only if $bsub_in_ibm = yes)
#--------------------------------------------------------------------------
set  bsub_in_ibm = yes  	# run on bluefire? yes or no.
set  PROJ_NUMBER = xxxxxxxx	# Account key to submit filter and mpas jobs in bluefire
set  time_filter = 00:30	# wall clock time for mpi filter runs
#--------------------------------------------------------------------------
# File naming convention
#--------------------------------------------------------------------------
set  mpas_fname = mpas_init.nc
set restart_in  = filter_ics
set restart_out = filter_restart
set    infl_in  = prior_inflate_ics
set    infl_out = prior_inflate_restart
set obs_seq_in  = obs_seq.out
set obs_seq_out = obs_seq.final

####################################################################################
# End of User Specificaton
####################################################################################

set REMOVE = '/bin/rm -rf'
set   COPY = 'cp -pf'
set   MOVE = 'mv -f'
set   LINK = 'ln -sf'
unalias cd
unalias ls

#------------------------------------------
# Check if we have all the necessary files.
#------------------------------------------
foreach fn ( input.nml.template namelist.input )
  if ( ! -r $fn ) then
      echo ABORT\: We cannot find required readable dependency $fn.
      exit -1
  endif 
end
foreach fn ( filter advance_time convertdate )
  if ( ! -x $fn ) then
       echo ABORT\: We cannot find required executable dependency $fn.
       exit -1
  endif
end 
if ( ! -e advance_model.csh) then
    echo ${COPY} $DART_DIR/../shell_scripts/advance_model.csh .
         ${COPY} $DART_DIR/../shell_scripts/advance_model.csh .
endif
if($bsub_in_ibm == yes) then
 foreach fn ( filter.template.lsf advance_model.template.lsf )
  if ( ! -r $fn ) then
      echo ABORT\: We cannot find required readable dependency $fn.
      exit -1
  endif 
 end
endif

if ( ! -d MPAS_RUN ) ${LINK} $MPAS_DIR MPAS_RUN		# for advance_model.csh
if ( ! -d logs ) mkdir logs				# to print out log files
${COPY} input.nml.template input.nml			# for advance_time

#------------------------------------------
# Time info
#------------------------------------------
set greg_ini = `echo $date_ini 0 -g | advance_time`
set greg_beg = `echo $date_beg 0 -g | advance_time`
set greg_end = `echo $date_end 0 -g | advance_time`
set  intv_hr = `expr $intv_sec \/ 3600`
set  intv_dh = `expr $intv_day \* 24`
   @ intv_hr += $intv_dh
set diff_day = `expr $greg_end[1] \- $greg_beg[1]`
set diff_sec = `expr $greg_end[2] \- $greg_beg[2]`
set diff_tot = `expr $diff_day \* 86400 \+ $diff_sec`
set n_cycles = `expr $diff_tot \/ $intv_sec \+ 1`
echo Total of $n_cycles cycles from $date_beg to $date_end will be run every $intv_hr hr.

set icyc = 1
if($date_beg != $date_ini) then
   set init_day = `expr $greg_beg[1] \- $greg_ini[1]`
   set init_sec = `expr $greg_beg[2] \- $greg_ini[2]`
   set init_dif = `expr $init_day \* 86400 \+ $init_sec`
   set icyc = `expr $init_dif \/ $intv_sec \+ 1`
endif
set ncyc = `expr $icyc \+ $n_cycles \- 1`
   
if($icyc == 1) ${LINK} $FG_DIR/$mpas_fname .

#--------------------------------------------------------
# Cycling gets started
#--------------------------------------------------------
set time_anl = `echo $date_beg 0 | advance_time`		#YYYYMMDDHH

while ( $icyc <= $ncyc )

  set time_pre = `echo $time_anl -$intv_hr | advance_time`	#YYYYMMDDHH
  set time_nxt = `echo $time_anl +$intv_hr | advance_time`	#YYYYMMDDHH
  set greg_obs = `echo $time_anl 0 -g | advance_time`
  set greg_obs_days = $greg_obs[1]
  set greg_obs_secs = $greg_obs[2]
  echo Cycle $icyc at $time_anl\: ${greg_obs_days}_${greg_obs_secs}

  #------------------------------------------------------
  # 1. Namelist setup
  #------------------------------------------------------
   cat >! script.sed << EOF
  /ens_size /c\
   ens_size                    = $nens,
  /cutoff  /c\
   cutoff                      = $cutoff,
  /horiz_dist_only/c\
   horiz_dist_only             = .${horiz_dist_only}.,
  /init_time_days/c\
   init_time_days           = -1,
  /init_time_seconds/c\
   init_time_seconds        = -1,
  /first_obs_days/c\
   first_obs_days           = -1,
  /first_obs_seconds/c\
   first_obs_seconds        = -1,
  /last_obs_days/c\
   last_obs_days            = -1,
  /last_obs_seconds/c\
   last_obs_seconds         = -1,
   /model_analysis_filename/c\
   model_analysis_filename      = '$mpas_fname',
   /grid_definition_filename/c\
   grid_definition_filename     = '$mpas_fname',
   /assimilation_period_days/c\
   assimilation_period_days     = $intv_day,
   /assimilation_period_seconds/c\
   assimilation_period_seconds  = $intv_sec,
   /single_restart_file_in/c\
   single_restart_file_in  = .$single_restart.,
   /single_restart_file_out/c\
   single_restart_file_out = .$single_restart.,
EOF

  if( $adaptive_inf == true ) then       # We are going to have a prior inflation only.

   cat >! script1.sed << EOF
   /inf_flavor /c\
   inf_flavor                  = 2,                       0,
EOF

   if($icyc == 1) then
   cat >! script2.sed << EOF
   /inf_initial_from_restart/c\
   inf_initial_from_restart    = .false.,                .false.,
   /inf_sd_initial_from_restart/c\
   inf_sd_initial_from_restart = .false.,                .false.,
EOF
   else
   cat >! script2.sed << EOF
   /inf_initial_from_restart/c\
   inf_initial_from_restart    = .true.,                .true.,
   /inf_sd_initial_from_restart/c\
   inf_sd_initial_from_restart = .true.,                .true.,
EOF
   endif

   cat script1.sed >> script.sed
   cat script2.sed >> script.sed

  else   # turn off the adaptive inflation in prior

   cat >! script1.sed << EOF
  /inf_flavor /c\
   inf_flavor                  = 0,                       0,
EOF
   cat script1.sed >> script.sed

  endif

   sed -f script.sed input.nml.template >! input.nml

   #------------------------------------------------------
   # 2. Link to a restart file to get filter started 
   # (assuming start_from_restart = .true. in input.nml)
   #------------------------------------------------------

   if($icyc == 1) then
      set dir_rst = $FG_DIR
   else
      set dir_rst = $RUN_DIR/${expname}/${time_anl}
   endif
   set fn_rst = $dir_rst/${restart_in}

   if($single_restart == true) then
      if(! -e $fn_rst) then
         echo $fn_rst does not exist. Stop.
         exit
      else
         ${LINK} ${fn_rst} ${restart_in}
      endif
   else
      set i = 1
      while ( $i <= $nens )
        set icnum = `echo $i + 10000 | bc | cut -b2-5`
        if(! -e ${fn_rst}.${icnum}) then
           echo ${fn_rst}.${icnum} does not exist. Stop.
           exit
        else
           ${LINK} ${fn_rst}.${icnum} ${restart_in}.${icnum}
        endif
        @ i++
      end 
   endif

   if( $adaptive_inf == true && $icyc > 1 ) then
       if(! -e $RUN_DIR/${expname}/${time_pre}/${infl_out}) then
          echo $RUN_DIR/${expname}/${time_pre}/${infl_out} does not exist. Stop.
          exit
       endif
       ${LINK} $RUN_DIR/${expname}/${time_pre}/${infl_out} ${infl_in}
   endif

   #------------------------------------------------------
   # 3. Obs sequence for this analysis cycle
   #------------------------------------------------------
   set fn_obs = $OBS_DIR/${time_anl}/obs_seq.out
   if(! -e $fn_obs) then
      echo $fn_obs does not exist. Stop.
      exit
   endif
   ${LINK} $fn_obs .

   #------------------------------------------------------
   # 4. Run filter
   #------------------------------------------------------
   set job_name = ${expname}.${icyc}
   echo Running filter: $job_name

  if($bsub_in_ibm == yes) then

   cat >! filter.sed << EOF
   s#JOB_NAME#${job_name}#g
   s#PROJ_NUMBER#${PROJ_NUMBER}#g
   s#NENS#${nens}#g
   s#JOB_TIME#${time_filter}#g
EOF

   sed -f filter.sed filter.template.lsf >! filter.lsf
   bsub < filter.lsf

   # Wait until the job is finished.
   set is_there = `bjobs -w | grep $job_name | wc -l`
   while ( $is_there != 0 )
    sleep 60
    set is_there = `bjobs -w | grep $job_name | wc -l`
   end

  else

   echo `date +%s` >&! filter_started
   ./filter >! filter.log
   if ( -e obs_seq.final )  touch filter_done

  endif

   # Check errors in filter.
   if ( -e filter_started && ! -e filter_done ) then
        echo "Filter was not normally finished. Exiting."
        ${REMOVE} filter_started
        exit
   endif

   ${REMOVE} filter_started filter_done
   echo Filter is done for Cycle ${icyc}\: $time_anl

   #------------------------------------------------------
   # 5. Target time for model advance
   #------------------------------------------------------
   @ greg_obs_secs += $intv_sec
   if ($greg_obs_secs >= 86400) then
       @ greg_obs_secs -= 86400
       @ greg_obs_days += 1
   endif
   set datenew = `convertdate | tail -1 |cut -d: -f2` << EOF
2
$greg_obs_days $greg_obs_secs
EOF
   set date_nxt = `echo ${datenew[1]}-${datenew[2]}-${datenew[3]}_${datenew[4]}:${datenew[5]}:${datenew[6]}`
   echo Target date: $date_nxt $greg_obs_days $greg_obs_secs
   ${MOVE} input.nml input.nml.filter.${icyc}

   #------------------------------------------------------
   # 6. Convert the analysis to the dart format for each member
   #------------------------------------------------------
   ${REMOVE} assim_model_state_ic.*
   cat >! input.nml.ic.restart_file_tool << EOF
   input_file_name              = "$restart_out",
   output_file_name             = "assim_model_state_ic",
   ens_size                     = $nens,
   single_restart_file_in       = .$single_restart.,
   single_restart_file_out      = .false.,
   write_binary_restart_files   = .true.,
   overwrite_data_time          = .false.,
   new_data_days                = -1,
   new_data_secs                = -1,
   input_is_model_advance_file  = .false.,
   output_is_model_advance_file = .true.,
   overwrite_advance_time       = .true.,
   new_advance_days             = $greg_obs_days,
   new_advance_secs             = $greg_obs_secs,
EOF
   sed -e '/^   input_file_name/,/^   new_advance_secs/d' \
       -e '/&restart_file_tool_nml/r ./input.nml.ic.restart_file_tool' input.nml.filter.${icyc} >! input.nml
 
   $RUN_DIR/restart_file_tool >! logs/restart_file_tool.ic_${icyc}.log

   set n_ics = `ls -1 assim_model_state_ic.* | wc -l`
   if( $n_ics != $nens ) then
       echo Not enough assim_model_state_ic files: $n_ics. Stop.
       exit
   endif

   #------------------------------------------------------
   # 7. Advance model for each member
   #------------------------------------------------------
   # Run forecast for ensemble members until next analysis time
   echo Advance models for $nens members now...
   ${REMOVE} assim_model_state_ud.*

   set n = 1
   while ( $n <= $nens )

   if($bsub_in_ibm == yes) then

      set job_ensemble = ${expname}_${icyc}_ens${n}

      cat >! advance.sed << EOF
      s#JOB_NAME#${job_ensemble}#g
      s#PROJ_NUMBER#${PROJ_NUMBER}#g
      s#ENS_MEM#${n}#g
EOF

      sed -f advance.sed advance_model.template.lsf >! advance_model.lsf
      bsub < advance_model.lsf
      sleep 3

   else

      set icnum = `echo $n + 10000 | bc | cut -b2-5`
      if(-e filter_control${icnum}) ${REMOVE} filter_control${icnum}

      echo $n                            >!   filter_control${icnum}
      echo assim_model_state_ic.${icnum} >>   filter_control${icnum}
      echo assim_model_state_ud.${icnum} >>   filter_control${icnum}
      ./advance_model.csh $n 1 filter_control${icnum}

   endif

   @ n++
   end

   #------------------------------------------------------
   # 8. Store output files
   #------------------------------------------------------
   set sav_dir = $OUT_DIR/${time_anl}
   if(! -d $sav_dir) mkdir -p $sav_dir
   echo Saving output files...
   ls -lrt >! ${sav_dir}/list
   ${COPY} input.nml.filter.${icyc} ${sav_dir}/
   ${COPY} input.nml ${sav_dir}/input.nml.ic.restart_file_tool.${icyc}

   # Save for adaptive inflation for next cycle.
   if( $adaptive_inf == true ) then
     if ( -e $infl_out && ! -z $infl_out ) then
         ${MOVE}  $infl_out ${sav_dir}/
         if ( ! $status == 0 ) then
              echo "Failed moving $infl_out to ${sav_dir}"
              exit
         endif
     else
         echo $infl_out does not exist. Stop.
         exit
     endif
   endif

   foreach FILE ( preassim.nc analysis.nc obs_seq.final )
     ${MOVE} $FILE $sav_dir/
     if( ! $status == 0 ) then
         echo Failed moving $FILE to ${sav_dir}.
         exit
     endif
   end
   ${COPY} $mpas_fname $sav_dir/
   if($icyc == 1) ${MOVE} ${restart_in}* $sav_dir/


   if($bsub_in_ibm == no) then

   ${MOVE} filter.log $sav_dir/

   else

   # Check if all members are done advancing model.
   set is_all_done = `bjobs -w | grep $job_ensemble | wc -l`
   while ( $is_all_done > 0 )
     sleep 60
     set is_all_done = `bjobs -w | grep $job_ensemble | wc -l`
   end
   sleep 60

   endif

   #------------------------------------------------------
   # 9. Check if model was completed for all members.
   #------------------------------------------------------
   cd $RUN_DIR
   set n_prior = `ls -1 assim_model_state_ud.* | wc -l`
   if($n_prior != $nens) then
      set n = 1
      set ifailed = 0
      set ens_failed = ""
      while ( $n <= $nens )
       set iens = `echo $n + 10000 | bc | cut -b2-5`
       if(! -e assim_model_state_ud.${iens}) then
          set ens_failed = `echo $ens_failed $n`
          @ ifailed++
       endif
       @ n++
      end
      echo $ifailed members failed: $ens_failed. Exit.
      exit
   endif

   #------------------------------------------------------
   # 10. Update the mpas template file for the next cycle.
   #------------------------------------------------------
   if( -e $mpas_fname ) ${REMOVE} $mpas_fname
   echo ${COPY} advance_temp1/$mpas_fname .
        ${COPY} advance_temp1/$mpas_fname .

   #------------------------------------------------------
   # 11. Get ready to run filter for next cycle.
   #------------------------------------------------------
   if(! -d $OUT_DIR/${time_nxt}) mkdir -p $OUT_DIR/${time_nxt}

   if($single_restart == true) then

   if( -e $restart_in ) ${REMOVE} $restart_in

   cat >! input.nml.ud.restart_file_tool << EOF
   input_file_name              = "assim_model_state_ud",
   output_file_name             = "$restart_in",
   ens_size                     = $nens,
   single_restart_file_in       = .false.,
   single_restart_file_out      = .$single_restart.,
   write_binary_restart_files   = .true.,
   overwrite_data_time          = .false.,
   new_data_days                = -1,
   new_data_secs                = -1,
   input_is_model_advance_file  = .false.,
   output_is_model_advance_file = .false.,
   overwrite_advance_time       = .false.,
   new_advance_days             = -1,
   new_advance_secs             = -1
EOF
   sed -e '/^   input_file_name/,/^   new_advance_secs/d' \
       -e '/&restart_file_tool_nml/r ./input.nml.ud.restart_file_tool' input.nml.filter.${icyc} >! input.nml

   $RUN_DIR/restart_file_tool >! logs/restart_file_tool.ud_${icyc}.log
   ${COPY} input.nml ${sav_dir}/input.nml.ud.restart_file_tool.${icyc}
   ${MOVE} ${restart_in} $OUT_DIR/${time_nxt}/

   else
 
   set n = 1
   while ( $n <= $nens )
    set iens = `echo $n + 10000 | bc | cut -b2-5`
    ${MOVE} assim_model_state_ud.${iens} $OUT_DIR/${time_nxt}/${restart_in}.${iens}
    @ n++
   end

   endif


   echo Filter is ready to go for the next cycle now.
   set time_anl = $time_nxt
   @ icyc++

end

echo Cycling is done for ${expname}. Script exiting normally.

exit 0


