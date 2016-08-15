#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
##############################################################################################
#
#  driver_filter_mpas.csh
#
#  THIS IS A TOP-LEVEL DRIVER SCRIPT TO RUN FILTER AND THE MPAS MODEL BOTH IN AN MPI VERSION.
#
#  This is a sample script for a retrospective case study, and was tested on
#  NCAR IBM Supercomputer (yellowstone and yellowstone) using a "bsub" command.
#  Any back-up option to the hpss storage is not supported here. 
#  Instead, all the outputs will be locally stored. 
#  Check your disk space before running this script.
#
#  Required file to run this script:
#  1. namelist.atmosphere     (for mpas) - a namelist template for the mpas model run.
#  2. advance_model.csh       (for mpas) - which decides how to update the model states.
#  3. mpas_init.nc            (for mpas) - you can change the filename below (mpas_fname).
#  4. input.nml               (for filter) - a namelist template for filter. 
#  5. filter_ics              (for filter) - an initial filter file to start with.
#  6. filter.template.lsf     (for an mpi filter run; with async >= 2)
#  7. advance_model.template.lsf (for an mpi mpas run; using separate nodes)
#
#  Note1: For your own complete filter design, you need to edit your input.nml
#         for the parameters which are not set up in "Assimilation parameters" below.
#         At least double-check your &filter_nml, &obs_kind_nml, &model_nml, &location_nml 
#         and &mpas_vars_nml before running this script.
#  Note2: The mpas model configurations are not edited in this script.
#         Thus, namelist.atmosphere is assumed to be all set up before being called here.
#  Note3: advance_model.csh takes all the necessary files for the model run from MPAS_RUN/
#         except for an initial condition (e.g., analysis from DART).
#  Note4: All the logical parameters are case-sensitive. They should be either true or false.
# 
#  Note5: Before running this script, we assume that
#         - both mpas and dart are successfully compiled.
#         - Initial ensemble forecasts are ready in both mpas netcdf format and DART binary format in IENS_DIR/
#         - All the observations are ready in OBS_DIR/ for all cycles.
#
#  Input files to run filter:
#  1. FG_DIR/${restart_in} - first guess in dart binary format for the initial cycle (ex. filter_ics)                           
#  2. OBS_DIR/YYYYMMDDHH/${obs_seq_in} - an obs sequence file for each analysis cycle.
#
#  Output from this script will be all stored in $RUN_DIR/YYYYMMDDHH/ for each analysis cycle.
#
#  Written by Soyoung Ha (MMM/NCAR) Sep-30-2011
#  Updated by Soyoung Ha (MMM/NCAR) Apr-16-2013 for yellowstone
##############################################################################################
# Experiment name and the cycle period
#--------------------------------------------------------------------------
set echo

set  expname = OSSE			# experiment name
set init_cyc = 2008080112               # initial cycle (YYYYMMDDHH)
set date_ini = 2008-08-01_12:00:00	# initial cycle for this experiment in UTC format
set date_beg = 2008-08-01_12:00:00	# start date to run this script in UTC format
set date_end = 2008-08-31_00:00:00	# end date to run this script in UTC format
set intv_day = 0			# cycling frequency - assimilation_period_days    in input.nml
set intv_sec = 21600			# cycling frequency - assimilation_period_seconds in input.nml
#--------------------------------------------------------------------------
# Assimilation parameters (Logical parameters should be true or false (case-sensitive.))
#--------------------------------------------------------------------------
set nens            = 96	# ensemble size for ens_size in input.nml
set cutoff          = 0.20	# horizontal location - cutoff in input.nml
set horiz_dist_only = true	# horizontal localization only (true or false)
set adaptive_inf    = false 	# adaptive_inflation - If true, this script only supports 
 				# spatially-varying state space prior inflation.
				# And you also need to edit inf_sd_initial, inf_damping,
				# inf_lower_bound, and inf_sd_lower_bound in &filter_nml.
set single_restart = false	# true if all copies read from and written to a single file in filter.
set use_u_for_wind = false	# Use normal velocity ('u') on edges for wind assimilation
set output_ens_obs = true	# Print out ensemble observations in the output obs sequence file.
#--------------------------------------------------------------------------
# Configuration for cycling
#--------------------------------------------------------------------------
set   sav_ens_anal = true	# true if you want to save the ensemble analysis locally (in ENS_DIR).
set   sav_ens_fcst = false	# true if you want to save the ensemble forecast locally (in ENS_DIR).
#--------------------------------------------------------------------------
# Directories
#--------------------------------------------------------------------------
set DART_DIR = /mmm/mmmtmp/syha/MPAS/DART/models/mpas_atm/work	# where dart executables exist
set MPAS_DIR = /mmm/mmmtmp/syha/MPAS             		# where all the auxiliary files exist to run mpas
set  OBS_DIR = /mmm/mmmtmp/syha/OBS_SEQ				# where all obs sequence files exist
set  RUN_DIR = `pwd`						# top-level working directory
set   FG_DIR = $RUN_DIR/${init_cyc}                             # where filter_ics exists for initial cycle (2008-08-01_12 UTC)
set  OUT_DIR = $RUN_DIR      				        # storage for output - a subdirectory for each cycle will be
                                                                # automatically generated during cycles
set  ENS_DIR = advance_temp                                     # each ensemble member directory (under $RUN_DIR)
set IENS_DIR = /mmm/mmmtmp/syha/INITIAL_ENSEMBLE                # initial ensemble directory (for the initial cycle)
#--------------------------------------------------------------------------
# File naming convention - you may not need to change dart file names.
#--------------------------------------------------------------------------
set restart_in  = filter_ics
set restart_out = filter_restart
set    infl_in  = prior_inflate_ics
set    infl_out = prior_inflate_restart
set obs_seq_in  = obs_seq.out
set obs_seq_out = obs_seq.final

set    mpas_exe = atmosphere_model      # mpas model executable
set    mpas_grd = x1.10242              # mpas grid info    
set   mpas_fini = mpas_init	        # can be used for both the analysis and the grid info files (without .nc)
set   mpas_frst = x1.10242.restart	# prefix of a template file (in ENS_DIR for the initial cycle)
					# should be matched with either config_output_name or 
					# config_restart_name in namelist.atmosphere (without .nc)

#--------------------------------------------------------------------------
# LSF setup on NCAR yellowstone (only if $bsub_in_ibm = yes)
#--------------------------------------------------------------------------
set  bsub_in_ibm = yes  	# run on yellowstone? yes or no.
set  PROJ_NUMBER = NMMM000x	# Account key to submit filter and mpas jobs in yellowstone
set  time_filter = 00:20	# wall clock time for mpi filter runs

####################################################################################
# End of User Specificaton
####################################################################################

set REMOVE = '/bin/rm -rf'
set   COPY = 'cp -pf'
set   MOVE = 'mv -f'
set   LINK = 'ln -sf'
unalias cd
unalias ls

echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'    
echo Running $0 at $RUN_DIR
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'  
echo We edit input.nml for DART thru this script,
echo but the MPAS model configuration is assumed to be all setup thru
echo namelist.atmosphere and stream files (except for the cycling part).    
echo
cat namelist.atmosphere
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'  
echo

#------------------------------------------
# Check if we have all the necessary files.
#------------------------------------------
foreach fn ( input.nml namelist.atmosphere )
  if ( ! -r $fn ) then
      echo ABORT\: We cannot find required readable dependency $fn.
  endif 
end

if ( ! -e advance_model.csh) then
    echo ${COPY} $DART_DIR/../shell_scripts/advance_model.csh .
         ${COPY} $DART_DIR/../shell_scripts/advance_model.csh .
endif

# Update advance_model.csh
cat >! config_restart.sed << EOF
/set save_forecast/c\
set save_forecast = ${sav_ens_fcst}
/set save_analysis/c\
set save_analysis = ${sav_ens_anal}
EOF
if(-e advance_new.csh) ${REMOVE} advance_new.csh
sed -f config_restart.sed advance_model.csh >! advance_new.csh
${MOVE} advance_new.csh advance_model.csh
chmod +x advance_model.csh


#------------------------------------------
# Initialize the first cycle
    
set ie = 1
while ( $ie <= $nens )
if($date_ini == $date_beg) then

  # Get all the initial ensemble files ready.
  if( $ie == 1 && ! -d ${FG_DIR}) mkdir -p ${FG_DIR}
  set icnum = `echo $ie + 10000 | bc | cut -b2-5`
  set mpas_fname = ${mpas_fini}
  if(! -e ${FG_DIR}/${mpas_fname}.e${ie}.nc) ${LINK} ${IENS_DIR}/${mpas_fname}.e${ie}.nc ${FG_DIR}/	|| exit
  if(! -e ${FG_DIR}/${restart_in}.${icnum})  ${LINK} ${IENS_DIR}/${restart_in}.${icnum} ${FG_DIR}/	|| exit
  ${LINK} ${FG_DIR}/${mpas_fname}.e${ie}.nc .		|| exit
  ${LINK} ${FG_DIR}/${restart_in}.${icnum} .      	|| exit
  if(-d ${ENS_DIR}1) ${REMOVE} ${ENS_DIR}*

else

  echo "Skip initializing ensemble..."
  set mpas_temp = ${ENS_DIR}${ie}/${mpas_fname}.nc
  set time_ie = `ncdump -v xtime ${mpas_temp} | tail -2 | head -1 | awk '{print $1}' | cut -c2-`
  if($time_ie != $date_beg) then
     echo $mpas_temp should have the time $date_beg, but has $time_ie. Stop.
     exit -1
  endif

endif
@ ie++
end

foreach fn ( filter advance_time convertdate restart_file_tool dart_to_model model_to_dart )
  if ( ! -x $fn ) then
       echo ${LINK} $DART_DIR/$fn .
            ${LINK} $DART_DIR/$fn .
       #echo ABORT\: We cannot find required executable dependency $fn.
       #exit -1
  endif
end 

if($bsub_in_ibm == yes) then
 foreach fn ( filter.template.lsf advance_model.template.lsf )
  if ( ! -r $fn ) then
      echo ${COPY} $DART_DIR/../shell_scripts/$fn .
           ${COPY} $DART_DIR/../shell_scripts/$fn .
      #echo ABORT\: We cannot find required readable dependency $fn.
      #exit -1
  endif 
 end
endif

cd $MPAS_DIR    
set mfiles = ( ${mpas_exe} ${mpas_grd}.graph.info* *BL *DATA stream* )
cd -

if ( ! -d MPAS_RUN ) ${LINK} $MPAS_DIR MPAS_RUN		# for advance_model.csh
if ( ! -d logs ) mkdir logs				# to print out log files
${COPY} input.nml input.nml.template 			# for advance_time

foreach f ( $mfiles )    
    ${LINK} MPAS_RUN/$f .
end
      
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
else
   if( -d advance_temp1) then
       echo We start new experiment now. 
       echo Cleaning up advance_temp directories first...
      \rm -fR advance_temp*
   endif
endif
set ncyc = `expr $icyc \+ $n_cycles \- 1`   


#--------------------------------------------------------
# MAIN LOOP FOR CYCLING
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
    
  # MPAS namelist configuration
  #--------------------------------------------------------
  if($icyc == 1) then      # First cycle is run as cold start.
     set mpas_fname = ${mpas_fini}
     set cycling    = .false.
     set do_restart = .false.
     if(! -e ${mpas_fname}.nc) ${COPY} ${FG_DIR}/${mpas_fname}.e1.nc ${mpas_fname}.nc  || exit 1
  else
     set mpas_fname = ${mpas_frst}
     set cycling    = .true.
     set do_restart = .true.
     set mpas_template = ${RUN_DIR}/advance_temp1/${mpas_fname}.nc
    #set mpas_template = ${IENS_DIR}/ENS_1/${mpas_fname}.`echo $date_ini | sed -e 's/:/./g'`.nc
     if(! -e ${mpas_fname}.nc) ln -sf $mpas_template ${mpas_fname}.nc		|| exit
  endif

  ${REMOVE} init.sed  
  cat >! init.sed << EOF
  /config_do_DAcycling /c\
   config_do_DAcycling = $cycling
  /config_do_restart /c\
   config_do_restart = ${do_restart}
EOF
  mv namelist.atmosphere namelist.input.temp                    || exit
  sed -f init.sed namelist.input.temp >! namelist.atmosphere    || exit
            
  # DART namelist configuration
  #--------------------------------------------------------    
   cat >! script.sed << EOF
  /ens_size /c\
   ens_size                    = $nens,
  /cutoff  /c\
   cutoff                      = $cutoff,
  /horiz_dist_only/c\
   horiz_dist_only             = .${horiz_dist_only}.,
  /start_from_restart/c\
   start_from_restart       = .true.,
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
   model_analysis_filename      = '${mpas_fname}.nc',
   /grid_definition_filename/c\
   grid_definition_filename     = '${mpas_fname}.nc',
   /output_state_vector /c\
   output_state_vector          = .false.,
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

  if( $use_u_for_wind == true ) then
      set is_use_u_there = `grep use_u_for_wind input.nml.template | wc -l`
      set is_u_there = `grep KIND_EDGE_NORMAL_SPEED input.nml.template | wc -l`
	
      if($is_use_u_there == 0) then
         echo No use_u_for_wind in your input.nml.template.
         exit -1
      endif
      if($is_u_there == 0) then
         echo No KIND_EDGE_NORMAL_SPEED in your input.nml.template.
         exit -1
      endif

      cat >! script3.sed << EOF
  /use_u_for_wind/c\
   use_u_for_wind               = .true.,
EOF

      cat script3.sed >> script.sed
  endif

  set nobs = 0
  if( $output_ens_obs == true ) set nobs = $nens
  cat >! script4.sed << EOF
  /num_output_obs_members/c\
   num_output_obs_members   = $nobs,
EOF
  cat script4.sed >> script.sed

  sed -f script.sed input.nml.template >! input.nml

   #------------------------------------------------------
   # 2. Link to a restart file to get filter started 
   # (assuming start_from_restart = .true. in input.nml)
   #------------------------------------------------------

   if($icyc == 1) then
      set dir_rst = $FG_DIR
   else
      set dir_rst = $RUN_DIR/${time_anl}
   endif
   set fn_rst = $dir_rst/${restart_in}
   if(! -e ${mpas_fname}.nc) ${LINK} $dir_rst/${mpas_fname}.nc .

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
       if(! -e $RUN_DIR/${time_pre}/${infl_out}) then
          echo $RUN_DIR/${time_pre}/${infl_out} does not exist. Stop.
          exit
       endif
       ${LINK} $RUN_DIR/${time_pre}/${infl_out} ${infl_in}
   endif

   #------------------------------------------------------
   # 3. Obs sequence for this analysis cycle 
   #------------------------------------------------------
   #set fn_obs = $OBS_DIR/${time_anl}/obs_seq.out
   set fn_obs = $OBS_DIR/obs_seq.${time_anl}.out
   if(! -e $fn_obs) then
      echo $fn_obs does not exist. Stop.
      exit
   endif
   ${LINK} $fn_obs ${obs_seq_in}

   #------------------------------------------------------
   # 4. Run filter
   #------------------------------------------------------
   set wdir = `pwd`
   set job_name = ${expname}.${icyc}
   echo Running filter: $job_name at $wdir

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
 
   $DART_DIR/restart_file_tool >! logs/restart_file_tool.ic_${icyc}.log

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
   echo Saving output files in ${sav_dir}...
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

   foreach FILE ( Prior_Diag.nc Posterior_Diag.nc obs_seq.final )
     ${MOVE} $FILE $sav_dir/
     if( ! $status == 0 ) then
         echo Failed moving $FILE to ${sav_dir}.
         exit
     endif
   end
   ${COPY} ${mpas_fname}.nc $sav_dir/
   if($icyc != 1) ${MOVE} ${restart_in}* $sav_dir/


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
   # 10. Get ready to run filter for the next cycle
   #------------------------------------------------------
   if(! -d $OUT_DIR/${time_nxt}) mkdir -p $OUT_DIR/${time_nxt}

   if($single_restart == true) then

   #if( -e $restart_in ) ${REMOVE} $restart_in

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

   $DART_DIR/restart_file_tool >! logs/restart_file_tool.ud_${icyc}.log
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

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

