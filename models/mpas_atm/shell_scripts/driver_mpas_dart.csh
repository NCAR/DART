#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
##############################################################################################
#  driver_mpas_dart.csh
#
#  THIS IS A TOP-LEVEL DRIVER SCRIPT FOR CYCLING RUN. 
#  BOTH THE ENSEMBLE KALMAN FILTER AND THE MPAS FORECAST ARE RUN IN AN MPI VERSION.
#
#  This is a sample script for cycling runs in a retrospective case study, 
#  and was tested on NCAR IBM Supercomputer (yellowstone) using a "bsub" command.
#
#  Note:
#  1. This script does NOT specify all the options available for EnKF data assimilation.
#     For your own complete filter design, you need to edit your template input.nml
#     for the parameters which are not set up in the "Assimilation parameters" section below.
#     You may want to edit at least &filter_nml, &obs_kind_nml, &model_nml, &location_nml 
#     and &mpas_vars_nml sections to confirm your filter configuration before running this script.
#     For adaptive inflation, we only allow the choice for spatialy-varying prior inflation here.
#     For more options, check DART/filter/filter.html and edit this script accordingly.
#  2. An option for back-up in the hpss storage is supported here. 
#     If failed in backup, all the outputs will be locally stored. 
#     For such a case, check if you have enough disk space before running this script.
#     Depending on the sav_ options, you may end up storing huge ensemble data during the cycling.
#  3. This script assumes that the initial ensemble (both in the DART binary format and the mpas 
#     restart format), an mpas restart template file, and obs sequence files for the whole period
#     are available for the cycling run. More descriptions below.
#
#  Required scripts to run this driver:
#  1. namelist.input             (for mpas) - a namelist template for the mpas model run.
#  2. advance_model.csh          (for mpas) - which makes the mpas forecast run during the cycling.
#  3. input.nml                  (for filter) - a namelist template for filter. 
#  4. filter.template.lsf        (for an mpi filter run; with async >= 2)
#  5. advance_model.template.lsf (for an mpi mpas run; using separate nodes for each ensemble member)
#
#  Input files to run this script:
#  A. RUN_DIR/${mpas_fname}.nc   - an mpas template file (for mpas grid fields)
#  B. FG_DIR/${mpas_fname}.e#.nc - initial ensemble in mpas restart format.
#                                  Can be generated from mpas ensemble forecast valid at the initial cycle.
#  C. FG_DIR/${restart_in} - initial ensemble in dart binary format.
#                            Converted from the initial mpas ensemble in B by create_filter_ics.csh.
#  D. OBS_DIR/obs_seq${YYYYMMDDHH} - obs sequence files for each analysis cycle (YYYYMMDDHH) for the 
#                                    whole period (from ${date_ini} to ${date_end}).
#
#  Written by So-Young Ha (MMM/NCAR) Sep-30-2011
#  Updated and tested on yellowstone by So-Young Ha (MMM/NCAR) Feb-20-2013
#
#  For any questions or comments, contact: syha@ucar.edu (+1-303-497-2601)
#
##############################################################################################
# USER SPECIFIED PARAMETERS
##############################################################################################
# Experiment name and the cycle period
#--------------------------------------------------------------------------
set  expname = MPAS_DART_test           # experiment name
set init_dir = 2008080112		# initial date directory for filter_ics
set date_ini = 2008-08-01_12:00:00	# initial cycle for this experiment
set date_beg = 2008-08-01_12:00:00	# start date to run this script
set date_end = 2008-08-31_12:00:00	# end date to run this script
set intv_day = 0			# cycling frequency - assimilation_period_days    in input.nml
set intv_sec = 21600			# cycling frequency - assimilation_period_seconds in input.nml
#--------------------------------------------------------------------------
# Assimilation parameters (Logical parameters should be true or false (case-sensitive.))
#--------------------------------------------------------------------------
set nens            = 80	# ensemble size for ens_size in input.nml
set cutoff          = 0.20	# horizontal location - cutoff in input.nml
set horiz_dist_only = false	# horizontal localization only - if true, edit &location_nml as well.
set adaptive_inf    = true  	# adaptive_inflation - If true, this script only supports 
 				# spatially-varying state space prior inflation.
				# And you also need to edit inf_sd_initial, inf_damping,
				# inf_lower_bound, and inf_sd_lower_bound in &filter_nml.
set single_restart = false	# true if all copies read from and written to a single file in filter.
set use_u_for_wind = false	# Use normal velocity ('u') on edges for wind assimilation
set update_u_from_reconstruct = true  # Use reconstructed winds at cell center, 
                                      # then update normal velocity by the wind increments 
                                      # at cell center (not by filter).
set output_ens_obs = true	# Print out ensemble observations in the output obs sequence file.
#--------------------------------------------------------------------------
# Configuration for MPAS cycling (for advance_model.csh and namelist.input)
#--------------------------------------------------------------------------
set   sav_analysis = false	# true if you want to save the ensemble analysis locally (in ENS_DIR).
set   sav_forecast = false	# true if you want to save the ensemble forecast locally (in ENS_DIR).
set   sav_diag_prs = false	# true if you want to save diagnostic variables at pressure levels (in ENS_DIR).
set     sst_update = true       # true if config_sst_update = true in the model simulation
set   sst_interval = 01_00:00:00	# sst update interval when sst_update = true
set   sst_fname = sfc_update.nc		# sst input file name when sst_update = true
#--------------------------------------------------------------------------
# Directories
#--------------------------------------------------------------------------
set DART_DIR = /glade/scratch/syha/DART/branch_dev/models/mpas_atm/work # where dart executables exist
set MPAS_DIR = /glade/scratch/syha/MPAS/run      	# where all the aux files exist to run mpas forecast
set  RUN_DIR = /glade/scratch/syha/MPAS_DART/$expname	# top-level working directory
set   FG_DIR = $RUN_DIR/$init_dir                       # where filter_ics exists for initial cycle
set  OBS_DIR = /glade/scratch/syha/OBS_SEQ/data        	# where all obs sequence files exist
set  ENS_DIR = advance_temp          # where the background forecast exists for each member 
#--------------------------------------------------------------------------
# File naming convention (for input.nml)
#--------------------------------------------------------------------------
set restart_in  = filter_ics
set restart_out = filter_restart
set    infl_in  = prior_inflate_ics
set    infl_out = prior_inflate_restart
set  obs_seq_in = obs_seq.out
set obs_seq_out = obs_seq.final
set  mpas_fname = mpas_init	# both for the analysis and the grid files (without .nc)
#--------------------------------------------------------------------------
# BSUB setup on NCAR yellowstone (only if $bsub_in_ibm = yes)
#--------------------------------------------------------------------------
set  bsub_in_ibm = yes  	# run on yellowstone? yes or no.
set  PROJ_NUMBER = P64000101	# Account key to submit filter and mpas jobs on yellowstone
set  time_filter = 01:10	# wall clock time for mpi filter runs
set        queue = small	# quene for the batch job

set    hpss_save = yes          # Backup in HPSS? yes or no. If yes, edit below.
set       hsidir = MPAS/CYCLE_TEST/$expname
set       HSICMD = 'hsi put -P'
####################################################################################
# END OF USER SPECIFIED PARAMETERS
####################################################################################

set REMOVE = '/bin/rm -rf'
set   COPY = 'cp -pf'
set   MOVE = 'mv -f'
set   LINK = 'ln -sf'
unalias cd
unalias ls

if(! -e $RUN_DIR) mkdir -p $RUN_DIR
cd $RUN_DIR

#------------------------------------------
# Check if we have all the necessary files.
#------------------------------------------
if(! -r input.nml) ${COPY} $DART_DIR/input.nml .
if(! -r namelist.input) ${COPY} $MPAS_DIR/namelist.input .

if ( ! -e advance_model.csh) then
    echo ${COPY} $DART_DIR/../shell_scripts/advance_model.csh .
         ${COPY} $DART_DIR/../shell_scripts/advance_model.csh .
    if( ! $status == 0 ) then
        echo ABORT\: We cannot find required script $fn.
    endif
endif

# Update advance_model.csh
cat >! config_restart.sed << EOF
/set save_analysis/c\
set save_analysis = ${sav_analysis}
/set save_forecast /c\
set save_forecast = ${sav_forecast}
/set save_diag_plev /c\
set save_diag_plev = ${sav_diag_prs}
EOF
if(-e advance_new.csh) ${REMOVE} advance_new.csh
sed -f config_restart.sed advance_model.csh >! advance_new.csh
${MOVE} advance_new.csh advance_model.csh
chmod +x advance_model.csh

# For an initial cycle, we already linked the mpas template files 
# thru create_filter_ics.csh. But if you run this script in the middle
# of the cycling, you need to make sure that a correct background
# forecast is linked in each ensemble member directory.

set ie = 1
while ( $ie <= $nens )
if($date_ini == $date_beg) then
  if(! -e ${FG_DIR}/${mpas_fname}.e${ie}.nc) then
     echo We cannot find an initial ensemble for member ${ie}. Stop.
     exit -1
  endif
  ${LINK} ${FG_DIR}/${mpas_fname}.e${ie}.nc
else
  set mpas_temp = ${ENS_DIR}${ie}/${mpas_fname}.nc
  if(! -e $mpas_temp) then
     echo We cannot find ${mpas_temp} to start with. Stop.
     exit -1
  endif
  set time_ie = `ncdump -v xtime ${mpas_temp} | tail -2 | head -1 | awk '{print $1}' | cut -c2-`
  if($time_ie != $date_beg) then
     echo $mpas_temp should have the time $date_beg, but has $time_ie. Stop.
     exit -1
  endif
endif
@ ie++
end

if(! -e ${mpas_fname}.nc) ${COPY} ${FG_DIR}/${mpas_fname}.e1.nc ${mpas_fname}.nc


foreach fn ( filter advance_time convertdate restart_file_tool dart_to_model model_to_dart )
  if ( ! -x $fn ) then
       echo ${LINK} $DART_DIR/$fn .
            ${LINK} $DART_DIR/$fn .
       if( ! $status == 0 ) then
           echo ABORT\: We cannot find required executable dependency $fn.
           exit -1
       endif
  endif
end 

if($bsub_in_ibm == yes) then
 foreach fn ( filter.template.lsf advance_model.template.lsf )
  if ( ! -r $fn ) then
      echo ${COPY} $DART_DIR/../shell_scripts/$fn .
           ${COPY} $DART_DIR/../shell_scripts/$fn .
      if( ! $status == 0 ) then
          echo ABORT\: We cannot find required readable dependency $fn.
          exit -1
      endif
  endif 
 end
endif

# Preparation for the model run
if ( -e MPAS_RUN ) ${REMOVE} MPAS_RUN
if (! -d $MPAS_DIR ) then
    echo $MPAS_DIR does not exist. Stop.
    exit -1
endif
${LINK} $MPAS_DIR MPAS_RUN

if ( $sst_update == true ) then
     cat >! namelist.sed << EOF
  /config_sfc_update_interval /c\
   config_sfc_update_interval = '${sst_interval}'
  /config_sfc_update_name /c\
   config_sfc_update_name    = '${sst_fname}'
/config_sst_update /c\
 config_sst_update           = .true.
EOF

sed -f namelist.sed namelist.input >! namelist.input.sst_update
if(! -z namelist.input.sst_update) then
   ${MOVE} namelist.input.sst_update namelist.input
else
   echo ABORT\: Failed in updating namelist.input for config_sst_update.
   exit -1
endif
${LINK} ${MPAS_DIR}/${sst_fname} .		|| exit 1

else
   echo NO SST_UPDATE...
endif

if ( ! -d logs ) mkdir logs			# to print out log files
${COPY} input.nml input.nml.template 		# Need to update input.nml with user-specified options

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
if($n_cycles < 0) then
   echo Cannot figure out how many cycles to run. Check the time setup.
   exit
endif

echo Running at $RUN_DIR
echo " "

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
# Cycling gets started
#--------------------------------------------------------
set time_ini = `echo $date_ini 0 | advance_time`		#YYYYMMDDHH
set time_anl = `echo $date_beg 0 | advance_time`		#YYYYMMDDHH
set time_end = `echo $date_end 0 | advance_time`		#YYYYMMDDHH

set first_cycle = $icyc
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
   /advance_time_present/c\
   advance_time_present     = .true.,
EOF

  if( $adaptive_inf == true ) then       # For a spatially-varying prior inflation.

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

  set is_update_u_there = `grep update_u_from_reconstruct input.nml.template | wc -l`

  if($is_update_u_there == 0) then
     echo No update_u_from_reconstruct in your input.nml.template.
     exit -1
  endif

  cat >! scriptu.sed << EOF
  /update_u_from_reconstruct/c\
   update_u_from_reconstruct    = .${update_u_from_reconstruct}.,
EOF
  cat scriptu.sed >> script.sed

  set nobs = 0
  if( $output_ens_obs == true ) set nobs = $nens
  cat >! script4.sed << EOF
  /num_output_obs_members/c\
   num_output_obs_members   = $nobs,
EOF
  cat script4.sed >> script.sed

  ${REMOVE} input.nml
  sed -f script.sed input.nml.template >! input.nml 		|| exit 2

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
   # 3. Obs sequence for this analysis cycle - one obs time at each analysis cycle
   #------------------------------------------------------
   #set fn_obs = $OBS_DIR/${time_anl}/obs_seq.out
   set fn_obs = $OBS_DIR/obs_seq${time_anl}
   if(! -e $fn_obs) then
      echo $fn_obs does not exist. Stop.
      exit
   endif
   ${LINK} $fn_obs ${obs_seq_in}

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
   s#QUEUE#${queue}#g
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
      s#QUEUE#${queue}#g
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
   set sav_dir = $RUN_DIR/${time_anl}
   if(! -d ${sav_dir}) mkdir -p ${sav_dir}
   echo Saving output files...
   echo ${sav_dir} >! ${sav_dir}/list
   ls -lrt >> ${sav_dir}/list
   ls -al MPAS_RUN/*.exe >> ${sav_dir}/list

   ${COPY} input.nml.filter.${icyc} ${sav_dir}/
   ${COPY} input.nml ${sav_dir}/input.nml.ic.restart_file_tool.${icyc}

   # Check HPSS first.
   set hd = ${hsidir}/${time_anl}
   if ( ${hpss_save} == yes ) then
        hsi mkdir -p ${hd}
        if ( ! $status == 0 ) then
        echo We cannot access to ${hd}. We back up files locally.
        set hpss_save = no
        endif
   endif
   echo "HSI BACKUP? ${hpss_save} in ${hd}/"

   # We can back up these files in
   if ( $hpss_save == yes ) then
        #set hd = $hsidir/${time_anl}
        #echo Backup to HPSS:$hd ...

        if ( $icyc == 1) then
             ${COPY} ${sav_dir}/list .
             tar cvf scripts.tar *.csh *template* input.nml namelist.input list
             ${HSICMD} -d scripts.tar : ${hsidir}/scripts.tar
        endif

        foreach FILE ( input.nml.filter.${icyc} list )
          echo ${HSICMD} ${sav_dir}/$FILE : ${hd}/$FILE
               ${HSICMD} ${sav_dir}/$FILE : ${hd}/$FILE
        end
 
        if( $single_restart == false ) then
            set frst = ${restart_in}.tar
            echo tar cfh ${frst} ${restart_in} files...
                 tar cfh ${frst} ${restart_in}*
        else
            set frst = ${restart_in}
        endif
        gzip -f $frst
            
        foreach FILE ( Prior_Diag.nc Posterior_Diag.nc $obs_seq_out ${frst}.gz )
          echo ${HSICMD} $FILE : ${hd}/$FILE
               ${HSICMD} $FILE : ${hd}/$FILE
          if ( ! $status == 0 ) then
               echo "Failed in back-up $FILE in $hd/"
               echo ${MOVE} $FILE ${sav_dir}/ instead.
               ${MOVE} $FILE ${sav_dir}/
               touch BOMBED
          endif
        end 
        ${REMOVE} ${frst}.gz

        # We need to keep obs_seq.final for post-processing later.
        ${MOVE} $obs_seq_out ${sav_dir}/
         
        if( $adaptive_inf == true ) then
          if ( -e $infl_out && ! -z $infl_out ) then
               echo ${HSICMD} $infl_out : ${hd}/$infl_out
                    ${HSICMD} $infl_out : ${hd}/$infl_out
          else
              echo $infl_out does not exist. Stop.
              exit
          endif
        endif


   else		# local back-up

   foreach FILE ( Prior_Diag.nc Posterior_Diag.nc $obs_seq_out )
     ${MOVE} $FILE $sav_dir/
     if( ! $status == 0 ) then
         echo Failed moving $FILE to ${sav_dir}.
         exit
     endif
   end

   endif	# hpss_save

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
   set ifailed = 0
   set ic_bsub = 0
   if($n_prior != $nens) then
      set n = 1
      set ens_failed = ""
      while ( $n <= $nens )
       set iens = `echo $n + 10000 | bc | cut -b2-5`
          
       if( -e assim_model_state_ic.${iens}) then
           set job_ensemble = ${expname}_${icyc}_ens${n}
           cat >! advance.sed << EOF
           s#JOB_NAME#${job_ensemble}#g
           s#ENS_MEM#${n}#g
EOF
           sed -f advance.sed advance_model.template.lsf >! advance_model.lsf
           bsub < advance_model.lsf
           sleep 10
           @ ic_bsub++
       endif

       if(! -e assim_model_state_ud.${iens}) then
          set ens_failed = `echo $ens_failed $n`
          @ ifailed++
       endif
       @ n++
      end

      if($ic_bsub > 0) then
       set is_it_done = `bjobs -w | grep ${expname}_${icyc} | wc -l`
       while ( $is_it_done > 0 )
        sleep 60
        set is_it_done = `bjobs -w | grep ${expname}_${icyc} | wc -l`
       end
       sleep 60
      else
       echo $ifailed members failed: $ens_failed. Exit.
       exit
      endif

   endif

   #------------------------------------------------------
   # 10. Get ready to run filter for next cycle.
   #------------------------------------------------------
   if(! -d $RUN_DIR/${time_nxt}) mkdir -p $RUN_DIR/${time_nxt}

   if($single_restart == true) then

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
   ${MOVE} ${restart_in} $RUN_DIR/${time_nxt}/

   else
 
   set n = 1
   while ( $n <= $nens )
    set iens = `echo $n + 10000 | bc | cut -b2-5`
    ${MOVE} assim_model_state_ud.${iens} $RUN_DIR/${time_nxt}/${restart_in}.${iens}
    @ n++
   end

   endif

   if( ! $status == 0 ) then
   echo Filter is ready to go for the next cycle now.
   endif

   set time_anl = $time_nxt
   @ icyc++

end

echo Final backup...
set sav_dir = ${time_anl}
if( $single_restart == false ) then
    set frst = ${restart_in}.tar
    echo tar cfh ${frst} ${sav_dir}/${restart_in} files...
         tar cfh ${frst} ${sav_dir}/${restart_in}*
else
    set frst = ${restart_in}
endif
gzip -f $frst

gzip -f ${ENS_DIR}*/${mpas_fname}.nc
tar cvf FG.${time_anl}.tar ${ENS_DIR}*/${mpas_fname}.nc.gz

if ( $hpss_save == yes ) then
     set hd = ${hsidir}/${time_anl}
     hsi mkdir -p ${hd}
     ${HSICMD} -d ${frst}.gz ${hd}/${frst}.gz
     ${HSICMD} -d FG.${time_anl}.tar ${hd}/FG.${time_anl}.tar  
else
     ${MOVE} ${frst}.gz $sav_dir/
     ${MOVE} FG.${time_anl}.tar $sav_dir/
endif

echo Cycling is done for $n_cycles cycles in ${expname}.
echo Last ensemble forecasts are valid at ${time_anl}.
echo Script exiting normally.

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

