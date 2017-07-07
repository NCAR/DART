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
#  THIS IS A TOP-LEVEL DRIVER SCRIPT FOR CYCLING RUNS. 
#  BOTH THE ENSEMBLE KALMAN FILTER AND THE MPAS FORECAST ARE RUN IN MPI.
#
#  This is a sample script for cycling runs in a retrospective case study, 
#  and was tested on the NCAR IBM Supercomputer (cheyenne) using a "qsub" command.
#
#  Note:
#  1. For the general configuration including all the file names, edit setup_parameters.csh.
#  2. For the model configuration, our general policy is that we only edit the parameters that
#     affect the I/O stream here and leave all the rest unchanged (ex. physics options). 
#     This means that it is the user's responsibility to edit all other namelist parameters 
#     before running this script. One exception is the time info, which will be updated inside 
#     advance_model.csh for each cycle.
#  3. This script does NOT specify all the options available for EnKF data assimilation either.
#     For your own complete filter design, you need to edit your input.nml
#     - at least &filter_nml, &obs_kind_nml, &model_nml, &location_nml and &mpas_vars_nml sections 
#     to set up your filter configuration before running this script.
#  4. For adaptive inflation, we only support the choice of prior adaptive inflation in the state
#     space here. For more options, check DART/assimilation_code/modules/assimilation/filter_mod.html.
#  5. All the logical parameters are case-sensitive. They should be either true or false.
#  6. An option for back-up in the hpss storage is supported here. 
#     If failed in backup, all the output files will be locally stored. 
#     For such a case, check if you have enough disk space before running this script.
#
#  Required scripts to run this driver:
#  (All the template files are available in DART/models/mpas_atm/shell_scripts/.)
#  1. setup_params.csh           (for the general configuration of this experiment)
#  2. namelist.atmosphere        (for mpas)   - a namelist template for mpas.
#  3. input.nml                  (for filter) - a namelist template for filter. 
#  4. filter.template.pbs        (for an mpi filter run; with async >= 2)
#  5. advance_model.template     (for an mpi mpas run; using separate nodes for each ensemble member)
#  6. advance_model.csh          (for mpas/filter) - a driver script to run mpas forecast at each cycle
#
#  Input files to run this script:
#  A. input_state_file_list  - a list of input ensemble netcdf files for DART/filter
#  B. output_state_file_list - a list of output ensemble netcdf files from DART/filter
#  C. RUN_DIR/member#/${mpas_filename}    - the input file listed in input_state_file_list for each member
#  D. OBS_DIR/${obs_seq_in}.${YYYYMMDDHH} - obs sequence files for each analysis cycle (YYYYMMDDHH) 
#     for the entire period.
# 
#  Written by Soyoung Ha (MMM/NCAR)
#  Updated and tested on yellowstone (Feb-20-2013)
#  Updated for MPAS V5 and DART/Manhattan; tested on cheyenne (Jun-27-2017)
#  Updated for a better streamline and consistency: Ryan Torn (Jul-6-2017)
#
#  For any questions or comments, contact: syha@ucar.edu (+1-303-497-2601)
##############################################################################################
# USER SPECIFIED PARAMETERS
##############################################################################################

set echo

if ( $#argv >= 1 ) then
   set fn_param = ${1}
else
   set fn_param = `pwd`/setup_params.csh
endif

if(! -e $fn_param ) then
   echo $fn_param does not exist. Cannot proceed.
   exit
endif

source $fn_param

#--------------------------------------------------------------------------
# Experiment name and the cycle period
#--------------------------------------------------------------------------
echo Experiment name: $EXPERIMENT_NAME

####################################################################################
# END OF USER SPECIFIED PARAMETERS
####################################################################################

if( ! -e $RUN_DIR ) mkdir -p $RUN_DIR
cd $RUN_DIR
echo Running at $RUN_DIR
${COPY} ${fn_param} .

#------------------------------------------
# Check if we have all the necessary files.
#------------------------------------------

foreach fn ( advance_model.csh )
   if ( ! -r $fn || -z $fn ) then
      echo ${COPY} ${DART_DIR}/../shell_scripts/${fn} .
           ${COPY} ${DART_DIR}/../shell_scripts/${fn} .
      if( ! $status == 0 ) then
         echo ABORT\: We cannot find required script $fn
         exit
      endif
   endif
end

foreach fn ( filter advance_time update_mpas_states )
   if ( ! -x $fn ) then
      echo ${COPY} $DART_DIR/${fn} .
           ${COPY} $DART_DIR/${fn} .
      if ( ! $status == 0 ) then
         echo ABORT\: We cannot find required executable dependency $fn.
         exit
      endif
   endif
end 

if( $RUN_IN_PBS == yes ) then
   foreach fn ( filter.template.pbs advance_model.template )
      if ( ! -r $fn ) then
         echo ${COPY} ${DART_DIR}/../shell_scripts/${fn} .
              ${COPY} ${DART_DIR}/../shell_scripts/${fn} .
         if ( ! $status == 0 ) then
            echo ABORT\: We cannot find required script $fn.
            exit
         endif
      endif 
   end
endif

if ( ! -d MPAS_RUN ) then

   if ( ! -d $MPAS_DIR ) then
      echo $MPAS_DIR does not exist. Stop.
      exit
   endif
   ${LINK} $MPAS_DIR MPAS_RUN

endif

#  Check to see if MPAS and DART namelists exist.  If not, copy them from model source
foreach fn ( ${NML_MPAS} ${NML_INIT} )
  if ( ! -r $fn ) then
    ${COPY} ${MPAS_DIR}/${fn} .
  endif
end

foreach fn ( ${STREAM_ATM} ${STREAM_INIT} )
  if ( ! -r $fn || -z $fn ) then
    ${COPY} ${MPAS_DIR}/${fn} .
  endif
end

if ( ! -r ${NML_DART} ) then
  ${COPY} ${DART_DIR}/${NML_DART} .
endif


#--------------------------------------------------------------------------
#  Take file names from input.nml, check to make sure there is consistency in variables.
#--------------------------------------------------------------------------
set  input_list = `grep input_state_file_list  ${NML_DART} | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g" | sed -e 's/"//g'`
set output_list = `grep output_state_file_list ${NML_DART} | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g" | sed -e 's/"//g'`
set  obs_seq_in = `grep obs_sequence_in_name   ${NML_DART} | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g" | sed -e 's/"//g'`
set obs_seq_out = `grep obs_sequence_out_name  ${NML_DART} | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g" | sed -e 's/"//g'`

set input_var = `grep ens_size ${NML_DART} | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g" | sed -e 's/"//g'`
if ( $input_var != $ENS_SIZE ) then
  echo "ERROR ens_size in ${NML_DART} does not match script ${ENS_SIZE}.  Exiting"
  exit
endif

set input_var = `grep assimilation_period_days ${NML_DART} | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g" | sed -e 's/"//g'`
if ( $input_var != $INTV_DAY ) then
  echo "ERROR assimilation_period_days in ${NML_DART} does not match script ${INTV_DAY}.  Exiting"
  exit
endif

set input_var = `grep assimilation_period_seconds ${NML_DART} | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g" | sed -e 's/"//g'`
if ( $input_var != $INTV_SEC ) then
  echo "ERROR assimilation_period_days in ${NML_DART} does not match script ${INTV_SEC}.  Exiting"
  exit
endif

set fn_grid_def = `grep grid_definition_filename ${NML_DART} | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g" | sed -e 's/"//g'`
if ( ! -r $fn_grid_def ) then
  echo "ERROR!  $fn_grid_def does not exist in ${RUN_DIR}, but is used for grid_definition_filename.  Exiting"
  exit
endif

set fn_model_anal = `grep model_analysis_filename ${NML_DART} | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g" | sed -e 's/"//g'`
if ( ! -r $fn_model_anal ) then
  echo "ERROR!  $fn_model_anal does not exist in ${RUN_DIR}, but is used for model_analysis_filename.  Exiting"
  exit
endif

#--------------------------------------------------------------------------
# Check for MPAS-related files and namelist entries
#--------------------------------------------------------------------------

#  Read dt and grid spacing from file.  Echo them so user is aware.
set dt = `grep config_dt $NML_MPAS | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g" | sed -e 's/"//g'`
set dx = `grep config_len_disp $NML_MPAS | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g" | sed -e 's/"//g'`

echo "MPAS is running with $dx m horizontal spacing and a timestep of $dt s"
echo

@ ndecomp = $MODEL_NODES * $N_PROCS
set fgraph = ${MPAS_GRID}.graph.info.part.${ndecomp}
if ( ! -e ${fgraph} ) then
   ${LINK} ${GRID_DIR}/${fgraph} ${fgraph}
   if(! -e ${fgraph}) then
      echo "Cannot find ${fgraph} for MODEL_NODES * N_PROCS (= $MODEL_NODES * $N_PROCS)"
      exit
   endif
endif

# Sanity checks for input files
set file_decomp = `grep config_block_decomp_file_prefix $NML_MPAS | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g" | sed -e 's/"//g'`
if ( $file_decomp != ${MPAS_GRID}.graph.info.part. ) then
  echo "config_block_decomp_file_prefix in $NML_MPAS does not match grid information provided.  Exiting"
  exit
endif

set file_sst_update = `grep config_sst_update $NML_MPAS | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g" | sed -e 's/"//g'`
if ( $SST_UPDATE == true && $file_sst_update != "true" ) then
  echo "Configuration is to update SST, but config_sst_update is not true in ${NML_MPAS}. Exiting"
  exit
endif

if ( $SST_UPDATE == true ) then
  ${LINK} ${SST_DIR}/${SST_FNAME} .		|| exit 1
else
  echo NO SST_UPDATE...
endif

ls -l ${NML_MPAS} ${STREAM_ATM} 

if ( ! -d logs ) mkdir logs			# to print out log files

#------------------------------------------
# Time info
#------------------------------------------
set greg_ini = `echo $DATE_INI 0 -g | ./advance_time`
set greg_beg = `echo $DATE_BEG 0 -g | ./advance_time`
set greg_end = `echo $DATE_END 0 -g | ./advance_time`
set  intv_hr = `expr $INTV_SEC \/ 3600`
set  intv_dh = `expr $INTV_DAY \* 24`
   @ intv_hr += $intv_dh
set diff_day = `expr $greg_end[1] \- $greg_beg[1]`
set diff_sec = `expr $greg_end[2] \- $greg_beg[2]`
set diff_tot = `expr $diff_day \* 86400 \+ $diff_sec`
set n_cycles = `expr $diff_tot \/ $INTV_SEC \+ 1`
echo "Total of ${n_cycles} cycles from $DATE_BEG to $DATE_END will be run every $intv_hr hr."
if($n_cycles < 0) then
   echo Cannot figure out how many cycles to run. Check the time setup.
   exit
endif

echo Running at $RUN_DIR
echo " "

set icyc = 1
if( $DATE_BEG != $DATE_INI ) then
   set init_day = `expr $greg_beg[1] \- $greg_ini[1]`
   set init_sec = `expr $greg_beg[2] \- $greg_ini[2]`
   set init_dif = `expr $init_day \* 86400 \+ $init_sec`
   set icyc = `expr $init_dif \/ $INTV_SEC \+ 1`
endif
set ncyc = `expr $icyc \+ $n_cycles \- 1`
echo
   
#------------------------------------------
# Initial ensemble for $DATE_INI
#------------------------------------------

set fini = `sed -n '/<immutable_stream name=\"output\"/,/\/>/{/Scree/{p;n};/##/{q};p}' ${STREAM_INIT} | \
            grep filename_template | awk -F= '{print $2}' | sed -e 's/"//g'`
set frst = `sed -n '/<immutable_stream name=\"restart\"/,/\/>/{/Scree/{p;n};/##/{q};p}' ${STREAM_ATM} | \
            grep filename_template | awk -F= '{print $2}' | awk -F$ '{print $1}' | sed -e 's/"//g'`

# ${DART_DIR}/../shell_scripts/driver_init_ens.csh 

#--------------------------------------------------------
# Cycling gets started
#--------------------------------------------------------
set time_ini = `echo $DATE_INI 0 | ./advance_time`		#YYYYMMDDHH
set time_anl = `echo $DATE_BEG 0 | ./advance_time`		#YYYYMMDDHH
set time_end = `echo $DATE_END 0 | ./advance_time`		#YYYYMMDDHH

while ( $icyc <= $ncyc )

  set time_pre = `echo $time_anl -$intv_hr | ./advance_time`	#YYYYMMDDHH
  set time_nxt = `echo $time_anl +$intv_hr | ./advance_time`	#YYYYMMDDHH
  set anal_utc = `echo $time_anl 0 -w | ./advance_time`
  set greg_obs = `echo $time_anl 0 -g | ./advance_time`
  set greg_obs_days = $greg_obs[1]
  set greg_obs_secs = $greg_obs[2]
  echo Cycle $icyc at ${time_anl}\: ${greg_obs_days}_${greg_obs_secs}

  set sav_dir = ${OUTPUT_DIR}/${time_anl}
  mkdir -p ${sav_dir}

  #------------------------------------------------------
  # 1. Namelist setup
  #------------------------------------------------------
  if($icyc == 1) then
     set cycling    = .false.
     set do_restart = .false.
  else
     set cycling    = .true.
     set do_restart = .true.
  endif
  ${REMOVE} init.sed script*.sed

  cat >! init.sed << EOF3
  /config_do_DAcycling /c\
   config_do_DAcycling = ${cycling}
  /config_do_restart /c\
   config_do_restart = ${do_restart}
EOF3
  mv ${NML_MPAS} namelist.temp
  sed -f init.sed namelist.temp >! ${NML_MPAS}
  ${REMOVE} init.sed namelist.temp

  if ( $ADAPTIVE_INF == true ) then       # For a spatially-varying prior inflation.

    if ($icyc == 1) then
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
    cat script2.sed >> script.sed

  endif

  ${MOVE} ${NML_DART} ${NML_DART}.temp
  sed -f script.sed ${NML_DART}.temp >! ${NML_DART} 		|| exit 2
  ${REMOVE} script.sed script2.sed ${NML_DART}.temp

  #------------------------------------------------------
  # 2. Update input files to get filter started 
  # (assuming start_from_restart = .true. in input.nml)
  #------------------------------------------------------
  set f_rst = ${frst}`echo ${anal_utc} | sed -e 's/:/\./g'`.nc
  set f_anl = analysis.`echo ${f_rst} | cut -d . -f2`.nc

  echo "Input ensemble for ${time_anl}"
  if( -e ${input_list})  ${REMOVE} ${input_list}
  if( -e ${output_list}) ${REMOVE} ${output_list}
  echo ${input_list} ${output_list}
  echo
  set i = 1
  while ( $i <= ${ENS_SIZE} )
    if (! -e ${ENS_DIR}${i}/${f_rst}) then
      echo "Cannot find ${ENS_DIR}${i}/${f_rst}".
    else
      echo ${ENS_DIR}${i}/${f_rst} >> ${input_list}
    endif
    echo ${ENS_DIR}${i}/${f_anl} >> ${output_list}
    @ i++
  end
  tail -1 ${input_list}
  tail -1 ${output_list}

  set ne = `cat ${input_list} | wc -l `
  if ( $ne != $ENS_SIZE ) then
     echo "We need ${ENS_SIZE} initial ensemble members, but found ${ne} only."
     exit
  endif 
  echo

  if ( $ADAPTIVE_INF == true && $icyc > 1 ) then
    if ( ! -e ${OUTPUT_DIR}/${time_pre}/${INFL_OUT}_mean.nc ) then
      echo ${OUTPUT_DIR}/${time_pre}/${INFL_OUT}_mean.nc does not exist. Stop.
      exit
    endif
    ${LINK} ${OUTPUT_DIR}/${time_pre}/${INFL_OUT}_mean.nc ${INFL_IN}_mean.nc
    ${LINK} ${OUTPUT_DIR}/${time_pre}/${INFL_OUT}_sd.nc   ${INFL_IN}_sd.nc
  endif

  #------------------------------------------------------
  # 3. Obs sequence for this analysis cycle - one obs time at each analysis cycle
  #------------------------------------------------------
  set fn_obs = ${OBS_DIR}/${obs_seq_in}.${time_anl}
  if ( ! -e ${fn_obs} ) then
     echo ${fn_obs} does not exist. Stop.
     exit
  endif
  ${LINK} ${fn_obs} ${obs_seq_in}

  #------------------------------------------------------
  # 4. Run filter
  #------------------------------------------------------
  set job_name = ${EXPERIMENT_NAME}.${icyc}
  echo Running filter: $job_name

  if ( $RUN_IN_PBS == yes ) then

    cat >! filter.sed << EOF
    s#JOB_NAME#${job_name}#g
    s#PROJ_NUMBER#${PROJ_NUMBER}#g
    s#NODES#${FILTER_NODES}#g
    s#NCPUS#${N_CPUS}#g
    s#JOB_TIME#${TIME_FILTER}#g
    #s#QUEUE#${QUEUE}#g
EOF

    sed -f filter.sed filter.template.pbs >! filter.pbs
    qsub filter.pbs
    ${REMOVE} filter.sed

    # Wait until the job is finished.
    set is_there = `qstat | grep $job_name | wc -l`
    while ( $is_there != 0 )
      sleep 60
      set is_there = `qstat | grep $job_name | wc -l`
    end
    ${MOVE} ${job_name}.o* logs/.

  else

    echo `date +%s` >&! filter_started
    ./filter >! filter.log
    if ( -e ${obs_seq_out})  touch filter_done

  endif

  # Check errors in filter.
  if ( -e filter_started && ! -e filter_done ) then
    echo "Filter was not normally finished. Exiting."
    ${REMOVE} filter_started
    exit
  endif

  ${REMOVE} filter_started filter_done
  echo Filter is done for Cycle ${icyc}\: ${time_anl}

  #------------------------------------------------------
  # 5. Target time for model advance
  #------------------------------------------------------
  set greg_obs = `echo $time_anl ${INTV_DAY}d${INTV_SEC}s -g | ./advance_time`
  set greg_obs_days = $greg_obs[1]
  set greg_obs_secs = $greg_obs[2]
  echo Target date: $time_nxt ${greg_obs_days}_${greg_obs_secs}
  ${COPY} input.nml ${OUTPUT_DIR}/${time_anl}/input.nml.filter.${icyc}

  #------------------------------------------------------
  # 6. Run update_mpas_states for all ensemble members
  #------------------------------------------------------
  ${DART_DIR}/update_mpas_states >! logs/update_mpas_states.ic_${icyc}.log

  #------------------------------------------------------
  # 7. Advance model for each member
  #------------------------------------------------------
  # Run forecast for ensemble members until next analysis time
  echo Advance models for ${ENS_SIZE} members now...

  set n = 1
  while ( $n <= $ENS_SIZE )

    if ( $RUN_IN_PBS == yes ) then

      set job_ensemble = ${EXPERIMENT_NAME}_${icyc}_ens${n}

      cat >! advance.sed << EOF
      s#JOB_NAME#${job_ensemble}#g
      s#PROJ_NUMBER#${PROJ_NUMBER}#g
      s#ENS_MEM#${n}#g
      s#QUEUE#${QUEUE}#g
      s#NCPUS#${N_CPUS}#g
      s#NODES#${MODEL_NODES}#g
      s#NPROC#${N_PROCS}#g
      s#JOB_TIME#${TIME_MPAS}#g
EOF

      sed -f advance.sed advance_model.template >! advance_model.pbs
      qsub advance_model.pbs
      ${REMOVE} advance.sed
      sleep 1

    else

      ./advance_model.csh $n $n >! logs/advance_model.${icyc}.${n}.log

    endif

    @ n++
  
  end

  if ( $RUN_IN_PBS == yes ) then

    # Check if all members are done advancing model.
    set is_all_done = `qstat | grep $job_ensemble | wc -l`
    while ( $is_all_done > 0 )
      sleep 60
      set is_all_done = `qstat | grep $job_ensemble | wc -l`
    end
    sleep 60

  endif

  #------------------------------------------------------
  # 8. Store output files
  #------------------------------------------------------
  echo Saving output files for ${time_anl}.
  ls -lrt >> ${sav_dir}/list
  set fstat = `grep stages_to_write input.nml | awk -F= '{print $2}'`
  set fs = `echo $fstat | sed -e 's/,/ /g' | sed -e "s/'//g"`

  set hdir = ${HPSS_DIR}/${time_anl}
  if ( ${HPSS_SAVE} == yes ) then
    echo $hdir
    hsi mkdir -p $hdir
    if ( ! $status == 0 ) then
      echo We cannot access to ${hdir}. We back up files locally.
      set HPSS_SAVE = no
    endif
  endif
  echo "HSI BACKUP? ${HPSS_SAVE} in ${hdir}/"

  foreach f ( $fs )
    ${MOVE} ${f}*.nc ${sav_dir}/
  end
  ${MOVE} ${obs_seq_out} ${sav_dir}/.
 
  if ( $HPSS_SAVE == yes ) then
    cd ${sav_dir}
    foreach f ( * )
      gzip -f $f
      ${HSICMD} ${f}.gz : ${hdir}/${f}.gz
    end
  endif

  #------------------------------------------------------
  # 9. Get ready to run filter for next cycle.
  #------------------------------------------------------
  cd $RUN_DIR
  mv ${input_list} ${sav_dir}/.

  set fcst_utc = `echo $time_nxt 0 -w | ./advance_time`
  set f_fcst = `head -1 list.${time_nxt}.txt`
  set fout = `basename ${f_fcst}`

  set n = 1
  while ( $n <= ${ENS_SIZE} )
    if( ! -e ${ENS_DIR}${n}/${fout} ) then
      echo Missing ${ENS_DIR}${n}/${fout}
      echo ${ENS_DIR}${n}/${fout} >> missing.${time_nxt}.txt
    else 
      echo ${ENS_DIR}${n}/${fout} >> ${input_list}
    endif
    if ( -e logs/start_member${n} )  ${REMOVE} logs/start_member${n}
    @ n++
  end

  echo Filter is ready to go for the next cycle now.
  set time_anl = $time_nxt
  @ icyc++

  exit

end

echo Cycling is done for $n_cycles cycles in ${EXPERIMENT_NAME}.
echo Script exiting normally.

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

