#!/bin/csh
########################################################################
#
#   init_mpas_grib.csh - shell script that can be used to create an 
#                        initial MPAS ensemble from an external grib
#                        file, then allows the user to integrate the
#                        forecast forward in time, so the forecasts are
#                        available for cycling.
#
########################################################################

  set ensemble_member = ${1}     #  ensemble member
  set paramfile       = ${2}	 #  the parameter file with variables

  source $paramfile

  set statfile = `sed -n '/\"input\"/,/\/>/{/Scree/{p;n};/##/{q};p}' ${RUN_DIR}/${STREAM_INIT} | grep filename_template | awk -F\" '{print $(NF-1)}'`
  
  foreach fn ( ${RUN_DIR}/MPAS_RUN/${statfile} )
     if ( ! -r ${fn} ) then
        echo ABORT\: init_mpas_grib.csh could not find required readable dependency ${fn}
        exit 1
     endif
  end

  set temp_dir = ${INIT_DIR}/${ENS_DIR}${ensemble_member}

  if ( -d ${temp_dir} )  ${REMOVE} ${temp_dir}
  mkdir -p ${temp_dir}
  ${LINK} ${temp_dir} ${RUN_DIR}/${temp_dir}
  cd ${RUN_DIR}/${temp_dir}

  # Get the program and necessary files for the model
  ${COPY} ${RUN_DIR}/input.nml                         .   || exit 1                    
  ${LINK} ${RUN_DIR}/MPAS_RUN/atmosphere_model         .   || exit 1
  ${LINK} ${RUN_DIR}/MPAS_RUN/init_atmosphere_model    .   || exit 1
  ${LINK} ${RUN_DIR}/MPAS_RUN/ungrib.exe               .   || exit 1
  ${LINK} ${RUN_DIR}/MPAS_RUN/stream*                  .   || exit 1
  ${LINK} ${RUN_DIR}/MPAS_RUN/*BL                      .   || exit 1
  ${LINK} ${RUN_DIR}/MPAS_RUN/*DATA                    .   || exit 1
  ${LINK} ${RUN_DIR}/advance_time                      .   || exit 1
  ${LINK} ${RUN_DIR}/*graph*                           .   || exit 1             
  ${LINK} ${RUN_DIR}/Vtable                            .   || exit 1

  #  Determine the initial, final and run times for the MPAS integration
  set curr_utc = `echo $DATE_INI -${INIT_FORECAST_LENGTH}h -w | ./advance_time`
  set targ_utc = `echo $DATE_INI 0 -w | ./advance_time`

  set fdays = 0
  set fhours = $INIT_FORECAST_LENGTH
  while ( $fhours >= 24 )
     @ fdays++
     @ fhours = $fhours - 24
  end
  set intv_utc = `echo $fdays + 100 | bc | cut -b2-3`_`echo $fhours + 100 | bc | cut -b2-3`:00:00

  #  Link grib file, run WPS ungrib program to convert into format that MPAS can use
  ${LINK} `head -n $ensemble_member ${INIT_GRIB_FILE_LIST} | tail -1` GRIBFILE.AAA

  cat >! script.sed << EOF
  /start_date/c\
  start_date = '${curr_utc}',
  /end_date/c\
  end_date   = '${curr_utc}',
EOF
  sed -f script.sed ${RUN_DIR}/namelist.wps >! namelist.wps

  ./ungrib.exe >& ungrib.out

  #  create namelist file for init version of MPAS
  cat >! script.sed << EOF
  /config_start_time/c\
  config_start_time = '$curr_utc'
  /config_stop_time/c\
  config_stop_time = '$curr_utc'
EOF
  sed -f script.sed ${RUN_DIR}/MPAS_RUN/${NML_INIT} >! ${NML_INIT}

  # Get the grid info files - now for PIO
  set is_grid_info = `grep config_block_decomp_file_prefix ${RUN_DIR}/${NML_MPAS} | wc -l`
  if( $is_grid_info != 1 ) then
     echo Cannot find grid info. Stop.
     exit
  endif
  set fs_grid = `grep config_block_decomp_file_prefix ${RUN_DIR}/${NML_MPAS} | awk '{print $3}' | sed -e "s/'//g"`
  ${LINK} ${RUN_DIR}/MPAS_RUN/${fs_grid}* .
  
  ${LINK} ${RUN_DIR}/MPAS_RUN/${statfile} .

  #  Run init version of MPAS to create initial condition file
#  mpirun.lsf ./init_atmosphere_model  || exit 2
  mpiexec_mpt dplace -s 1 ./init_atmosphere_model  || exit 2
  ${REMOVE} FILE:*

  #  Generate MPAS namelist file
  cat >! script.sed << EOF
  /config_start_time/c\
  config_start_time = '$curr_utc'
  /config_run_duration/c\
  config_run_duration = '$intv_utc'
  /config_do_restart/c\
  config_do_restart = false
  /config_do_DAcycling/c\
  config_do_DAcycling = false
EOF
  sed -f script.sed ${RUN_DIR}/${NML_MPAS} >! ${NML_MPAS}

  # clean out any old rsl files
  if ( -e log.0000.out ) ${REMOVE} log.*

  #  Run MPAS for the specified amount of time 
#  mpirun.lsf ./atmosphere_model        || exit 3
  mpiexec_mpt dplace -s 1 ./atmosphere_model   || exit 3

  # Check the output file
  set frst = `sed -n '/<immutable_stream name=\"restart\"/,/\/>/{/Scree/{p;n};/##/{q};p}' ${STREAM_ATM} | \
              grep filename_template | awk -F= '{print $2}' | awk -F$ '{print $1}' | sed -e 's/"//g'`

  set fout = ${frst}`echo ${targ_utc} | sed -e 's/:/\./g'`.nc
 
  foreach rfile ( `ls -1 ${frst}*` )
    if ( $rfile != $fout )  ${REMOVE} $rfile
  end
 
  cd $RUN_DIR

