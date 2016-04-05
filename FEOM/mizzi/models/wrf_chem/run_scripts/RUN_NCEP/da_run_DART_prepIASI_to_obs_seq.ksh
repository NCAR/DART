#!/bin/ksh -aeux 
#
  set echo
#
# MODIFIED VERSION OF /DART/models/WRF/regression/CONUS-V3/icbc_real.ksh
# TO SETUP AN ENVIRONMENT TO CONVERT OBSERVATIONS TO obs_seq AND RUN dart-wrf-filter FOR A SINGLE STEP.
#
# SET SWITCHES FOR SCRIT OPTIONS
  export RUN_PREP_TO_ASCII=true
  export RUN_PREP_ASCII_TO_DART=true
  export RUN_IASI_ASCII_TO_DART=true
  export RUN_MDART_TO_SDART=true
#
# SET TIME INFORMATION
  export START_DATE=2008061106
  export END_DATE=2008062100
  export TIME_INC=6
  export ASIM_WINDOW=3
#
# SYSTEM SPECIFIC SETTINGS
  export PROCS=8
#  export OB_TYPE=ob_reanal
  export OB_TYPE=obs
#
# INITIAL CONDITION FILES
# Set to '1' if you want a single IC file, any other # if you want separate files (the 
# latter is suggested if you have large grids and lots of members)
  export single_file=1
#
# PATHS
  export WRFDA_VER=WRFDAv3.4_dmpar
  export WRF_VER=WRFv3.4_dmpar
  export DART_VER=DART_DEVEL
#
# INDEPENDENT DIRECTORIES
  export CODE_DIR=/glade/home/mizzi/TRUNK
  export DATA_DIR=/ptmp/mizzi/AVE_TEST_DATA
  export ASIM_DIR=/ptmp/mizzi/dart_assim
#
# DEPENDENT DIRECTORIES
  export HYBRID_DIR=${CODE_DIR}/HYBRID_TRUNK
  export WRF_DIR=${CODE_DIR}/${WRF_VER}
  export VAR_DIR=${CODE_DIR}/${WRFDA_VER}
  export BUILD_DIR=${VAR_DIR}/var/build
  export DART_DIR=${CODE_DIR}/${DART_VER}
  export TOOL_DIR=${VAR_DIR}/var/da
  export ICBC_DIR=${ASIM_DIR}
  export HYBRID_SCRIPTS_DIR=${HYBRID_DIR}/hybrid_scripts
#
# MAKE ASSIMILATION DIRECTORY AND GO TO IT
  if [[ ! -d ${ASIM_DIR} ]]; then mkdir -p ${ASIM_DIR}; fi
  mkdir -p ${ASIM_DIR}
#
# BEGIN DAY AND TIME LOOP
  export L_DATE=${START_DATE}
  while [[ ${L_DATE} -le ${END_DATE} ]]; do
     export YYYY=$(echo $L_DATE | cut -c1-4)
     export YY=$(echo $L_DATE | cut -c3-4)
     export MM=$(echo $L_DATE | cut -c5-6)
     export DD=$(echo $L_DATE | cut -c7-8)
     export HH=$(echo $L_DATE | cut -c9-10)
     export PAST_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} -${ASIM_WINDOW} 2>/dev/null)  
     export PAST_YYYY=$(echo $PAST_DATE | cut -c1-4)
     export PAST_MM=$(echo $PAST_DATE | cut -c5-6)
     export PAST_DD=$(echo $PAST_DATE | cut -c7-8)
     export PAST_HH=$(echo $PAST_DATE | cut -c9-10)
#
     export IASI_PAST_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} -24 2>/dev/null)  
     export IASI_PAST_YYYY=$(echo $IASI_PAST_DATE | cut -c1-4)
     export IASI_PAST_MM=$(echo $IASI_PAST_DATE | cut -c5-6)
     export IASI_PAST_DD=$(echo $IASI_PAST_DATE | cut -c7-8)
     export IASI_PAST_HH=$(echo $IASI_PAST_DATE | cut -c9-10)
     export IASI_NEXT_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} +24 2>/dev/null)  
     export IASI_NEXT_YYYY=$(echo $IASI_NEXT_DATE | cut -c1-4)
     export IASI_NEXT_MM=$(echo $IASI_NEXT_DATE | cut -c5-6)
     export IASI_NEXT_DD=$(echo $IASI_NEXT_DATE | cut -c7-8)
     export IASI_NEXT_HH=$(echo $IASI_NEXT_DATE | cut -c9-10)
#
# DART TIME INFO (NO LEADING ZEROS)
     export DT_YYYY=${YYYY}
     export DT_YY=$(echo $L_DATE | cut -c3-4)
     export DT_MM=${MM} 
     export DT_DD=${DD} 
     export DT_HH=${HH} 
     (( DT_MM = ${DT_MM} + 0 ))
     (( DT_DD = ${DT_DD} + 0 ))
     (( DT_HH = ${DT_HH} + 0 ))
#    
     export YEAR_INIT=${DT_YYYY}
     export MONTH_INIT=${DT_MM}
     export DAY_INIT=${DT_DD}
     export HOUR_INIT=${DT_HH}
     export YEAR_END=${DT_YYYY}
     export MONTH_END=${DT_MM}
     export DAY_END=${DT_DD}
     export HOUR_END=${DT_HH}
     export DA_TIME_WINDOW=0
#
# RUN PREP_TO_ASCII
     if ${RUN_PREP_TO_ASCII}; then
#
# APM: USE MODIFIED VERSION OF input.nml and prepbufr.csh from 
# DART/observations/NCEP/prep_bufr/work
#
# APM: USE MODIFIED VERSION create_real_obs.f90 from
# DART/observations/NCEP/ascii_to_obs
        cd ${ASIM_DIR}
        export OBS_PREP_DIR=${DATA_DIR}/${OB_TYPE}/${L_DATE}
        if [[ ! -d ${ASIM_DIR}/prep_bufr/${YYYY}${MM} ]]; then mkdir -p ${ASIM_DIR}/prep_bufr/${YYYY}${MM}; fi
        cd ${ASIM_DIR}/prep_bufr/${YYYY}${MM}
#
        if [[ ${OB_TYPE} = "prepqm" ]]; then
           cp ${OBS_PREP_DIR}/prepqm${YY}${MM}${DD}${HH} .
        else
           cp ${OBS_PREP_DIR}/ob.bufr prepqm${YY}${MM}${DD}${HH}
        fi
        cp ${DART_DIR}/observations/NCEP/prep_bufr/work/input.nml .
        ${DART_DIR}/observations/NCEP/prep_bufr/work/prepbufr.csh
        cd ${ASIM_DIR}
     fi
#
# RUN PREP_ASCII_TO_DART
     if ${RUN_PREP_ASCII_TO_DART}; then 
        cd ${ASIM_DIR}
        if [[ ! -d ${ASIM_DIR}/prep_ascii_to_dart/${YYYY}${MM} ]]; then mkdir -p ${ASIM_DIR}/prep_ascii_to_dart/${YYYY}${MM}; fi
        cd ${ASIM_DIR}/prep_ascii_to_dart/${YYYY}${MM}
        if [[ ${HH} -eq 0 ]] then
           export L_YYYY=${PAST_YYYY}
           export L_MM=${PAST_MM}
           export L_DD=${PAST_DD}
           export L_HH=24
           export D_DATE=${L_YYYY}${L_MM}${L_DD}${L_HH}
        else
           export L_YYYY=${YYYY}
           export L_MM=${MM}
           export L_DD=${DD}
           export L_HH=${HH}
           export D_DATE=${L_YYYY}${L_MM}${L_DD}${L_HH}
        fi
        cp ${ASIM_DIR}/prep_bufr/${YYYY}${MM}/temp_obs.${D_DATE} .
        ${HYBRID_SCRIPTS_DIR}/da_create_dart_ncep_ascii_to_obs_input_nml.ksh
        ${DART_DIR}/observations/NCEP/ascii_to_obs/work/create_real_obs
        cd ${ASIM_DIR}
     fi
#
# RUN_IASI_ASCII_TO_DART
     if ${RUN_IASI_ASCII_TO_DART}; then 
        cd ${ASIM_DIR}
        if [[ ! -d ${ASIM_DIR}/IASI_ascii_to_dart/${YYYY}${MM} ]]; then mkdir -p ${ASIM_DIR}/IASI_ascii_to_dart/${YYYY}${MM}; fi
        cd ${ASIM_DIR}/IASI_ascii_to_dart/${YYYY}${MM}
        if [[ ${HH} -eq 0 ]] then
           export L_YYYY=${PAST_YYYY}
           export L_MM=${PAST_MM}
           export L_DD=${PAST_DD}
           export L_HH=24
           export D_DATE=${L_YYYY}${L_MM}${L_DD}${L_HH}
        else
           export L_YYYY=${YYYY}
           export L_MM=${MM}
           export L_DD=${DD}
           export L_HH=${HH}
           export D_DATE=${L_YYYY}${L_MM}${L_DD}${L_HH}
        fi
#
# USE IASI DATA FROM PAST, PRESENT, AND NEXT DATES TO ENSURE FULL COVERAGE OF $ASIM_WINDOW 
        if [[ -e input.nml ]]; then rm -rf input.nml; fi
        export OBS_IASI_DIR=${DATA_DIR}/obs_IASI/${IASI_PAST_YYYY}${IASI_PAST_MM}${IASI_PAST_DD}00
        cp ${OBS_IASI_DIR}/IASIO3PROF_OBSSEQ_method2_${IASI_PAST_YYYY}${IASI_PAST_MM}${IASI_PAST_DD}.dat iasi_asciidata.input
        cp ${DART_DIR}/observations/IASI/work/input.nml .
        ${DART_DIR}/observations/IASI/work/iasi_ascii_to_obs
        export IASI_PAST_FILE=iasi_obs_seq_${IASI_PAST_YYYY}${IASI_PAST_MM}${IASI_PAST_DD}
        cp iasi_obs_seq.out ${IASI_PAST_FILE}
#
        if [[ -e input.nml ]]; then rm -rf input.nml; fi
        export OBS_IASI_DIR=${DATA_DIR}/obs_IASI/${YYYY}${MM}${DD}00
        cp ${OBS_IASI_DIR}/IASIO3PROF_OBSSEQ_method2_${YYYY}${MM}${DD}.dat iasi_asciidata.input
        cp ${DART_DIR}/observations/IASI/work/input.nml .
        ${DART_DIR}/observations/IASI/work/iasi_ascii_to_obs
        export IASI_PRES_FILE=iasi_obs_seq_${YYYY}${MM}${DD}
        cp iasi_obs_seq.out ${IASI_PRES_FILE}
#
        if [[ -e input.nml ]]; then rm -rf input.nml; fi
        export OBS_IASI_DIR=${DATA_DIR}/obs_IASI/${IASI_NEXT_YYYY}${IASI_NEXT_MM}${IASI_NEXT_DD}00
        cp ${OBS_IASI_DIR}/IASIO3PROF_OBSSEQ_method2_${IASI_NEXT_YYYY}${IASI_NEXT_MM}${IASI_NEXT_DD}.dat iasi_asciidata.input
        cp ${DART_DIR}/observations/IASI/work/input.nml .
        ${DART_DIR}/observations/IASI/work/iasi_ascii_to_obs
        export IASI_NEXT_FILE=iasi_obs_seq_${IASI_NEXT_YYYY}${IASI_NEXT_MM}${IASI_NEXT_DD}
        cp iasi_obs_seq.out ${IASI_NEXT_FILE}
#
        cd ${ASIM_DIR}
     fi
#
# RUN_MDART_TO_SDART (CONVERT MANY OBS_SEQ FILES TO A SINGLE OBS_SEQ FILE AND FILTER)
     if ${RUN_MDART_TO_SDART}; then
        cd ${ASIM_DIR}
        if [[ ! -d ${ASIM_DIR}/mdart_to_sdart/${YYYY}${MM} ]]; then mkdir -p ${ASIM_DIR}/mdart_to_sdart/${YYYY}${MM}; fi
        cd ${ASIM_DIR}/mdart_to_sdart/${YYYY}${MM}
        if [[ ${HH} -eq 0 ]] then
           export L_YYYY=${PAST_YYYY}
           export L_MM=${PAST_MM}
           export L_DD=${PAST_DD}
           export L_HH=24
           export D_DATE=${L_YYYY}${L_MM}${L_DD}${L_HH}
        else
           export L_YYYY=${YYYY}
           export L_MM=${MM}
           export L_DD=${DD}
           export L_HH=${HH}
           export D_DATE=${L_YYYY}${L_MM}${L_DD}${L_HH}
        fi
        cp ${DART_DIR}/models/wrf/work/advance_time ./.
        cp ${DART_DIR}/models/wrf/work/obs_sequence_tool ./.
        cp ${DART_DIR}/models/wrf/work/input.nml ./.
        cp ${ASIM_DIR}/prep_ascii_to_dart/${YYYY}${MM}/obs_seq${D_DATE} ./.
        cp ${ASIM_DIR}/IASI_ascii_to_dart/${YYYY}${MM}/${IASI_PAST_FILE} ./.
        cp ${ASIM_DIR}/IASI_ascii_to_dart/${YYYY}${MM}/${IASI_PRES_FILE} ./.
        cp ${ASIM_DIR}/IASI_ascii_to_dart/${YYYY}${MM}/${IASI_NEXT_FILE} ./.
#
# CALCULATE GREGORIAN TIMES FOR START AND END OF ASSIMILAtION WINDOW
        export ASIM_MIN_DATE=$($BUILD_DIR/da_advance_time.exe $L_DATE -$ASIM_WINDOW 2>/dev/null)
        export ASIM_MAX_DATE=$($BUILD_DIR/da_advance_time.exe $L_DATE +$ASIM_WINDOW 2>/dev/null)
        set -A temp `echo $ASIM_MIN_DATE 0 -g | ./advance_time`
        export ASIM_MIN_DAY_GREG=${temp[0]}
        export ASIM_MIN_SEC_GREG=${temp[1]}
        set -A temp `echo $ASIM_MAX_DATE 0 -g | ./advance_time` 
        export ASIM_MAX_DAY_GREG=${temp[0]}
        export ASIM_MAX_SEC_GREG=${temp[1]}
#
# SETUP OBS_SEQUENCE_TOOL INPUT.NML
        export NL_NUM_INPUT_FILES=4
        export NL_FILENAME_SEQ="'obs_seq${D_DATE}','${IASI_PAST_FILE}','${IASI_PRES_FILE}','${IASI_NEXT_FILE}'"
        export NL_FILENAME_OUT="'obs_seq_${L_DATE}.out'"
        export NL_FIRST_OBS_DAYS=${ASIM_MIN_DAY_GREG}
        export NL_FIRST_OBS_SECONDS=${ASIM_MIN_SEC_GREG}
        export NL_LAST_OBS_DAYS=${ASIM_MAX_DAY_GREG}
        export NL_LAST_OBS_SECONDS=${ASIM_MAX_SEC_GREG}
        export NL_SYNONYMOUS_COPY_LIST="'NCEP BUFR observation','observation'"
        export NL_SYNONYMOUS_QC_LIST="'NCEP QC index','Data QC'"
        ${HYBRID_SCRIPTS_DIR}/da_create_dart_input_nml.ksh       
        ./obs_sequence_tool
        if [[ ! -d ${DATA_DIR}/obs_IASI/${L_DATE} ]]; then mkdir -p ${DATA_DIR}/obs_IASI/${L_DATE}; fi
        cp obs_seq_${L_DATE}.out ${DATA_DIR}/obs_IASI/${L_DATE}/.
        cd ${ASIM_DIR}
     fi 
#
# LOOP TO NEXT DAY AND TIME 
     export P_DATE=${L_DATE}
     export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${TIME_INC} 2>/dev/null)  
  done 
exit
