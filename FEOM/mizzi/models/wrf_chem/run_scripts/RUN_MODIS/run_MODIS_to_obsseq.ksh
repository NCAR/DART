#!/bin/ksh -x 
#
  set echo
#
# MODIFIED VERSION OF /DART/models/WRF/regression/CONUS-V3/icbc_real.ksh
# TO SETUP AN ENVIRONMENT TO CONVERT OBSERVATIONS TO obs_seq.
#
# SET SWITCHES FOR SCRIT OPTIONS
  export RUN_PREP_TO_ASCII=false
  export RUN_PREP_ASCII_TO_DART=false
  export RUN_MODIS_ASCII_TO_DART=true
  export RUN_MDART_TO_SDART=true
#
# SET TIME INFORMATION
  export START_DATE=2008060106
  export END_DATE=2008063018
  export TIME_INC=6
  export ASIM_WINDOW=3
  export START_DATE_DATA=2008060100
  export END_DATE_DATA=2008063000
#
# SYSTEM SPECIFIC SETTINGS
  export PROCS=8
#  export OB_TYPE=ob_reanal
  export OB_TYPE=obs
  export NL_APM_SCALE=1.
  export NL_APM_SCALE_SW=.FALSE.
#
# INITIAL CONDITION FILES
# Set to '1' if you want a single IC file, any other # if you want separate files (the 
# latter is suggested if you have large grids and lots of members)
  export single_file=1
#
# PATHS
  export WRFDA_VER=WRFDAv3.4_dmpar
  export WRF_VER=WRFv3.4_dmpar
  export DART_VER=DART_CHEM
#
# INDEPENDENT DIRECTORIES
  export CODE_DIR=/glade/p/work/mizzi
  export ACD_DIR=/glade/p/acd/mizzi
  export SCRATCH_DIR=/glade/scratch/mizzi
#
# DEPENDENT DIRECTORIES
  export ASIM_DIR=${SCRATCH_DIR}/MODIS_to_OBSSEQ
  export TRUNK_DIR=${CODE_DIR}/TRUNK
  export HYBRID_DIR=${CODE_DIR}/HYBRID_TRUNK
  export WRF_DIR=${TRUNK_DIR}/${WRF_VER}
  export VAR_DIR=${TRUNK_DIR}/${WRFDA_VER}
  export BUILD_DIR=${VAR_DIR}/var/build
  export DART_DIR=${TRUNK_DIR}/${DART_VER}
  export TOOL_DIR=${VAR_DIR}/var/da
  export ICBC_DIR=${ASIM_DIR}
  export HYBRID_SCRIPTS_DIR=${HYBRID_DIR}/hybrid_scripts
  export DATA_DIR=${ACD_DIR}/AVE_TEST_DATA
  export OBS_MODIS_DIR=${DATA_DIR}/obs_MODIS_AOD_ascii
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
     export MODIS_PAST_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} -24 2>/dev/null)  
     export MODIS_PAST_YYYY=$(echo $MODIS_PAST_DATE | cut -c1-4)
     export MODIS_PAST_MM=$(echo $MODIS_PAST_DATE | cut -c5-6)
     export MODIS_PAST_DD=$(echo $MODIS_PAST_DATE | cut -c7-8)
     export MODIS_PAST_HH=$(echo $MODIS_PAST_DATE | cut -c9-10)
     export MODIS_NEXT_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} +24 2>/dev/null)  
     export MODIS_NEXT_YYYY=$(echo $MODIS_NEXT_DATE | cut -c1-4)
     export MODIS_NEXT_MM=$(echo $MODIS_NEXT_DATE | cut -c5-6)
     export MODIS_NEXT_DD=$(echo $MODIS_NEXT_DATE | cut -c7-8)
     export MODIS_NEXT_HH=$(echo $MODIS_NEXT_DATE | cut -c9-10)
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
# RUN_MODIS_ASCII_TO_DART
     if ${RUN_MODIS_ASCII_TO_DART}; then 
        cd ${ASIM_DIR}
        if [[ ! -d ${ASIM_DIR}/MODIS_ascii_to_dart/${YYYY}${MM} ]]; then mkdir -p ${ASIM_DIR}/MODIS_ascii_to_dart/${YYYY}${MM}; fi
        cd ${ASIM_DIR}/MODIS_ascii_to_dart/${YYYY}${MM}
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
# USE MODIS DATA FROM PAST, PRESENT, AND NEXT DATES TO ENSURE FULL COVERAGE OF $ASIM_WINDOW 
#
# RUN FOR PAST
        if [[ ${MODIS_PAST_DATE} -ge ${START_DATE_DATA} ]]; then
           rm -rf input.nml
           rm -rf modis_asciidata.input
           rm -rf modis_obs_seq.out
           if [[ -e ${OBS_MODIS_DIR}/modis_ascii_${MODIS_PAST_YYYY}${MODIS_PAST_MM}${MODIS_PAST_DD}00.input ]]; then
              cp ${OBS_MODIS_DIR}/modis_ascii_${MODIS_PAST_YYYY}${MODIS_PAST_MM}${MODIS_PAST_DD}00.input modis_asciidata.input
              cp ${DART_DIR}/observations/MODIS/work/input.nml ./.
              ${DART_DIR}/observations/MODIS/work/modis_ascii_to_obs
              export MODIS_PAST_FILE=modis_obs_seq_${MODIS_PAST_YYYY}${MODIS_PAST_MM}${MODIS_PAST_DD}
              mv modis_obs_seq.out ${MODIS_PAST_FILE}
           fi
        fi
#
# RUN FOR PRESENT
        rm -rf input.nml
        rm -rf modis_asciidata.input
        if [[ -e ${OBS_MODIS_DIR}/modis_ascii_${YYYY}${MM}${DD}00.input ]]; then
           cp ${OBS_MODIS_DIR}/modis_ascii_${YYYY}${MM}${DD}00.input modis_asciidata.input
           cp ${DART_DIR}/observations/MODIS/work/input.nml .
           ${DART_DIR}/observations/MODIS/work/modis_ascii_to_obs
           export MODIS_PRES_FILE=modis_obs_seq_${YYYY}${MM}${DD}
           mv modis_obs_seq.out ${MODIS_PRES_FILE}
        fi
#
# RUN FOR NEXT
        if [[ ${MODIS_NEXT_DATE} -le ${END_DATE_DATA} ]]; then
           rm -rf input.nml
           rm -rf modis_asciidata.input
           if [[ -e ${OBS_MODIS_DIR}/modis_ascii_${MODIS_NEXT_YYYY}${MODIS_NEXT_MM}${MODIS_NEXT_DD}00.input ]]; then
              cp ${OBS_MODIS_DIR}/modis_ascii_${MODIS_NEXT_YYYY}${MODIS_NEXT_MM}${MODIS_NEXT_DD}00.input modis_asciidata.input
              cp ${DART_DIR}/observations/MODIS/work/input.nml .
              ${DART_DIR}/observations/MODIS/work/modis_ascii_to_obs
              export MODIS_NEXT_FILE=modis_obs_seq_${MODIS_NEXT_YYYY}${MODIS_NEXT_MM}${MODIS_NEXT_DD}
              mv modis_obs_seq.out ${MODIS_NEXT_FILE}
           fi           
        fi
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
        cp ${DART_DIR}/models/wrf_chem/work/advance_time ./.
        cp ${DART_DIR}/models/wrf_chem/work/obs_sequence_tool ./.
        cp ${DART_DIR}/models/wrf_chem/work/input.nml ./.

        export MODIS_PAST_FILE=modis_obs_seq_${MODIS_PAST_YYYY}${MODIS_PAST_MM}${MODIS_PAST_DD}
        export MODIS_PRES_FILE=modis_obs_seq_${YYYY}${MM}${DD}
        export MODIS_NEXT_FILE=modis_obs_seq_${MODIS_NEXT_YYYY}${MODIS_NEXT_MM}${MODIS_NEXT_DD}
        if ${RUN_PREP_ASCII_TO_DART}; then 
           cp ${ASIM_DIR}/prep_ascii_to_dart/${YYYY}${MM}/obs_seq${D_DATE} ./.
        fi
        if ${RUN_MODIS_ASCII_TO_DART}; then
           if [[ ${MODIS_PAST_DATE} -ge ${START_DATE_DATA} ]]; then 
              cp ${ASIM_DIR}/MODIS_ascii_to_dart/${YYYY}${MM}/${MODIS_PAST_FILE} ./.
           fi
           cp ${ASIM_DIR}/MODIS_ascii_to_dart/${YYYY}${MM}/${MODIS_PRES_FILE} ./.
           if [[ ${MODIS_NEXT_DATE} -le ${END_DATE_DATA} ]]; then 
              cp ${ASIM_DIR}/MODIS_ascii_to_dart/${YYYY}${MM}/${MODIS_NEXT_FILE} ./.
           fi
        fi
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
        if ${RUN_PREP_ASCII_TO_DART}; then 
           export NL_NUM_INPUT_FILES=4
           export NL_FILENAME_SEQ="'obs_seq${D_DATE}','${MODIS_PAST_FILE}','${MODIS_PRES_FILE}','${MODIS_NEXT_FILE}'"
        elif [[  ${MODIS_PAST_DATE} -lt ${START_DATE_DATA} ]]; then
           export NL_NUM_INPUT_FILES=2
           export NL_FILENAME_SEQ="'${MODIS_PRES_FILE}','${MODIS_NEXT_FILE}'"
        elif [[  ${MODIS_NEXT_DATE} -gt ${END_DATE_DATA} ]]; then
           export NL_NUM_INPUT_FILES=2
           export NL_FILENAME_SEQ="'${MODIS_PAST_FILE}','${MODIS_PRES_FILE}'"
        else 
           export NL_NUM_INPUT_FILES=3
           export NL_FILENAME_SEQ="'${MODIS_PAST_FILE}','${MODIS_PRES_FILE}','${MODIS_NEXT_FILE}'"
        fi
        export NL_FILENAME_OUT="'obs_seq_${L_DATE}.out'"
        export NL_FIRST_OBS_DAYS=${ASIM_MIN_DAY_GREG}
        export NL_FIRST_OBS_SECONDS=${ASIM_MIN_SEC_GREG}
        export NL_LAST_OBS_DAYS=${ASIM_MAX_DAY_GREG}
        export NL_LAST_OBS_SECONDS=${ASIM_MAX_SEC_GREG}
        export NL_SYNONYMOUS_COPY_LIST="'MODIS_AOD_RETRIEVAL'"
        export NL_SYNONYMOUS_QC_LIST="'MODIS QC index'"
        ${HYBRID_SCRIPTS_DIR}/da_create_dart_input_nml.ksh       
#
# Make obs_def_apm_nml for apm_scale to adjust observation error variance
        rm -rf obs_def_apm.nml
        cat <<EOF > obs_def_apm.nml
&obs_def_apm_nml
apm_scale=${NL_APM_SCALE}
apm_scale_sw=${NL_APM_SCALE_SW}
/
EOF
#
        ./obs_sequence_tool
        if [[ -e obs_seq_${L_DATE}.out ]]; then
           cp obs_seq_${L_DATE}.out ${DATA_DIR}/obs_MODIS_AOD/obs_seq_modis_aod_${D_DATE}
        else
           touch ${DATA_DIR}/obs_MODIS_AOD/NO_OBS_SEQ.OUT_DATA_${D_DATE}
        fi
        cd ${ASIM_DIR}
     fi 
#
# LOOP TO NEXT DAY AND TIME 
     export P_DATE=${L_DATE}
     export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${TIME_INC} 2>/dev/null)  
  done 
exit
