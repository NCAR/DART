#!/bin/ksh -aeux 
#
  set echo
#
# TO SETUP AN ENVIRONMENT TO RUN obs_sequence_tool.
#
  export OBS_SEQ_EXT=out
  export OBS_DIR=obs_MOPCOMB
  export OBS_DIR=obs_MOPITT_Cx_Std_SVD
# SET TIME INFORMATION
  export START_DATE=2008060118
  export END_DATE=2008060500
  export END_DATE=2008060118
  export TIME_INC=6
  export ASIM_WINDOW=3
#
# USE SWITCHES
  export USE_HSI=true
  export USE_INPUT_DATA=true
#
# PATHS
  export WRFDA_VER=WRFDAv3.4_dmpar
  export WRF_VER=WRFv3.4_dmpar
  export DART_VER=DART_DEVEL_COM
#
# INDEPENDENT DIRECTORIES
  export CODE_DIR=/glade/p/work/mizzi/TRUNK
  export HSI_INPUT_DIR=/MIZZI/AVE_TEST_DATA
  export INPUT_DIR=/glade/scratch/mizzi/AVE_TEST_DATA
  export HSI_OUTPUT_DIR=/MIZZI/DART_TEST_AVE/MOPITT_ONLY
  export OUTPUT_DIR=/glade/scratch/mizzi/DART_TEST_AVE/MOPITT_ONLY
  export ASIM_DIR=/glade/scratch/mizzi/RUN_OBS_SEQUENCE_TOOL
#
# DEPENDENT DIRECTORIES
  export WRF_DIR=${CODE_DIR}/${WRF_VER}
  export VAR_DIR=${CODE_DIR}/${WRFDA_VER}
  export BUILD_DIR=${VAR_DIR}/var/build
  export DART_DIR=${CODE_DIR}/${DART_VER}
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
# RUN OBS_SEQUENCE_TOOL
     cd ${ASIM_DIR}
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
     if ${USE_INPUT_DATA}; then
        if ${USE_HSI}; then
           hsi get obs_seq.${OBS_SEQ_EXT} : ${HSI_INPUT_DIR}/${OBS_DIR}/${L_DATE}/obs_seq_${L_DATE}.${OBS_SEQ_EXT}
        else
           cp ${INPUT_DIR}/${OBS_DIR}/${L_DATE}/obs_seq_${L_DATE}.${OBS_SEQ_EXT} obs_seq.${OBS_SEQ_EXT}
        fi
     else
        if ${USE_HSI}; then
           hsi get obs_seq.${OBS_SEQ_EXT} : ${HSI_OUTPUT_DIR}/${L_DATE}/dart_filter/obs_seq.${OBS_SEQ_EXT}
        else
           cp ${OUTPUT_DIR}/${L_DATE}/dart_filter/obs_seq.${OBS_SEQ_EXT} obs_seq.${OBS_SEQ_EXT}
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
     export NL_NUM_INPUT_FILES=1
     export NL_FILENAME_SEQ="'"obs_seq.${OBS_SEQ_EXT}"'"
     export NL_FILENAME_OUT="'obs_seq.proc'"
     export NL_FIRST_OBS_DAYS=${ASIM_MIN_DAY_GREG}
     export NL_FIRST_OBS_SECONDS=${ASIM_MIN_SEC_GREG}
     export NL_LAST_OBS_DAYS=${ASIM_MAX_DAY_GREG}
     export NL_LAST_OBS_SECONDS=${ASIM_MAX_SEC_GREG}
     export NL_PRINT_ONLY=.false.
# min wrf
     export NL_MIN_LAT=24.
     export NL_MAX_LAT=36.
     export NL_MIN_LON=-154
     export NL_MAX_LON=-72.
# global
#     export NL_MIN_LAT=-90.
#     export NL_MAX_LAT=90.
#     export NL_MIN_LON=-180
#     export NL_MAX_LON=180.
     export NL_OBS_TYPES="'MOPITT_CO_RETRIEVAL'"
#     export NL_OBS_TYPES="''"
     export NL_KEEP_TYPES=.true.
     export NL_QC_METADATA="'DART quality control'"
     export NL_MIN_QC=0
     export NL_MAX_QC=2
     rm input.nml
     ${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_input.nml.ksh       
     ./obs_sequence_tool
     if ${USE_INPUT_DATA}; then
        if ${USE_HSI}; then
           hsi "cd ${HSI_INPUT_DIR}/${OBS_DIR}/${L_DATE}; put obs_seq.proc"
        else
           cp obs_seq.proc ${INPUT_DIR}/${OBS_DIR}/${L_DATE}/obs_seq.proc        
        fi
     else
        cp obs_seq.proc ${OUTPUT_DIR}/${L_DATE}/dart_filter/obs_seq.proc
     fi
#
# LOOP TO NEXT DAY AND TIME 
     export P_DATE=${L_DATE}
     export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${TIME_INC} 2>/dev/null)  
  done 
exit
