#!/bin/ksh -aeux 
#
set echo
# 
# Script to combine multiple obs_seq files into a single obs_seq file
#
# SET TIME INFORMATION
  export START_DATE=2008060106
#  export START_DATE=2008060112
  export END_DATE=2008063018
#  export END_DATE=2008060112
  export TIME_INC=6
  export ASIM_WINDOW=3
#
# SYSTEM SPECIFIC SETTINGS
  export PROCS=8
  export NL_APM_SCALE=1.
  export NL_APM_SCALE_SW=.FALSE.
#
# PATHS
  export WRFDA_VER=WRFDAv3.4_dmpar
  export WRF_VER=WRFv3.4_dmpar
  export DART_VER=DART_CHEM
#
# INDEPENDENT DIRECTORIES
  export ROOT_DIR=/glade/p/work/mizzi
  export CODE_DIR=/glade/p/work/mizzi/TRUNK
  export DATA_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA
  export ASIM_DIR=/glade/scratch/mizzi/MET_MOP_IAS_COMB
  export RET_MOPITT_OBS_DIR=${DATA_DIR}/obs_MOPITT_CO_DnN_Mig_DA_DBL_bloc
  export RET_IASI_OBS_DIR=${DATA_DIR}/obs_IASI_CO_DnN_Mig_DA_DBL_bloc
  export WRITE_OUT_NAME=obs_MOPnIAS_CO_Mig_DA_DBL_bloc
#
# DEPENDENT DIRECTORIES
  export HYBRID_DIR=${ROOT_DIR}/HYBRID_TRUNK
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
  cd ${ASIM_DIR}
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
# Use obs_sequence_tool to combine multiple obs_seq files
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
        export MET_FLG=0
        export MOP_FLG=0
        export IAS_FLG=0
        if [[ -e ${DATA_DIR}/obs_MET/${L_DATE}/obs_seq_${L_DATE}.out ]]; then 
           export MET_FLG=1
           cp ${DATA_DIR}/obs_MET/${L_DATE}/obs_seq_${L_DATE}.out ./obs_seq_MET_${L_DATE}.out
        fi
        if [[ -e ${RET_MOPITT_OBS_DIR}/obs_seq_mopitt_${D_DATE} ]];  then
           export MOP_FLG=1
           cp ${RET_MOPITT_OBS_DIR}/obs_seq_mopitt_${D_DATE} ./obs_seq_MOP_${L_DATE}.out
        fi
        if [[ -e ${RET_IASI_OBS_DIR}/obs_seq_iasi_${D_DATE} ]]; then
           export IAS_FLG=1
           cp ${RET_IASI_OBS_DIR}/obs_seq_iasi_${D_DATE} ./obs_seq_IAS_${L_DATE}.out
        fi
        export FILE_FLG=2
        if [[ ${MET_FLG} -eq 1 && ${MOP_FLG} -eq 1 && ${IAS_FLG} -eq 1 ]]; then
           export FILE_FLG=3
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
# min wrf
        if [[ ${FILE_FLG} -eq 2 ]]; then
           export NL_NUM_INPUT_FILES=2
           export NL_FILENAME_SEQ="'obs_seq_MET_${L_DATE}.out','obs_seq_MOP_${L_DATE}.out'"
        fi
        if [[ ${FILE_FLG} -eq 3 ]]; then
           export NL_NUM_INPUT_FILES=3
           export NL_FILENAME_SEQ="'obs_seq_MET_${L_DATE}.out','obs_seq_MOP_${L_DATE}.out','obs_seq_IAS_${L_DATE}.out'"
        fi
        export NL_FILENAME_OUT="'obs_seq.proc'"
        export NL_FIRST_OBS_DAYS=${ASIM_MIN_DAY_GREG}
        export NL_FIRST_OBS_SECONDS=${ASIM_MIN_SEC_GREG}
        export NL_LAST_OBS_DAYS=${ASIM_MAX_DAY_GREG}
        export NL_LAST_OBS_SECONDS=${ASIM_MAX_SEC_GREG}
        export NL_SYNONYMOUS_COPY_LIST="'NCEP BUFR observation','MOPITT CO observation','IASI CO observation'"
        export NL_SYNONYMOUS_QC_LIST="'NCEP QC index','MOPITT CO QC index','IASI CO QC index'"
        export NL_MIN_LAT=7.
        export NL_MAX_LAT=54.
        export NL_MIN_LON=184.
        export NL_MAX_LON=310.
        rm input.nml
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
        ./obs_sequence_tool
#        if [[ ! -d ${DATA_DIR}/obs_MOPITT/${L_DATE} ]]; then mkdir -p ${DATA_DIR}/obs_MOPITT/${L_DATE}; fi
        mkdir -p ${DATA_DIR}/${WRITE_OUT_NAME}/${L_DATE}
        cp obs_seq.proc ${DATA_DIR}/${WRITE_OUT_NAME}/${L_DATE}/obs_seq_comb_${L_DATE}.out
        cd ${ASIM_DIR}
#
# LOOP TO NEXT DAY AND TIME 
     export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${TIME_INC} 2>/dev/null)  
  done 
exit
