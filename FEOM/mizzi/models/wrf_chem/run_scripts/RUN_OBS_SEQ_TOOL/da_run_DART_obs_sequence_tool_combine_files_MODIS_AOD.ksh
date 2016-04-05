#!/bin/ksh -aeux 
#
  set echo
# 
# Script to combine multiple obs_seq files into a single obs_seq file
#
# SET TIME INFORMATION
  export START_DATE=2008060106
  export END_DATE=2008063018
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
  export ASIM_DIR=/glade/scratch/mizzi/MODIS_OBSSEQ_COMB
  export MET_OBS_DIR=${DATA_DIR}/obs_MET
  export RET_MODIS_AOD_OBS_DIR=${DATA_DIR}/obs_MODIA_AOD
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
#
# DART TIME INFO
     if [[ ${HH} -eq 0 ]]; then
        export DD_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} -1 2>/dev/null)  
        export DD_YYYY=$(echo $DD_DATE | cut -c1-4)
        export DD_YY=$(echo $DD_DATE | cut -c3-4)
        export DD_MM=$(echo $DD_DATE | cut -c5-6)
        export DD_DD=$(echo $DD_DATE | cut -c7-8)
        export DD_HH=24
     else
        export DD_YYYY=$(echo $L_DATE | cut -c1-4)
        export DD_YY=$(echo $L_DATE | cut -c3-4)
        export DD_MM=$(echo $L_DATE | cut -c5-6)
        export DD_DD=$(echo $L_DATE | cut -c7-8)
        export DD_HH=$(echo $L_DATE | cut -c9-10)
     fi
     export DT_DATE=${DD_YYYY}${DD_MM}${DD_DD}${DD_HH}
#    
# Use obs_sequence_tool to combine multiple obs_seq files
     if [[ ! -d ${ASIM_DIR}/mdart_to_sdart/${YYYY}${MM} ]]; then 
        mkdir -p ${ASIM_DIR}/mdart_to_sdart/${YYYY}${MM} 
     fi
     cd ${ASIM_DIR}/mdart_to_sdart/${YYYY}${MM}
#
     cp ${DART_DIR}/models/wrf_chem/work/advance_time ./.
     cp ${DART_DIR}/models/wrf_chem/work/obs_sequence_tool ./.
     cp ${DART_DIR}/models/wrf_chem/work/input.nml ./.
#
# GET OBS FILES TO BE COMBINED
     rm -rf ./obs_seq_MET
     rm -rf ./obs_seq_MODIS_AOD
     export MET_FLG=0
     if [[ -e ${MET_OBS_DIR}/${L_DATE}/obs_seq_${L_DATE}.out ]]; then
        export MET_FLG=1
        cp ${MET_OBS_DIR}/${L_DATE}/obs_seq_${L_DATE}.out ./obs_seq_MET
     fi
     export IAS_FLG=0
     if [[ -e ${RET_MODIS_AOD_OBS_DIR}/obs_seq_modis_aod_${DT_DATE} ]]; then
        export IAS_FLG=1
        cp ${RET_IASI_O3_OBS_DIR}/obs_seq_modis_aod_${DT_DATE} ./obs_seq_MODIS_AOD
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
     if [[ ${MET_FLG} -eq 1 && ${IAS_FLG} -eq 0 ]]; then
        export NL_NUM_INPUT_FILES=1
        export NL_FILENAME_SEQ="'obs_seq_MET'"
     elif [[ ${MET_FLG} -eq 0 ]]; then
        echo APM: No MET data
        exit
     elif [[ ${MET_FLG} -eq 1 && ${IAS_FLG} -eq 1 ]]; then
        export NL_NUM_INPUT_FILES=2
        export NL_FILENAME_SEQ="'obs_seq_MET','obs_seq_MODIS_AOD'"
     fi    
     export NL_FILENAME_OUT="'obs_seq.proc'"
     export NL_FIRST_OBS_DAYS=${ASIM_MIN_DAY_GREG}
     export NL_FIRST_OBS_SECONDS=${ASIM_MIN_SEC_GREG}
     export NL_LAST_OBS_DAYS=${ASIM_MAX_DAY_GREG}
     export NL_LAST_OBS_SECONDS=${ASIM_MAX_SEC_GREG}
     export NL_SYNONYMOUS_COPY_LIST="'NCEP BUFR observation','MODIS_AOD_RETRIEVAL'"
     export NL_SYNONYMOUS_QC_LIST="'NCEP QC index','MODIS QC index'"
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
#
     ./obs_sequence_tool
     mkdir -p ${DATA_DIR}/obs_MODCOMB_AOD_Mig_DA/${L_DATE}
     cp obs_seq.proc ${DATA_DIR}/obs_MODCOMB_AOD_Mig_DA/${L_DATE}/obs_seq_comb_${L_DATE}.out
     cd ${ASIM_DIR}
#
# LOOP TO NEXT DAY AND TIME 
     export P_DATE=${L_DATE}
     export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${TIME_INC} 2>/dev/null)  
  done 
exit
