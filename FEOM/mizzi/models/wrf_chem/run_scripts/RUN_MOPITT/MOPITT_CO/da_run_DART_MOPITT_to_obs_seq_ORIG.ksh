#!/bin/ksh -aeux 
#
  set echo
#
# MODIFIED VERSION OF /DART/models/WRF/regression/CONUS-V3/icbc_real.ksh
# TO SETUP AN ENVIRONMENT TO CONVERT OBSERVATIONS TO obs_seq.
#
# SET TIME INFORMATION
  export START_DATE=2008060106
  export END_DATE=2008070100
  export TIME_INC=6
  export ASIM_WINDOW=3
#
# SYSTEM SPECIFIC SETTINGS
  export PROCS=8
  export OB_TYPE=obs
#
# PATHS
  export WRFDA_VER=WRFDAv3.4_dmpar
  export WRF_VER=WRFv3.4_dmpar
  export DART_VER=DART_DEVEL_COM
#
# INDEPENDENT DIRECTORIES
  export CODE_DIR=/glade/p/work/mizzi/TRUNK
  export DATA_DIR=/glade/scratch/mizzi/AVE_TEST_DATA
  export ASIM_DIR=/glade/scratch/mizzi/MOPITT_to_OBSSEQ
  export OBS_MOPITT_DIR=/glade/p/acd/mizzi/for_arthur/MOPITT/proc_files
#
# DEPENDENT DIRECTORIES
  export HYBRID_DIR=${CODE_DIR}/HYBRID_TRUNK
  export WRF_DIR=${CODE_DIR}/${WRF_VER}
  export VAR_DIR=${CODE_DIR}/${WRFDA_VER}
  export BUILD_DIR=${VAR_DIR}/var/build
  export DART_DIR=${CODE_DIR}/${DART_VER}
  export TOOL_DIR=${VAR_DIR}/var/da
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
# RUN_MOPITT_ASCII_TO_DART
        cd ${ASIM_DIR}
        if [[ ! -d ${ASIM_DIR}/MOPITT_ascii_to_dart/${YYYY}${MM} ]]; then mkdir -p ${ASIM_DIR}/MOPITT_ascii_to_dart/${YYYY}${MM}; fi
        cd ${ASIM_DIR}/MOPITT_ascii_to_dart/${YYYY}${MM}
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
        export NL_YEAR=${L_YYYY}
        export NL_MONTH=${L_MM}
        export NL_DAY=${L_DD}
        export NL_HOUR=${L_HH}
        if [[ ${L_HH} -eq 24 ]]; then
           NL_BIN_BEG=21.01
           NL_BIN_END=3.00
        elif [[ ${L_HH} -eq 6 ]]; then
           NL_BIN_BEG=3.01
           NL_BIN_END=9.00
        elif [[ ${L_HH} -eq 12 ]]; then
           NL_BIN_BEG=9.01
           NL_BIN_END=15.00
        elif [[ ${L_HH} -eq 18 ]]; then
           NL_BIN_BEG=15.01
           NL_BIN_END=21.00
        fi
        export NL_FILEDIR="'"${ASIM_DIR}/MOPITT_ascii_to_dart/${YYYY}${MM}/"'" 
        export NL_FILENAME="'"${D_DATE}.dat"'" 
#
# USE MOPITT DATA 
        if [[ -e input.nml ]]; then rm -rf input.nml; fi
        ${HYBRID_SCRIPTS_DIR}/da_create_dart_mopitt_input_nml.ksh
        cp ${OBS_MOPITT_DIR}/${D_DATE}.dat ./.

        cp ${DART_DIR}/observations/MOPITT/work/mopitt_ascii_to_obs ./.
        ./mopitt_ascii_to_obs
        export MOPITT_FILE=mopitt_obs_seq${D_DATE}
        cp ${MOPITT_FILE} ${OBS_MOPITT_DIR}/obs_seq_mopitt_${D_DATE}
#
# LOOP TO NEXT DAY AND TIME 
     export P_DATE=${L_DATE}
     export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${TIME_INC} 2>/dev/null)  
  done 
exit
