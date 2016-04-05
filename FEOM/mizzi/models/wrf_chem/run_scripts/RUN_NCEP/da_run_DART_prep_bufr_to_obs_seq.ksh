#!/bin/ksh -aeux 
#
# set switches
  export DO_PREP_BUFR=true
  export DO_CREATE_ASCII_TO_DART=true
#
# set time data
#  export START_DATE=2008060106
#  export END_DATE=2008063000
  export START_DATE=2008070106
  export START_DATE=2008072906
  export END_DATE=2008073118
  export TIME_INC=6
#
# system settings
  export PROCS=8
  export USE_HSI=false
#
# set paths
  export WRFDA_VER=WRFDAv3.4_dmpar
  export WRF_VER=WRFv3.4_dmpar
  export DART_VER=DART_CHEM
#
# independent directories
  export CODE_DIR=/glade/p/work/mizzi/TRUNK
  export DATA_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA/obs_MET
  export HSI_DATA_DIR=/MIZZI/AVE_TEST_DATA/obs_MET
  export ASIM_DIR=/glade/scratch/mizzi/dart_assim
  export RUN_ROOT=/glade/scratch/mizzi
#
# dependent directories
  export HYBRID_DIR=/glade/p/work/mizzi/HYBRID_TRUNK
  export WRF_DIR=${CODE_DIR}/${WRF_VER}
  export VAR_DIR=${CODE_DIR}/${WRFDA_VER}
  export BUILD_DIR=${VAR_DIR}/var/build
  export DART_DIR=${CODE_DIR}/${DART_VER}
  export TOOL_DIR=${VAR_DIR}/var/da
  export ICBC_DIR=${ASIM_DIR}
  export OBS_DIR=${ASIM_DIR}
  export HYBRID_SCRIPTS_DIR=${HYBRID_DIR}/hybrid_scripts
#
# loop over IC times
  export L_DATE=${START_DATE}
  while [[ ${L_DATE} -le ${END_DATE} ]]; do
     export DATE=${L_DATE}
     export YYYY=$(echo $DATE | cut -c1-4)
     export YY=$(echo $DATE | cut -c3-4)
     export MM=$(echo $DATE | cut -c5-6)
     export DD=$(echo $DATE | cut -c7-8)
     export HH=$(echo $DATE | cut -c9-10)
     export FILE_DATE=${YYYY}-${MM}-${DD}_${HH}:00:00
     export PAST_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} -${TIME_INC} 2>/dev/null)  
     export PAST_YYYY=$(echo $PAST_DATE | cut -c1-4)
     export PAST_YY=$(echo $PAST_DATE | cut -c3-4)
     export PAST_MM=$(echo $PAST_DATE | cut -c5-6)
     export PAST_DD=$(echo $PAST_DATE | cut -c7-8)
     export PAST_HH=$(echo $PAST_DATE | cut -c9-10)
#
# DART time data (no leading zeros)
     export DT_YYYY=${YYYY}
     export DT_YY=$(echo $DATE | cut -c3-4)
     export DT_MM=${MM} 
     export DT_DD=${DD} 
     export DT_HH=${HH} 
     (( DT_MM = ${DT_MM} + 0 ))
     (( DT_DD = ${DT_DD} + 0 ))
     (( DT_HH = ${DT_HH} + 0 ))
#    
# make asimilation directory and go to it
     if ${DO_PREP_BUFR}; then
        cd ${RUN_ROOT}        
        if [[ -d ${ASIM_DIR}/prep_bufr ]]; then 
           cd ${ASIM_DIR}/prep_bufr
           rm -rf *
        else
           mkdir -p ${ASIM_DIR}/prep_bufr
           cd ${ASIM_DIR}/prep_bufr
        fi
#
# get prepbufr files
        export LL_DATE=${L_DATE}
        export LL_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${LL_DATE} +18 2>/dev/null)  
        while [[ ${LL_DATE} -le ${LL_END_DATE} ]]; do
           export LL_YYYY=$(echo $LL_DATE | cut -c1-4)
           export LL_YY=$(echo $LL_DATE | cut -c3-4)
           export LL_MM=$(echo $LL_DATE | cut -c5-6)
           export LL_DD=$(echo $LL_DATE | cut -c7-8)
           export LL_HH=$(echo $LL_DATE | cut -c9-10)
           if ${USE_HSI}; then
              hsi get prepqm${LL_YY}${LL_MM}${LL_DD}${LL_HH} : ${HSI_DATA_DIR}/${LL_YYYY}${LL_MM}${LL_DD}${LL_HH}/prepbufr.gdas.${LL_YYYY}${LL_MM}${LL_DD}${LL_HH}.wo40.be
           else
             cp ${DATA_DIR}/${LL_YYYY}${LL_MM}${LL_DD}${LL_HH}/prepbufr.gdas.${LL_YYYY}${LL_MM}${LL_DD}${LL_HH}.wo40.be prepqm${LL_YY}${LL_MM}${LL_DD}${LL_HH}
           fi  
           export P_DATE=${LL_DATE}
           export LL_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${TIME_INC} 2>/dev/null)  
        done
#
# get DART input.nml
        rm -rf input.nml
        cp ${DART_DIR}/observations/NCEP/prep_bufr/work/input.nml ./.
#
# run prepbufr to ascii converter
        ${DART_DIR}/observations/NCEP/prep_bufr/work/prepbufr.csh ${DT_YYYY} ${DT_MM} ${DT_DD} ${DT_DD} ${DART_DIR}/observations/NCEP/prep_bufr/exe
        cd ${RUN_ROOT}
     fi
#
     if ${DO_CREATE_ASCII_TO_DART}; then 
        cd ${RUN_ROOT}
        if [[ -d ${ASIM_DIR}/ascii_to_obs ]]; then 
           cd ${ASIM_DIR}/ascii_to_obs
           rm -rf *
        else
           mkdir -p ${ASIM_DIR}/ascii_to_obs
           cd ${ASIM_DIR}/ascii_to_obs
        fi
#
# set DART time data
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
# run ascii to obs_seq converter
        ${HYBRID_SCRIPTS_DIR}/da_create_dart_ncep_ascii_to_obs_input_nml.ksh
        ${DART_DIR}/observations/NCEP/ascii_to_obs/work/create_real_obs
#
# save obs_seq files
        export LL_DATE=${L_DATE}
        export LL_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${LL_DATE} +18 2>/dev/null)  
        while [[ ${LL_DATE} -le ${LL_END_DATE} ]]; do
           export LL_YYYY=$(echo $LL_DATE | cut -c1-4)
           export LL_YY=$(echo $LL_DATE | cut -c3-4)
           export LL_MM=$(echo $LL_DATE | cut -c5-6)
           export LL_DD=$(echo $LL_DATE | cut -c7-8)
           export LL_HH=$(echo $LL_DATE | cut -c9-10)
           export LL_PAST_DATE=$(${BUILD_DIR}/da_advance_time.exe ${LL_DATE} -${TIME_INC} 2>/dev/null)  
           export LL_PAST_YYYY=$(echo $PAST_DATE | cut -c1-4)
           export LL_PAST_YY=$(echo $PAST_DATE | cut -c3-4)
           export LL_PAST_MM=$(echo $PAST_DATE | cut -c5-6)
           export LL_PAST_DD=$(echo $PAST_DATE | cut -c7-8)
           export LL_PAST_HH=$(echo $PAST_DATE | cut -c9-10)
           if [[ ${LL_HH} -eq 0 ]] then
              export L_YYYY=${LL_PAST_YYYY}
              export L_MM=${LL_PAST_MM}
              export L_DD=${LL_PAST_DD}
              export L_HH=24
              export D_DATE=${L_YYYY}${L_MM}${L_DD}${L_HH}
           else
              export L_YYYY=${LL_YYYY}
              export L_MM=${LL_MM}
              export L_DD=${LL_DD}
              export L_HH=${LL_HH}
              export D_DATE=${L_YYYY}${L_MM}${L_DD}${L_HH}
           fi
           if [[ ! -d ${DATA_DIR}/${LL_DATE} ]]; then
              mkdir -p ${DATA_DIR}/${LL_DATE}
           fi
           cp -f obs_seq${D_DATE} ${DATA_DIR}/${LL_DATE}/obs_seq_${LL_DATE}.out
           export P_DATE=${LL_DATE}
           export LL_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${TIME_INC} 2>/dev/null)  
        done
        cd ${RUN_ROOT} 
     fi
     export P_DATE=${L_DATE}
     export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} 24 2>/dev/null)  
  done 
