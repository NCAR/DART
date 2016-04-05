#!/bin/ksh -x
###############################################################################
#
#  Script to run obs_diag for WRFCHEM in the DART framework
#
############################################################################### 
#
# Define experiment parameters
export START_DATE=2008061006
export START_DATE=2008060112
export END_DATE=2008061500
export END_DATE=2008062900
#export END_DATE=2008061118
#
#export START_DATE=2008060112
#export END_DATE=2008060200
#
export NL_APM_SCALE=1.
export NL_APM_SCALE_SW=.FALSE.
export DOMAIN=01
export NUM_MEMBERS=20
export CYCLE_PERIOD=6
export FCST_PERIOD=24
export FCST_PERIOD=6
export ASIM_PERIOD=3
export LBC_FREQ=3
(( INTERVAL_SEC=${LBC_FREQ}*60*60 ))
(( CYCLE_PERIOD_SEC=${CYCLE_PERIOD}*60*60 ))
#
# Define use options
export USE_HSI=false
export DELETE_FLG=false
#
# Define code versions
export DART_VER=DART_CHEM
export WRFCHEM_VER=WRFCHEMv3.4_dmpar
export WRF_VER=WRFv3.4_dmpar
export WRFDA_VER=WRFDAv3.4_dmpar
#
# Set job submission parameters
export PROJ_NUMBER=P19010000
export TIME_LIMIT_FILTER=1:40
export TIME_LIMIT_WRFCHEM=1:40
export NUM_TASKS=32
export TASKS_PER_NODE=16
export JOB_CLASS=small
#
# Define independent directory paths
#export DIR_NAME=MOPCOMB_Exp_2_RtDA_40M_p30p30_sp4
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_40M_p30p30_sp4
#export DIR_NAME=MOPCOMB_Exp_3_MgDA_40M_p30p30_sp4
#
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_p30p30
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_p10p00
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_p20p00
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_p30p00
#export DIR_NAME=MOPCOMB_Exp_3_MgDA_20M_p10p00
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_p10p30_loc
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_p10p30
#
#export DIR_NAME=MOPCOMB_Exp_3_MgDA_20M_100km_p10p00
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_loc_p10p30
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_loc_a_p10p30
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_loc_b_p10p30
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_bar_p10p30
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_bar_1_p10p30
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_bar_2_p10p30
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_NoRot_p10p30
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_NoRot_loc_p10p30
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_NoRot_bloc_p10p30
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_bloc_p10p30
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_TEST_p75
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_TEST_p50
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_TEST_p60
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_TEST_p55
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_DBL_p10p30
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_DBL_NV_p10p30
#export DIR_NAME=MOPCOMB_Exp_2_MgDA_20M_100km_DBL_bloc_p10p30
export DIR_NAME=MOPCOMB_Exp_2_Mg_XXXnIAS_20M_100km_p10p30
export DIR_NAME=MOPCOMB_Exp_2_Mg_MOPnIAS_20M_100km_p10p30
export DIR_NAME=MOPCOMB_Exp_2_Mg_MOPnXXX_20M_100km_p10p30
export DIR_NAME=MOPCOMB_Exp_2_MgDA_ROE_p75_20M_100km_p10p30
export DIR_NAME=IASCOMB_O3_Exp_2_MgDA_20M_100km_p10p30/
#
#
# Independent path settings
export SCRATCH_DIR=/glade/scratch/mizzi
export PROJECT_DIR=/glade/p/work/mizzi
export ACD_DIR=/glade/p/acd/mizzi
export HSI_DIR=/MIZZI
#
# Dependent path settings
export RUN_DIR=${ACD_DIR}/DART_OBS_DIAG/${DIR_NAME}
export EXP_DIR=${RUN_DIR}
export TRUNK_DIR=${PROJECT_DIR}/TRUNK
export DATA_DIR=${ACD_DIR}/AVE_TEST_DATA
export HSI_DATA_DIR=${HSI_DIR}/AVE_TEST_DATA
export SAVE_DIR=${SCRATCH_DIR}/DART_TEST_AVE/${DIR_NAME}
#export SAVE_DIR=${ACD_DIR}/DART_TEST_AVE/${DIR_NAME}
export HSI_SAVE_DIR=${HSI_DIR}/DART_TEST_AVE/${DIR_NAME}
#
export DART_DIR=${TRUNK_DIR}/${DART_VER}
export WRF_DIR=${TRUNK_DIR}/${WRF_VER}
export WRFCHEM_DIR=${TRUNK_DIR}/${WRFCHEM_VER}
export WRFDA_DIR=${TRUNK_DIR}/${WRFDA_VER}/var
#
# Copy necessary executables from DART to $RUN_DIR
if [[ ! -d ${RUN_DIR} ]]; then mkdir -p ${RUN_DIR}; fi
cd ${RUN_DIR}
cp ${DART_DIR}/models/wrf_chem/work/input.nml ./.
cp ${DART_DIR}/models/wrf_chem/work/advance_time ./.
cp ${DART_DIR}/models/wrf_chem/work/obs_diag ./.
#
# Build obs_seq.final file list
cd ${RUN_DIR}
rm -rf file_list.txt
export L_DATE=${START_DATE}
while [[ ${L_DATE} -le ${END_DATE} ]]; do
#
# Set date/time information
   export L_YY=`echo $L_DATE | cut -c1-4`
   export L_MM=`echo $L_DATE | cut -c5-6`
   export L_DD=`echo $L_DATE | cut -c7-8`
   export L_HH=`echo $L_DATE | cut -c9-10`
   export L_FILE_DATE=${L_YY}-${L_MM}-${L_DD}_${L_HH}:00:00
#
   export NEXT_DATE=`echo ${L_DATE} +${FCST_PERIOD}h | ./advance_time` 
   export NEXT_YY=`echo $NEXT_DATE | cut -c1-4`
   export NEXT_MM=`echo $NEXT_DATE | cut -c5-6`
   export NEXT_DD=`echo $NEXT_DATE | cut -c7-8`
   export NEXT_HH=`echo $NEXT_DATE | cut -c9-10`
   export NEXT_FILE_DATE=${NEXT_YY}-${NEXT_MM}-${NEXT_DD}_${NEXT_HH}:00:00
#
# Create obs_seq file list
   export FILE=obs_seq.final
   if [[ ! -d ${EXP_DIR}/${L_DATE}/dart_filter ]]; then
      mkdir -p ${EXP_DIR}/${L_DATE}/dart_filter
      cd ${EXP_DIR}/${L_DATE}/dart_filter
   else
      cd ${EXP_DIR}/${L_DATE}/dart_filter
   fi
   if [[ -f ${FILE} && ${DELETE_FLG} == true ]]; then
      rm -rf ${FILE}
   fi
   if [[ ! -f ${FILE} ]]; then
      if ${USE_HSI}; then
         hsi get ${FILE} : ${HSI_SAVE_DIR}/${L_DATE}/dart_filter/${FILE}
      else
         cp ${SAVE_DIR}/${L_DATE}/dart_filter/${FILE} ${FILE}
      fi
   fi
   cd ${RUN_DIR}
   if [[ -f ${EXP_DIR}/${L_DATE}/dart_filter/${FILE} ]]; then
      echo ${EXP_DIR}/${L_DATE}/dart_filter/${FILE} >> file_list.txt
   else
      echo APM: hsi get or cp failed for ${HSI_EXP_DIR}/${L_DATE}/dart_filter/${FILE}
      exit
   fi
#
# Loop to next cycle time   
   export L_DATE=${NEXT_DATE}
done
cd ${RUN_DIR}
#
###############################################################################
#
# CREATE DART NAMELIST
#
###############################################################################
export STR_YY=`echo $START_DATE | cut -c1-4`
export STR_MM=`echo $START_DATE | cut -c5-6`
export STR_DD=`echo $START_DATE | cut -c7-8`
export STR_HH=`echo $START_DATE | cut -c9-10`
#
export END_YY=`echo $END_DATE | cut -c1-4`
export END_MM=`echo $END_DATE | cut -c5-6`
export END_DD=`echo $END_DATE | cut -c7-8`
export END_HH=`echo $END_DATE | cut -c9-10`
#
export ASIM_MIN_DATE_STR=`echo ${START_DATE} -${ASIM_PERIOD}h | ./advance_time` 
export ASIM_MIN_YY_STR=`echo $ASIM_MIN_DATE_STR | cut -c1-4`
export ASIM_MIN_MM_STR=`echo $ASIM_MIN_DATE_STR | cut -c5-6`
export ASIM_MIN_DD_STR=`echo $ASIM_MIN_DATE_STR | cut -c7-8`
export ASIM_MIN_HH_STR=`echo $ASIM_MIN_DATE_STR | cut -c9-10`
export ASIM_MAX_DATE_STR=`echo ${START_DATE} +${ASIM_PERIOD}h | ./advance_time` 
export ASIM_MAX_YY_STR=`echo $ASIM_MAX_DATE_STR | cut -c1-4`
export ASIM_MAX_MM_STR=`echo $ASIM_MAX_DATE_STR | cut -c5-6`
export ASIM_MAX_DD_STR=`echo $ASIM_MAX_DATE_STR | cut -c7-8`
export ASIM_MAX_HH_STR=`echo $ASIM_MAX_DATE_STR | cut -c9-10`
#
export ASIM_MIN_DATE_END=`echo ${END_DATE} -${ASIM_PERIOD}h | ./advance_time` 
export ASIM_MIN_YY_END=`echo $ASIM_MIN_DATE_END | cut -c1-4`
export ASIM_MIN_MM_END=`echo $ASIM_MIN_DATE_END | cut -c5-6`
export ASIM_MIN_DD_END=`echo $ASIM_MIN_DATE_END | cut -c7-8`
export ASIM_MIN_HH_END=`echo $ASIM_MIN_DATE_END | cut -c9-10`
export ASIM_MAX_DATE_END=`echo ${END_DATE} +${ASIM_PERIOD}h | ./advance_time` 
export ASIM_MAX_YY_END=`echo $ASIM_MAX_DATE_END | cut -c1-4`
export ASIM_MAX_MM_END=`echo $ASIM_MAX_DATE_END | cut -c5-6`
export ASIM_MAX_DD_END=`echo $ASIM_MAX_DATE_END | cut -c7-8`
export ASIM_MAX_HH_END=`echo $ASIM_MAX_DATE_END | cut -c9-10`

(( STR_MM = ${STR_MM} + 0 ))
(( STR_DD = ${STR_DD} + 0 ))
(( STR_HH = ${STR_HH} + 0 ))     
(( END_MM = ${END_MM} + 0 ))
(( END_DD = ${END_DD} + 0 ))
(( END_HH = ${END_HH} + 0 ))     
(( ASIM_MIN_MM_STR = ${ASIM_MIN_MM_STR} + 0 ))
(( ASIM_MIN_DD_STR = ${ASIM_MIN_DD_STR} + 0 ))
(( ASIM_MIN_HH_STR = ${ASIM_MIN_HH_STR} + 0 ))
(( ASIM_MAX_MM_STR = ${ASIM_MAX_MM_STR} + 0 ))
(( ASIM_MAX_DD_STR = ${ASIM_MAX_DD_STR} + 0 ))
(( ASIM_MAX_HH_STR = ${ASIM_MAX_HH_STR} + 0 ))
(( ASIM_MIN_MM_END = ${ASIM_MIN_MM_END} + 0 ))
(( ASIM_MIN_DD_END = ${ASIM_MIN_DD_END} + 0 ))
(( ASIM_MIN_HH_END = ${ASIM_MIN_HH_END} + 0 ))
(( ASIM_MAX_MM_END = ${ASIM_MAX_MM_END} + 0 ))
(( ASIM_MAX_DD_END = ${ASIM_MAX_DD_END} + 0 ))
(( ASIM_MAX_HH_END = ${ASIM_MAX_HH_END} + 0 ))
#
# &obs_diag_nml
export NL_OBS_SEQUENCE_NAME="''"
export NL_OBS_SEQUENCE_LIST="'file_list.txt'"
export NL_FIRST_BIN_CENTER_YY=${STR_YY}
export NL_FIRST_BIN_CENTER_MM=${STR_MM}
export NL_FIRST_BIN_CENTER_DD=${STR_DD}
export NL_FIRST_BIN_CENTER_HH=${STR_HH}
export NL_FIRST_BIN_CENTER_MN=0
export NL_FIRST_BIN_CENTER_SS=0
export NL_LAST_BIN_CENTER_YY=${END_YY}
export NL_LAST_BIN_CENTER_MM=${END_MM}
export NL_LAST_BIN_CENTER_DD=${END_DD}
export NL_LAST_BIN_CENTER_HH=${END_HH}
export NL_LAST_BIN_CENTER_MN=0
export NL_LAST_BIN_CENTER_SS=0
export NL_BIN_SEPARATION_YY=0
export NL_BIN_SEPARATION_MM=0
export NL_BIN_SEPARATION_DD=0
export NL_BIN_SEPARATION_HH=6
export NL_BIN_SEPARATION_MN=0
export NL_BIN_SEPARATION_SS=0
export NL_BIN_WIDTH_YY=0
export NL_BIN_WIDTH_MM=0
export NL_BIN_WIDTH_DD=0
export NL_BIN_WIDTH_HH=6
export NL_BIN_WIDTH_MN=0
export NL_BIN_WIDTH_SS=0
export NL_TIME_TO_SKIP_YY=0
export NL_TIME_TO_SKIP_MM=0
export NL_TIME_TO_SKIP_DD=0
export NL_TIME_TO_SKIP_HH=0
export NL_TIME_TO_SKIP_MN=0
export NL_TIME_TO_SKIP_SS=0
export NL_MAX_NUM_BINS=1000
export NL_PLEVEL='1000., 850., 750., 650., 550., 450., 350., 250., 150., 75.'
export NL_NREGIONS=3
export NL_LONLIM1=239.,284.,0.
export NL_LONLIM2=247.,303.,360.
export NL_LATLIM1=24.,32.,20.
export NL_LATLIM2=34.,46.,80.
export NL_REG_NAMES="'Los Angles','New York','CONUS'"
export NL_PRINT_MISMATCHED_LOCS=.false.
export NL_PRINT_OBS_LOCATIONS=.false.
export NL_VERBOSE=.true.
#
# &utilities_nml
export NL_TERMLEVEL=2
#
# &schedule_nml
export NL_CALENDAR="'Gregorian'"
export NL_FIRST_BIN_START_YY=${ASIM_MIN_YY_STR}
export NL_FIRST_BIN_START_MM=${ASIM_MIN_MM_STR}
export NL_FIRST_BIN_START_DD=${ASIM_MIN_DD_STR}
export NL_FIRST_BIN_START_HH=${ASIM_MIN_HH_STR}
export NL_FIRST_BIN_START_MN=0
export NL_FIRST_BIN_START_SS=0
export NL_FIRST_BIN_END_YY=${ASIM_MAX_YY_STR}
export NL_FIRST_BIN_END_MM=${ASIM_MAX_MM_STR}
export NL_FIRST_BIN_END_DD=${ASIM_MAX_DD_STR}
export NL_FIRST_BIN_END_HH=${ASIM_MAX_HH_STR}
export NL_FIRST_BIN_END_MN=0
export NL_FIRST_BIN_END_SS=0
export NL_LAST_BIN_START_YY=${ASIM_MIN_YY_END}
export NL_LAST_BIN_START_MM=${ASIM_MIN_MM_END}
export NL_LAST_BIN_START_DD=${ASIM_MIN_DD_END}
export NL_LAST_BIN_START_HH=${ASIM_MIN_HH_END}
export NL_LAST_BIN_START_MN=0
export NL_LAST_BIN_START_SS=0
export NL_LAST_BIN_END_YY=${ASIM_MAX_YY_END}
export NL_LAST_BIN_END_MM=${ASIM_MAX_MM_END}
export NL_LAST_BIN_END_DD=${ASIM_MAX_DD_END}
export NL_LAST_BIN_END_HH=${ASIM_MAX_HH_END}
export NL_LAST_BIN_END_MN=0
export NL_LAST_BIN_END_SS=0
export NL_BIN_INTERVAL_DAYS=0
export NL_BIN_INTERVAL_SECONDS=21600
export NL_MAX_NUMBER_BINS=1000
export NL_PRINT_TABLE=.false.
#
# &assim_tools_nml
   export NL_CUTOFF=0.1
   export NL_SPECIAL_LOCALIZATION_OBS_TYPES="'MOPITT_CO_RETRIEVAL'"
   export NL_SPECIAL_LOCALIZATION_CUTOFFS=0.025
#
# &ensemble_manager_nml
   export NL_SINGLE_RESTART_FILE_IN=.false.       
   export NL_SINGLE_RESTART_FILE_OUT=.false.       
#
# &assim_model_nml
#
# &model_nml
   export NL_DEFAULT_STATE_VARIABLES=.false.
   export NL_WRF_STATE_VARIABLES="'U',     'KIND_U_WIND_COMPONENT',     'TYPE_U',  'UPDATE','999',
                              'V',     'KIND_V_WIND_COMPONENT',     'TYPE_V',  'UPDATE','999',
                              'W',     'KIND_VERTICAL_VELOCITY',    'TYPE_W',  'UPDATE','999',
                              'PH',    'KIND_GEOPOTENTIAL_HEIGHT',  'TYPE_GZ', 'UPDATE','999',
                              'T',     'KIND_POTENTIAL_TEMPERATURE','TYPE_T',  'UPDATE','999',
                              'MU',    'KIND_PRESSURE',             'TYPE_MU', 'UPDATE','999',
                              'QVAPOR','KIND_VAPOR_MIXING_RATIO',   'TYPE_QV', 'UPDATE','999',
                              'QRAIN', 'KIND_RAINWATER_MIXING_RATIO','TYPE_QRAIN', 'UPDATE','999',
                              'QCLOUD','KIND_CLOUD_LIQUID_WATER',   'TYPE_QCLOUD', 'UPDATE','999',
                              'QSNOW', 'KIND_SNOW_MIXING_RATIO',    'TYPE_QSNOW', 'UPDATE','999',
                              'QICE',  'KIND_CLOUD_ICE',            'TYPE_QICE', 'UPDATE','999',
                              'U10',   'KIND_U_WIND_COMPONENT',     'TYPE_U10','UPDATE','999',
                              'V10',   'KIND_V_WIND_COMPONENT',     'TYPE_V10','UPDATE','999',
                              'T2',    'KIND_TEMPERATURE',          'TYPE_T2', 'UPDATE','999',
                              'TH2',   'KIND_POTENTIAL_TEMPERATURE','TYPE_TH2','UPDATE','999',
                              'Q2',    'KIND_SPECIFIC_HUMIDITY',    'TYPE_Q2', 'UPDATE','999',
                              'PSFC',  'KIND_PRESSURE',             'TYPE_PS', 'UPDATE','999',
                              'o3',    'KIND_O3',                   'TYPE_O3', 'UPDATE','999',
                              'co',    'KIND_CO',                   'TYPE_CO', 'UPDATE','999',
                              'no',    'KIND_NO',                   'TYPE_NO', 'UPDATE','999',
                              'no2',   'KIND_NO2',                  'TYPE_NO2', 'UPDATE','999',
                              'hno3',  'KIND_HNO3',                 'TYPE_HNO3', 'UPDATE','999',
                              'hno4',  'KIND_HNO4',                 'TYPE_HNO4', 'UPDATE','999',
                              'n2o5',  'KIND_N2O5',                 'TYPE_N2O5', 'UPDATE','999',
                              'pan',   'KIND_PAN',                  'TYPE_PAN', 'UPDATE','999',
                              'mek',   'KIND_MEK',                  'TYPE_MEK', 'UPDATE','999',
                              'ald',   'KIND_ALD',                  'TYPE_ALD', 'UPDATE','999',
                              'ch3o2', 'KIND_CH3O2',                'TYPE_CH3O2', 'UPDATE','999',
                              'c3h8',  'KIND_C3H8',                 'TYPE_C3H8', 'UPDATE','999',
                              'c2h6',  'KIND_C2H6',                 'TYPE_C2H6', 'UPDATE','999',
                              'acet',  'KIND_ACET',                 'TYPE_ACET', 'UPDATE','999',
                              'hcho',  'KIND_HCHO',                 'TYPE_HCHO', 'UPDATE','999',
                              'c2h4',  'KIND_C2H4',                 'TYPE_C2H4', 'UPDATE','999',
                              'c3h6',  'KIND_C3H6',                 'TYPE_C3H6', 'UPDATE','999',
                              'tol',   'KIND_TOL',                  'TYPE_TOL', 'UPDATE','999',
                              'mvk',   'KIND_MVK',                  'TYPE_MVK', 'UPDATE','999',
                              'bigalk','KIND_BIGALK',               'TYPE_BIGALK', 'UPDATE','999',
                              'isopr', 'KIND_ISOPR',                'TYPE_ISOPR', 'UPDATE','999',
                              'macr',  'KIND_MACR',                 'TYPE_MACR', 'UPDATE','999',
                              'glyald','KIND_GLYALD',               'TYPE_GLYALD', 'UPDATE','999',
                              'c10h16','KIND_C10H16',               'TYPE_C10H16', 'UPDATE','999',"
   export NL_WRF_STATE_BOUNDS="'QVAPOR','0.0','NULL','CLAMP',
                           'QRAIN', '0.0','NULL','CLAMP',
                           'QCLOUD','0.0','NULL','CLAMP',
                           'QSNOW', '0.0','NULL','CLAMP',
                           'QICE',  '0.0','NULL','CLAMP',
                           'o3',    '0.0','NULL','CLAMP',
                           'co',    '0.0','NULL','CLAMP',
                           'no',    '0.0','NULL','CLAMP',
                           'no2',   '0.0','NULL','CLAMP',
                           'hno3',  '0.0','NULL','CLAMP',
                           'hno4',  '0.0','NULL','CLAMP',
                           'n2o5',  '0.0','NULL','CLAMP',
                           'pan',   '0.0','NULL','CLAMP',
                           'mek',   '0.0','NULL','CLAMP',
                           'ald',   '0.0','NULL','CLAMP',
                           'ch3o2', '0.0','NULL','CLAMP',
                           'c3h8',  '0.0','NULL','CLAMP',
                           'c2h6',  '0.0','NULL','CLAMP',
                           'acet'   '0.0','NULL','CLAMP',
                           'hcho'   '0.0','NULL','CLAMP',
                           'c2h4',  '0.0','NULL','CLAMP',
                           'c3h6',  '0.0','NULL','CLAMP',
                           'tol',   '0.0','NULL','CLAMP',
                           'mvk',   '0.0','NULL','CLAMP',
                           'bigalk','0.0','NULL','CLAMP',
                           'isopr', '0.0','NULL','CLAMP',
                           'macr',  '0.0','NULL','CLAMP',
                           'glyald','0.0','NULL','CLAMP',
                           'c10h16','0.0','NULL','CLAMP',"
   export NL_OUTPUT_STATE_VECTOR=.false.
   export NL_NUM_DOMAINS=${DOMAIN}
   export NL_CALENDAR_TYPE=3
   export NL_ASSIMILATION_PERIOD_SECONDS=${CYCLE_PERIOD_SEC}
   export NL_VERT_LOCALIZATION_COORD=3
   export NL_CENTER_SEARCH_HALF_LENGTH=500000.
   export NL_CENTER_SPLINE_GRID_SCALE=10
   export NL_SFC_ELEV_MAX_DIFF=100.0
   export NL_CIRCULATION_PRES_LEVEL=80000.0
   export NL_CIRCULATION_RADIUS=108000.0
#
# &restart_file_utility_nml
   export NL_SINGLE_RESTART_FILE_IN=.false.       
   export NL_SINGLE_RESTART_FILE_OUT=.false.       
#
# &dart_to_wrf_nml
   export NL_MODEL_ADVANCE_FILE=.false.
   export NL_ADV_MOD_COMMAND="'mpirun -np 64 ./wrf.exe'"
   export NL_DART_RESTART_NAME="'dart_wrf_vector'"
#
# &preprocess_nml
   export NL_INPUT_OBS_KIND_MOD_FILE="'${DART_DIR}/obs_kind/DEFAULT_obs_kind_mod.F90'"
   export NL_OUTPUT_OBS_KIND_MOD_FILE="'${DART_DIR}/obs_kind/obs_kind_mod.f90'"
   export NL_INPUT_OBS_DEF_MOD_FILE="'${DART_DIR}/obs_kind/DEFAULT_obs_def_mod.F90'"
   export NL_OUTPUT_OBS_DEF_MOD_FILE="'${DART_DIR}/obs_kind/obs_def_mod.f90'"
   export NL_INPUT_FILES="'${DART_DIR}/obs_def/obs_def_reanalysis_bufr_mod.f90',
                       '${DART_DIR}/obs_def/obs_def_radar_mod.f90',
                       '${DART_DIR}/obs_def/obs_def_metar_mod.f90',
                       '${DART_DIR}/obs_def/obs_def_dew_point_mod.f90',
                       '${DART_DIR}/obs_def/obs_def_altimeter_mod.f90',
                       '${DART_DIR}/obs_def/obs_def_gps_mod.f90',
                       '${DART_DIR}/obs_def/obs_def_gts_mod.f90',
                       '${DART_DIR}/obs_def/obs_def_vortex_mod.f90',
                       '${DART_DIR}/obs_def/obs_def_MOPITT_CO_mod.f90'"
#
# &obs_kind_nml
   export NL_ASSIMILATE_THESE_OBS_TYPES="'RADIOSONDE_TEMPERATURE',
                                      'RADIOSONDE_U_WIND_COMPONENT',
                                      'RADIOSONDE_V_WIND_COMPONENT',
                                      'ACARS_U_WIND_COMPONENT',
                                      'ACARS_V_WIND_COMPONENT',
                                      'ACARS_TEMPERATURE',
                                      'AIRCRAFT_U_WIND_COMPONENT',
                                      'AIRCRAFT_V_WIND_COMPONENT',
                                      'AIRCRAFT_TEMPERATURE',
                                      'SAT_U_WIND_COMPONENT',
                                      'SAT_V_WIND_COMPONENT',
                                      'MOPITT_CO_RETRIEVAL',
                                      'IASI_CO_RETRIEVAL',
                                      'IASI_O3_RETRIEVAL'"
#
# &replace_wrf_fields_nml
   export NL_FIELDNAMES="'SNOWC',
                      'ALBBCK',
                      'TMN',
                      'TSK',
                      'SH2O',
                      'SMOIS',
                      'SEAICE',
                      'HGT_d01',
                      'TSLB',
                      'SST',
                      'SNOWH',
                      'SNOW'"
   export NL_FIELDLIST_FILE="' '"
#
rm -rf input.nml
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_input.nml.ksh
#
###############################################################################
#
# RUN OBS_DIAG
#
###############################################################################
#
cd ${RUN_DIR}
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
./obs_diag
#
# Remove work/run files
rm -rf advance_time
rm -rf dart_log.nml
rm -rf dart_log.out
rm -rf file_list.txt
rm -rf input.nml
rm -rf LargeInnov.txt
rm -rf obs_diag
exit
