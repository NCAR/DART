#!/bin/ksh -aeux
#########################################################################
#
# Purpose: Set global environment variables for real_time_wrf_chem
#
#########################################################################
#
# CYCLE DATE-TIME:
#export DATE=2008072000
export DATE=2008072006
#export DATE=2008072012
#export DATE=2008072018
export INITIAL_DATE=2008072000
export FIRST_FILTER_DATE=2008072006
export CYCLE_PERIOD=6
#
export BUILD_DIR=/glade/p/work/mizzi/TRUNK/WRFDAv3.4_dmpar/var/da
export DART_DIR=/glade/p/work/mizzi/TRUNK/DART_CHEM_MY_BRANCH
cp ${DART_DIR}/models/wrf_chem/work/advance_time ./.
cp ${DART_DIR}/models/wrf_chem/work/input.nml ./.
export YYYY=$(echo $DATE | cut -c1-4)
export YY=$(echo $DATE | cut -c3-4)
export MM=$(echo $DATE | cut -c5-6)
export DD=$(echo $DATE | cut -c7-8)
export HH=$(echo $DATE | cut -c9-10)
export FILE_DATE=${YYYY}-${MM}-${DD}_${HH}:00:00
export PAST_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} -${CYCLE_PERIOD} 2>/dev/null)
export PAST_YYYY=$(echo $PAST_DATE | cut -c1-4)
export PAST_YY=$(echo $PAST_DATE | cut -c3-4)
export PAST_MM=$(echo $PAST_DATE | cut -c5-6)
export PAST_DD=$(echo $PAST_DATE | cut -c7-8)
export PAST_HH=$(echo $PAST_DATE | cut -c9-10)
export PAST_FILE_DATE=${PAST_YYYY}-${PAST_MM}-${PAST_DD}_${PAST_HH}:00:00
export NEXT_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} +${CYCLE_PERIOD} 2>/dev/null)
export NEXT_YYYY=$(echo $NEXT_DATE | cut -c1-4)
export NEXT_YY=$(echo $NEXT_DATE | cut -c3-4)
export NEXT_MM=$(echo $NEXT_DATE | cut -c5-6)
export NEXT_DD=$(echo $NEXT_DATE | cut -c7-8)
export NEXT_HH=$(echo $NEXT_DATE | cut -c9-10)
export NEXT_FILE_DATE=${NEXT_YYYY}-${NEXT_MM}-${NEXT_DD}_${NEXT_HH}:00:00
#
# DART TIME DATA
export DT_YYYY=${YYYY}
export DT_YY=${YY}
export DT_MM=${MM} 
export DT_DD=${DD} 
export DT_HH=${HH} 
(( DT_MM = ${DT_MM} + 0 ))
(( DT_DD = ${DT_DD} + 0 ))
(( DT_HH = ${DT_HH} + 0 ))
if [[ ${HH} -eq 0 ]]; then
   export TMP_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} -1 2>/dev/null)
   export TMP_YYYY=$(echo $TMP_DATE | cut -c1-4)
   export TMP_YY=$(echo $TMP_DATE | cut -c3-4)
   export TMP_MM=$(echo $TMP_DATE | cut -c5-6)
   export TMP_DD=$(echo $TMP_DATE | cut -c7-8)
   export TMP_HH=$(echo $TMP_DATE | cut -c9-10)
   export D_YYYY=${TMP_YYYY}
   export D_YY=${TMP_YY}
   export D_MM=${TMP_MM}
   export D_DD=${TMP_DD}
   export D_HH=24
   (( DD_MM = ${D_MM} + 0 ))
   (( DD_DD = ${D_DD} + 0 ))
   (( DD_HH = ${D_HH} + 0 ))
else
   export D_YYYY=${YYYY}
   export D_YY=${YY}
   export D_MM=${MM}
   export D_DD=${DD}
   export D_HH=${HH}
   (( DD_MM = ${D_MM} + 0 ))
   (( DD_DD = ${D_DD} + 0 ))
   (( DD_HH = ${D_HH} + 0 ))
fi
export D_DATE=${D_YYYY}${D_MM}${D_DD}${D_HH}
#
# CALCULATE GREGORIAN TIMES FOR START AND END OF ASSIMILATION WINDOW
set -A GREG_DATA `echo $DATE 0 -g | ${DART_DIR}/models/wrf_chem/work/advance_time`
export DAY_GREG=${GREG_DATA[0]}
export SEC_GREG=${GREG_DATA[1]}
set -A GREG_DATA `echo $NEXT_DATE 0 -g | ${DART_DIR}/models/wrf_chem/work/advance_time`
export NEXT_DAY_GREG=${GREG_DATA[0]}
export NEXT_SEC_GREG=${GREG_DATA[1]}
export ASIM_WINDOW=3
export ASIM_MIN_DATE=$($BUILD_DIR/da_advance_time.exe $DATE -$ASIM_WINDOW 2>/dev/null)
export ASIM_MAX_DATE=$($BUILD_DIR/da_advance_time.exe $DATE +$ASIM_WINDOW 2>/dev/null)
set -A temp `echo $ASIM_MIN_DATE 0 -g | ${DART_DIR}/models/wrf_chem/work/advance_time`
export ASIM_MIN_DAY_GREG=${temp[0]}
export ASIM_MIN_SEC_GREG=${temp[1]}
set -A temp `echo $ASIM_MAX_DATE 0 -g | ${DART_DIR}/models/wrf_chem/work/advance_time` 
export ASIM_MAX_DAY_GREG=${temp[0]}
export ASIM_MAX_SEC_GREG=${temp[1]}
#
# SELECT COMPOENT RUN OPTIONS:
export RUN_GEOGRID=false
export RUN_UNGRIB=true
export RUN_METGRID=true
export RUN_REAL=true
export RUN_PERT_WRFCHEM_MET_IC=true
export RUN_PERT_WRFCHEM_MET_BC=true
export RUN_EXO_COLDENS=true
export RUN_SEASON_WES=true
export RUN_WRFCHEM_BIO=true
export RUN_WRFCHEM_FIRE=true
export RUN_WRFCHEM_CHEMI=true
export RUN_PERT_WRFCHEM_CHEM_ICBC=true
export RUN_PERT_WRFCHEM_CHEM_EMISS=true
export RUN_MOPITT_CO_OBS=true
export RUN_IASI_CO_OBS=false
export RUN_IASI_O3_OBS=false
export RUN_MET_OBS=true
export RUN_COMBINE_OBS=true
export RUN_PREPROCESS_OBS=true
#
if [[ ${DATE} -eq ${INITIAL_DATE}  ]]; then
   export RUN_WRFCHEM_INITIAL=true
   export RUN_DART_FILTER=false
   export RUN_UPDATE_BC=false
   export RUN_WRFCHEM_CYCLE_CR=false
   export RUN_WRFCHEM_CYCLE_FR=false
else
   export RUN_WRFCHEM_INITIAL=false
   export RUN_DART_FILTER=true
   export RUN_UPDATE_BC=true
   export RUN_WRFCHEM_CYCLE_CR=true
   export RUN_WRFCHEM_CYCLE_FR=true
fi
#
export NL_APM_SCALE=1.
export NL_APM_SCALE_SW=.FALSE.
#
# FORECAST PARAMETERS:
export USE_DART_INFL=true
export FCST_PERIOD=24
(( CYCLE_PERIOD_SEC=${CYCLE_PERIOD}*60*60 ))
export NUM_MEMBERS=20
export MAX_DOMAINS=02
export CR_DOMAIN=01
export FR_DOMAIN=02
export NNXP_CR=180
export NNYP_CR=140
export NNZP_CR=37
export NNXP_FR=321
export NNYP_FR=291
export NNZP_FR=37
export ISTR_CR=1
export JSTR_CR=1
export ISTR_FR=86
export JSTR_FR=35
export DX_CR=15000
export DX_FR=3000
(( LBC_END=2*${FCST_PERIOD} ))
export LBC_FREQ=3
(( INTERVAL_SECONDS=${LBC_FREQ}*60*60 ))
export LBC_START=0
export START_DATE=${DATE}
export END_DATE=$($BUILD_DIR/da_advance_time.exe ${START_DATE} ${FCST_PERIOD} 2>/dev/null)
export START_YEAR=$(echo $START_DATE | cut -c1-4)
export START_MONTH=$(echo $START_DATE | cut -c5-6)
export START_DAY=$(echo $START_DATE | cut -c7-8)
export START_HOUR=$(echo $START_DATE | cut -c9-10)
export START_FILE_DATE=${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:00:00
export END_YEAR=$(echo $END_DATE | cut -c1-4)
export END_MONTH=$(echo $END_DATE | cut -c5-6)
export END_DAY=$(echo $END_DATE | cut -c7-8)
export END_HOUR=$(echo $END_DATE | cut -c9-10)
export END_FILE_DATE=${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:00:00
#
# LARGE SCALE FORECAST PARAMETERS:
export FG_TYPE=GFS
export GRIB_PART1=gfs_4_
export GRIB_PART2=.g2.tar
#
# COMPUTER PARAMETERS:
export PROJ_NUMBER=P19010000
export PROJ_NUMBER=NACD0002
export PROJ_NUMBER=NACD0002
export GEOGRID_TIME_LIMIT=0:10
export GEOGRID_NUM_TASKS=32
export GEOGRID_TASKS_PER_NODE=8
export GEOGRID_JOB_CLASS=regular
export WRFCHEM_TIME_LIMIT=6:00
export WRFCHEM_NUM_TASKS=256
export WRFCHEM_TASKS_PER_NODE=16
export WRFCHEM_JOB_CLASS=regular
export WRFDA_TIME_LIMIT=0:10
export WRFDA_NUM_TASKS=32
export WRFDA_TASKS_PER_NODE=8
export WRFDA_JOB_CLASS=premium
export FILTER_TIME_LIMIT=0:59
export FILTER_NUM_TASKS=32
export FILTER_TASKS_PER_NODE=8
export FILTER_JOB_CLASS=regular
#
# CODE VERSIONS:
export WPS_VER=WPSv3.6.1_dmpar
export WPS_GEOG_VER=WPSv3.6.1_GEOG_DATA
export WRFDA_VER=WRFDAv3.4_dmpar
export WRFDA_TOOLS_VER=WRFDA_TOOLSv3.4
export WRF_VER=WRFv3.6.1_dmpar
export WRFCHEM_VER=WRFCHEMv3.6.1_dmpar
export DART_VER=DART_CHEM_MY_BRANCH
#
# ROOT DIRECTORIES:
export SCRATCH_DIR=/glade/scratch/mizzi
export WORK_DIR=/glade/p/work/mizzi
export ACD_DIR=/glade/p/acd/mizzi
#
# DEPENDENT DIRECTORIES:
export RUN_DIR=${SCRATCH_DIR}/real_my_branch_FRAPPE_v3.6.1
export TRUNK_DIR=${WORK_DIR}/TRUNK
export WPS_DIR=${TRUNK_DIR}/${WPS_VER}
export WPS_GEOG_DIR=${TRUNK_DIR}/${WPS_GEOG_VER}/geog
export WRFCHEM_DIR=${TRUNK_DIR}/${WRFCHEM_VER}
export WRFVAR_DIR=${TRUNK_DIR}/${WRFDA_VER}
export DART_DIR=${TRUNK_DIR}/${DART_VER}
export BUILD_DIR=${WRFVAR_DIR}/var/da
export WRF_DIR=${TRUNK_DIR}/${WRF_VER}
export HYBRID_TRUNK_DIR=${WORK_DIR}/HYBRID_TRUNK
export HYBRID_SCRIPTS_DIR=${HYBRID_TRUNK_DIR}/hybrid_scripts
export SCRIPTS_DIR=${TRUNK_DIR}/${WRFDA_TOOLS_VER}/scripts
export REAL_TIME_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_REAL_TIME
#
export INPUT_DATA_DIR=${ACD_DIR}/AVE_TEST_DATA
export GRIB_DIR=${INPUT_DATA_DIR}/GFS_DATA_NOMADS
export OBSPROC_DIR=${WRFVAR_DIR}/var/obsproc
export COLDENS_INPUT_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_INITIAL_EMISS/WES-COLDENS
export VTABLE_DIR=${WPS_DIR}/ungrib/Variable_Tables
export OB_PREPBUFR_DIR=${INPUT_DATA_DIR}/obs_MET
export BE_DIR=${WRFVAR_DIR}/var/run
export PREPBUFR_DATA_DIR=${INPUT_DATA_DIR}/obs_MET
export MEGAN_BIO_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_INITIAL_EMISS/MEGAN-BIO
export MEGAN_INPUT_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_INITIAL_EMISS/MEGAN-DATA
export FIRE_INPUT_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_INITIAL_EMISS/WRF-FIRE
export PERT_CHEM_INPUT_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_PERT_CHEM/ICBC_PERT
export PERT_CHEM_EMISS_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_PERT_CHEM/EMISS_PERT
export GEOGRID_DIR=${RUN_DIR}/geogrid
export METGRID_DIR=${RUN_DIR}/${DATE}/metgrid
export REAL_DIR=${RUN_DIR}/${DATE}/real
export WRFCHEM_MET_IC_DIR=${RUN_DIR}/${DATE}/wrfchem_met_ic
export WRFCHEM_MET_BC_DIR=${RUN_DIR}/${DATE}/wrfchem_met_bc
export EXO_COLDENS_DIR=${RUN_DIR}/${DATE}/exo_coldens
export SEASONS_WES_DIR=${RUN_DIR}/${DATE}/seasons_wes
export WRFCHEM_BIO_DIR=${RUN_DIR}/${DATE}/wrfchem_bio
export WRFCHEM_FIRE_DIR=${RUN_DIR}/${DATE}/wrfchem_fire
export WRFCHEM_CHEMI_DIR=${RUN_DIR}/${DATE}/wrfchem_chemi
export WRFCHEM_CHEM_EMISS_DIR=${RUN_DIR}/${DATE}/wrfchem_chem_emiss
export WRFCHEM_INITIAL_DIR=${RUN_DIR}/${INITIAL_DATE}/wrfchem_initial
export WRFCHEM_CYCLE_CR_DIR=${RUN_DIR}/${DATE}/wrfchem_cycle_cr
export WRFCHEM_CYCLE_FR_DIR=${RUN_DIR}/${DATE}/wrfchem_cycle_fr
export WRFCHEM_LAST_CYCLE_CR_DIR=${RUN_DIR}/${PAST_DATE}/wrfchem_cycle_cr
export PREPBUFR_MET_OBS_DIR=${RUN_DIR}/${DATE}/prepbufr_met_obs
export MOPITT_CO_OBS_DIR=${RUN_DIR}/${DATE}/mopitt_co_obs
export IASI_CO_OBS_DIR=${RUN_DIR}/${DATE}/iasi_co_obs
export IASI_O3_OBS_DIR=${RUN_DIR}/${DATE}/iasi_o3_obs
export COMBINE_OBS_DIR=${RUN_DIR}/${DATE}/combine_obs
export PREPROCESS_OBS_DIR=${RUN_DIR}/${DATE}/preprocess_obs
export WRFCHEM_CHEM_ICBC_DIR=${RUN_DIR}/${DATE}/wrfchem_chem_icbc
export DART_FILTER_DIR=${RUN_DIR}/${DATE}/dart_filter
export UPDATE_BC_DIR=${RUN_DIR}/${DATE}/update_bc
#
# WPS PARAMETERS:
export SINGLE_FILE=false
export HOR_SCALE=1500
export VTABLE_TYPE=GFS
export METGRID_TABLE_TYPE=ARW
#
# WRF PREPROCESS PARAMETERS
export NL_MIN_LAT=7.
export NL_MAX_LAT=54.
export NL_MIN_LON=184.
export NL_MAX_LON=310.
#
# PERT CHEM PARAMETERS
export MOZ_SPREAD=0.30
export NL_MEAN=1.0
export NL_SPREAD=0.30
#
#########################################################################
#
#  NAMELIST PARAMETERS
#
#########################################################################
#
# WPS SHARE NAMELIST:
export NL_WRF_CORE="'"ARW"'"
export NL_MAX_DOM=${MAX_DOMAINS}
export NL_IO_FORM_GEOGRID=2
export NL_DEBUG_LEVEL=0
export NL_OPT_OUTPUT_FROM_GEOGRID_PATH="'"${GEOGRID_DIR}"'"
export NL_ACTIVE_GRID=".true.",".true."
#
# WPS GEOGRID NAMELIST:
export NL_S_WE=1,1
export NL_E_WE=${NNXP_CR},${NNXP_FR}
export NL_S_SN=1,1
export NL_E_SN=${NNYP_CR},${NNYP_FR}
export NL_S_VERT=1,1
export NL_E_VERT=${NNZP_CR},${NNZP_FR}
export NL_PARENT_ID="1,1"
export NL_PARENT_GRID_RATIO=1,5
export NL_I_PARENT_START=${ISTR_CR},${ISTR_FR}
export NL_J_PARENT_START=${JSTR_CR},${JSTR_FR}
export NL_GEOG_DATA_RES="'"30s"'","'"30s"'"
export NL_DX=${DX_CR}
export NL_DY=${DX_CR}
export NL_MAP_PROJ="'"lambert"'"
export NL_REF_LAT=40.0
export NL_REF_LON=-112.0
export NL_STAND_LON=-105.0
export NL_TRUELAT1=30.0
export NL_TRUELAT2=60.0
export NL_GEOG_DATA_PATH="'"${WPS_GEOG_DIR}"'"
export NL_OPT_GEOGRID_TBL_PATH="'"${WPS_DIR}/geogrid"'"
#
# WPS UNGRIB NAMELIST:
export NL_OUT_FORMAT="'"WPS"'"
#
# WPS METGRID NAMELIST:
export NL_IO_FORM_METGRID=2
#
# WRF NAMELIST:
# TIME CONTROL NAMELIST:
export NL_RUN_DAYS=0
export NL_RUN_HOURS=${FCST_PERIOD}
export NL_RUN_MINUTES=0
export NL_RUN_SECONDS=0
export NL_START_YEAR=${START_YEAR},${START_YEAR}
export NL_START_MONTH=${START_MONTH},${START_MONTH}
export NL_START_DAY=${START_DAY},${START_DAY}
export NL_START_HOUR=${START_HOUR},${START_HOUR}
export NL_START_MINUTE=00,00
export NL_START_SECOND=00,00
export NL_END_YEAR=${END_YEAR},${END_YEAR}
export NL_END_MONTH=${END_MONTH},${END_MONTH}
export NL_END_DAY=${END_DAY},${END_DAY}
export NL_END_HOUR=${END_HOUR},${END_HOUR}
export NL_END_MINUTE=00,00
export NL_END_SECOND=00,00
export NL_INTERVAL_SECONDS=${INTERVAL_SECONDS}
export NL_INPUT_FROM_FILE=".true.",".true."
export NL_HISTORY_INTERVAL=360,360
export NL_FRAMES_PER_OUTFILE=1,1
export NL_RESTART=".false."
export NL_RESTART_INTERVAL=1440
export NL_IO_FORM_HISTORY=2
export NL_IO_FORM_RESTART=2
export NL_FINE_INPUT_STREAM=0,2
export NL_IO_FORM_INPUT=2
export NL_IO_FORM_BOUNDARY=2
export NL_AUXINPUT2_INNAME="'"wrfinput_d\<domain\>"'"
export NL_AUXINPUT5_INNAME="'"wrfchemi_d\<domain\>_\<date\>"'"
export NL_AUXINPUT6_INNAME="'"wrfbiochemi_d\<domain\>_\<date\>"'"
export NL_AUXINPUT7_INNAME="'"wrffirechemi_d\<domain\>_\<date\>"'"
export NL_AUXINPUT2_INTERVAL_M=60480,60480
export NL_AUXINPUT5_INTERVAL_M=60,60
export NL_AUXINPUT6_INTERVAL_M=60480,60480
export NL_AUXINPUT7_INTERVAL_M=60,60
export NL_FRAMES_PER_AUXINPUT2=1,1
export NL_FRAMES_PER_AUXINPUT5=1,1
export NL_FRAMES_PER_AUXINPUT6=1,1
export NL_FRAMES_PER_AUXINPUT7=1,1
export NL_IO_FORM_AUXINPUT2=2
export NL_IO_FORM_AUXINPUT5=2
export NL_IO_FORM_AUXINPUT6=2
export NL_IO_FORM_AUXINPUT7=2
export NL_IOFIELDS_FILENAME="'"hist_io_flds"'","'"hist_io_flds"'"
export NL_WRITE_INPUT=".true."
export NL_INPUTOUT_INTERVAL=360
export NL_INPUT_OUTNAME="'"wrfapm_d\<domain\>_\<date\>"'"
export NL_DEBUG_LEVEL=0
#
# DOMAINS NAMELIST:
export NL_TIME_STEP=90
export NL_TIME_STEP_FRACT_NUM=0
export NL_TIME_STEP_FRACT_DEN=1
export NL_MAX_DOM=${MAX_DOMAINS}
export NL_S_WE=1,1
export NL_E_WE=${NNXP_CR},${NNXP_FR}
export NL_S_SN=1,1
export NL_E_SN=${NNYP_CR},${NNYP_FR}
export NL_S_VERT=1,1
export NL_E_VERT=${NNZP_CR},${NNZP_FR}
export NL_NUM_METGRID_LEVELS=27
export NL_NUM_METGRID_SOIL_LEVELS=4
export NL_DX=${DX_CR},${DX_FR}
export NL_DY=${DX_CR},${DX_FR}
export NL_GRID_ID=1,2
export NL_PARENT_ID=0,1
export NL_I_PARENT_START=${ISTR_CR},${ISTR_FR}
export NL_J_PARENT_START=${JSTR_CR},${JSTR_FR}
export NL_PARENT_GRID_RATIO=1,5
export NL_PARENT_TIME_STEP_RATIO=1,5
export NL_FEEDBACK=0
export NL_SMOOTH_OPTION=1
export NL_LAGRANGE_ORDER=2
export NL_INTERP_TYPE=2
export NL_EXTRAP_TYPE=2
export NL_T_EXTRAP_TYPE=2
export NL_USE_SURFACE=".true."
export NL_USE_LEVELS_BELOW_GROUND=".true."
export NL_LOWEST_LEV_FROM_SFC=".false."
export NL_FORCE_SFC_IN_VINTERP=1
export NL_ZAP_CLOSE_LEVELS=500
export NL_INTERP_THETA=".false."
export NL_HYPSOMETRIC_OPT=2
export NL_P_TOP_REQUESTED=5000.
export NL_ETA_LEVELS=1.000000,0.996200,0.989737,0.982460,0.974381,0.965422,\
0.955498,0.944507,0.932347,0.918907,0.904075,0.887721,0.869715,0.849928,\
0.828211,0.804436,0.778472,0.750192,0.719474,0.686214,0.650339,0.611803,\
0.570656,0.526958,0.480854,0.432582,0.382474,0.330973,0.278674,0.226390,\
0.175086,0.132183,0.096211,0.065616,0.039773,0.018113,0.000000,
#
# PHYSICS NAMELIST:
export NL_MP_PHYSICS=8,8
export NL_RA_LW_PHYSICS=4,4
export NL_RA_SW_PHYSICS=4,4
export NL_RADT=15,3
export NL_SF_SFCLAY_PHYSICS=1,1
export NL_SF_SURFACE_PHYSICS=2,2
export NL_BL_PBL_PHYSICS=1,1
export NL_BLDT=0,0
export NL_CU_PHYSICS=5,5
export NL_CUDT=1,1
export NL_CUGD_AVEDX=1
export NL_CU_RAD_FEEDBACK=".true.",".true."
export NL_CU_DIAG=1,1
export NL_ISFFLX=1
export NL_IFSNOW=0
export NL_ICLOUD=1
export NL_SURFACE_INPUT_SOURCE=1
export NL_NUM_SOIL_LAYERS=4
export NL_MP_ZERO_OUT=2
export NL_NUM_LAND_CAT=24
export NL_SF_URBAN_PHYSICS=0,1
export NL_MAXIENS=1
export NL_MAXENS=3
export NL_MAXENS2=3
export NL_MAXENS3=16
export NL_ENSDIM=144
#
# DYNAMICS NAMELIST:
export NL_ISO_TEMP=200.
export NL_TRACER_OPT=0,0
export NL_W_DAMPING=1
export NL_DIFF_OPT=1
export NL_DIFF_6TH_OPT=0,0
export NL_DIFF_6TH_FACTOR=0.12,0.12
export NL_KM_OPT=4
export NL_DAMP_OPT=1
export NL_ZDAMP=5000,5000
export NL_DAMPCOEF=0.15,0.15
export NL_NON_HYDROSTATIC=".true.",".true."
export NL_USE_BASEPARAM_FR_NML=".true."
export NL_MOIST_ADV_OPT=2,2
export NL_SCALAR_ADV_OPT=2,2
export NL_CHEM_ADV_OPT=2,2
export NL_TKE_ADV_OPT=2,2
export NL_H_MOM_ADV_ORDER=5,5
export NL_V_MOM_ADV_ORDER=3,3
export NL_H_SCA_ADV_ORDER=5,5
export NL_V_SCA_ADV_ORDER=3,3
#
# BDY_CONTROL NAMELIST:
export NL_SPEC_BDY_WIDTH=5
export NL_SPEC_ZONE=1
export NL_RELAX_ZONE=4
export NL_SPECIFIED=".true.",".false."
export NL_NESTED=".false.",".true."
#
# QUILT NAMELIST:
export NL_NIO_TASKS_PER_GROUP=0
export NL_NIO_GROUPS=1
#
# NAMELIST CHEM
export NL_KEMIT=11
export NL_CHEM_OPT=112,112
export NL_BIOEMDT=3,1
export NL_PHOTDT=3,1
export NL_CHEMDT=3,1
export NL_IO_STYLE_EMISSIONS=2
export NL_EMISS_INPT_OPT=111,111
export NL_EMISS_OPT=8,8
export NL_EMISS_OPT_VOL=0,0
export NL_CHEM_IN_OPT=0,0
export NL_PHOT_OPT=3,3
export NL_GAS_DRYDEP_OPT=1,1
export NL_AER_DRYDEP_OPT=1,1
export NL_BIO_EMISS_OPT=3,3
export NL_NE_AREA=118
export NL_GAS_BC_OPT=112,112
export NL_GAS_IC_OPT=112,112
export NL_GAS_BC_OPT=112,112
export NL_AER_BC_OPT=112,112
export NL_AER_IC_OPT=112,112
export NL_GASCHEM_ONOFF=1,1
export NL_AERCHEM_ONOFF=1,1
export NL_WETSCAV_ONOFF=1,1
export NL_CLDCHEM_ONOFF=0,0
export NL_VERTMIX_ONOFF=1,1
export NL_CHEM_CONV_TR=1,1
export NL_CONV_TR_WETSCAV=1,1
export NL_CONV_TR_AQCHEM=0,0
export NL_SEAS_OPT=0
export NL_DUST_OPT=0
export NL_DMSEMIS_OPT=1
export NL_BIOMASS_BURN_OPT=2,2
export NL_PLUMERISEFIRE_FRQ=15,15
export NL_SCALE_FIRE_EMISS=".true.",".true."
export NL_HAVE_BCS_CHEM=".true.",".true."
export NL_AER_RA_FEEDBACK=1,1
export NL_AER_OP_OPT=1
export NL_OPT_PARS_OUT=1
export NL_HAVE_BCS_UPPER=".false.",".false."
export NL_FIXED_UBC_PRESS=50.,50.
export NL_FIXED_UBC_INNAME="'"ubvals_b40.20th.track1_1996-2005.nc"'"
#
# WRFDA NAMELIST PARAMETERS
# WRFVAR1 NAMELIST:
export NL_PRINT_DETAIL_GRAD=false
export NL_VAR4D=false
export NL_MULTI_INC=0
#
# WRFVAR3 NAMELIST:
export NL_OB_FORMAT=1
export NL_NUM_FGAT_TIME=1
#
# WRFVAR4 NAMELIST:
export NL_USE_SYNOPOBS=true
export NL_USE_SHIPOBS=false
export NL_USE_METAROBS=true
export NL_USE_SOUNDOBS=true
export NL_USE_MTGIRSOBS=false
export NL_USE_PILOTOBS=true
export NL_USE_AIREOBS=true
export NL_USE_GEOAMVOBS=false
export NL_USE_POLARAMVOBS=false
export NL_USE_BOGUSOBS=false
export NL_USE_BUOYOBS=false
export NL_USE_PROFILEROBS=false
export NL_USE_SATEMOBS=false
export NL_USE_GPSPWOBS=false
export NL_USE_GPSREFOBS=false
export NL_USE_SSMIRETRIEVALOBS=false
export NL_USE_QSCATOBS=false
export NL_USE_AIRSRETOBS=false
#
# WRFVAR5 NAMELIST:
export NL_CHECK_MAX_IV=true
export NL_PUT_RAND_SEED=true
#
# WRFVAR6 NAMELIST:
export NL_NTMAX=100
#
# WRFVAR7 NAMELIST:
export NL_VAR_SCALING4=1.0
export NL_JE_FACTOR=1.0
export NL_CV_OPTIONS=3
#
# WRFVAR11 NAMELIST:
export NL_CV_OPTIONS_HUM=1
export NL_CHECK_RH=2
export NL_SEED_ARRAY1=$(${BUILD_DIR}/da_advance_time.exe ${DATE} 0 -f hhddmmyycc)
export NL_SEED_ARRAY2=`echo ${NUM_MEMBERS} \* 100000 | bc -l `
export NL_CALCULATE_CG_COST_FN=true
export NL_LAT_STATS_OPTION=false
#
# WRFVAR15 NAMELIST:
export NL_NUM_PSEUDO=0
export NL_PSEUDO_X=0
export NL_PSEUDO_Y=0
export NL_PSEUDO_Z=0
export NL_PSEUDO_ERR=0.0
export NL_PSEUDO_VAL=0.0
#
# WRFVAR16 NAMELIST:
export NL_ALPHACV_METHOD=2
export NL_ENSDIM_ALPHA=0
export NL_ALPHA_CORR_TYPE=3
export NL_ALPHA_CORR_SCALE=${HOR_SCALE}
export NL_ALPHA_STD_DEV=1.0
export NL_ALPHA_VERTLOC=false
export NL_ALPHA_TRUNCATION=1
#
# WRFVAR17 NAMELIST:
export NL_ANALYSIS_TYPE="'"RANDOMCV"'"
#
# WRFVAR18 NAMELIST:
export ANALYSIS_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} 0 -W 2>/dev/null)
#
# WRFVAR19 NAMELIST:
export NL_PSEUDO_VAR="'"t"'"
#
# WRFVAR21 NAMELIST:
export NL_TIME_WINDOW_MIN="'"$(${BUILD_DIR}/da_advance_time.exe ${DATE} -${ASIM_WINDOW} -W 2>/dev/null)"'"
#
# WRFVAR22 NAMELIST:
export NL_TIME_WINDOW_MAX="'"$(${BUILD_DIR}/da_advance_time.exe ${DATE} +${ASIM_WINDOW} -W 2>/dev/null)"'"
#
# WRFVAR23 NAMELIST:
export NL_JCDFI_USE=false
export NL_JCDFI_IO=false
#
# DART input.nml parameters
# &apm namelist parameters
export NL_APM_SCALE=1.
export NL_APM_SCALE_SW=.FALSE.
#
# &filter.nml
export NL_OUTLIER_THRESHOLD=3.
export NL_ENABLE_SPECIAL_OUTLIER_CODE=.false.
export NL_SPECIAL_OUTLIER_THRESHOLD=4.
export NL_ENS_SIZE=${NUM_MEMBERS}
export NL_OUTPUT_RESTART=.true.
export NL_START_FROM_RESTART=.true.
export NL_OBS_SEQUENCE_IN_NAME="'obs_seq.out'"       
export NL_OBS_SEQUENCE_OUT_NAME="'obs_seq.final'"
export NL_RESTART_IN_FILE_NAME="'filter_ic_old'"       
export NL_RESTART_OUT_FILE_NAME="'filter_ic_new'"       
set -A temp `echo ${ASIM_MIN_DATE} 0 -g | ${DART_DIR}/models/wrf_chem/work/advance_time`
(( temp[1]=${temp[1]}+1 ))
export NL_FIRST_OBS_DAYS=${temp[0]}
export NL_FIRST_OBS_SECONDS=${temp[1]}
set -A temp `echo ${ASIM_MAX_DATE} 0 -g | ${DART_DIR}/models/wrf_chem/work/advance_time`
export NL_LAST_OBS_DAYS=${temp[0]}
export NL_LAST_OBS_SECONDS=${temp[1]}
export NL_NUM_OUTPUT_STATE_MEMBERS=0
export NL_NUM_OUTPUT_OBS_MEMBERS=${NUM_MEMBERS}
if ${USE_DART_INFL}; then
   export NL_INF_FLAVOR_PRIOR=2
else 
   export NL_INF_FLAVOR_PRIOR=0
fi
export NL_INF_FLAVOR_POST=0  
if [[ ${START_DATE} -eq ${FIRST_FILTER_DATE} ]]; then
   export NL_INF_INITIAL_FROM_RESTART_PRIOR=.false.
   export NL_INF_SD_INITIAL_FROM_RESTART_PRIOR=.false.
   export NL_INF_INITIAL_FROM_RESTART_POST=.false.
   export NL_INF_SD_INITIAL_FROM_RESTART_POST=.false.
else
   export NL_INF_INITIAL_FROM_RESTART_PRIOR=.true.
   export NL_INF_SD_INITIAL_FROM_RESTART_PRIOR=.true.
   export NL_INF_INITIAL_FROM_RESTART_POST=.true.
   export NL_INF_SD_INITIAL_FROM_RESTART_POST=.true.
fi
export NL_INF_IN_FILE_NAME_PRIOR="'prior_inflate_ic_old'"
export NL_INF_IN_FILE_NAME_POST="'post_inflate_ics'"
export NL_INF_OUT_FILE_NAME_PRIOR="'prior_inflate_ic_new'"
export NL_INF_OUT_FILE_NAME_POST="'prior_inflate_restart'"
export NL_INF_DIAG_FILE_NAME_PRIOR="'prior_inflate_diag'"
export NL_INF_DIAG_FILE_NAME_POST="'post_inflate_diag'"
export NL_INF_INITIAL_PRIOR=1.0
export NL_INF_INITIAL_POST=1.0
export NL_INF_SD_INITIAL_PRIOR=0.6
export NL_INF_SD_INITIAL_POST=0.0
export NL_INF_DAMPING_PRIOR=0.9
export NL_INF_DAMPING_POST=1.0
export NL_INF_LOWER_BOUND_PRIOR=1.0
export NL_INF_LOWER_BOUND_POST=1.0
export NL_INF_UPPER_BOUND_PRIOR=100.0
export NL_INF_UPPER_BOUND_POST=100.0
export NL_INF_SD_LOWER_BOUND_PRIOR=0.6
export NL_INF_SD_LOWER_BOUND_POST=0.0
#
# &assim_tools_nml
export NL_CUTOFF=0.1
export NL_SPECIAL_LOCALIZATION_OBS_TYPES="'IASI_CO_RETRIEVAL','MOPITT_CO_RETRIEVAL'"
export NL_SAMPLING_ERROR_CORRECTION=.true.
export NL_SPECIAL_LOCALIZATION_CUTOFFS=0.1,0.1
export NL_ADAPTIVE_LOCALIZATION_THRESHOLD=2000
#
# &ensemble_manager_nml
export NL_SINGLE_RESTART_FILE_IN=.false.       
export NL_SINGLE_RESTART_FILE_OUT=.false.       
#
# &assim_model_nml
export NL_WRITE_BINARY_RESTART_FILE=.true.
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
                           'c10h16','KIND_C10H16',               'TYPE_C10H16', 'UPDATE','999'"
export NL_WRF_STATE_BOUNDS="'QVAPOR','0.0','NULL','CLAMP',
                        'QRAIN', '0.0','NULL','CLAMP',
                        'QCLOUD','0.0','NULL','CLAMP',
                        'QSNOW', '0.0','NULL','CLAMP',
                        'QICE',  '0.0','NULL','CLAMP',
                        'o3',    '0.0','NULL','CLAMP',
                        'co',    '1.e-4','NULL','CLAMP',
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
                        'c10h16','0.0','NULL','CLAMP'"
export NL_OUTPUT_STATE_VECTOR=.false.
export NL_NUM_DOMAINS=${CR_DOMAIN}
export NL_CALENDAR_TYPE=3
export NL_ASSIMILATION_PERIOD_SECONDS=${CYCLE_PERIOD_SEC}
export NL_VERT_LOCALIZATION_COORD=3
export NL_CENTER_SEARCH_HALF_LENGTH=500000.
export NL_CENTER_SPLINE_GRID_SCALE=10
export NL_SFC_ELEV_MAX_DIFF=100.0
export NL_CIRCULATION_PRES_LEVEL=80000.0
export NL_CIRCULATION_RADIUS=108000.0
export NL_ALLOW_OBS_BELOW_VOL=.false.
#
# &obs_diag_nml
export NL_FIRST_BIN_CENTER_YY=${DT_YYYY}
export NL_FIRST_BIN_CENTER_MM=${DT_MM}
export NL_FIRST_BIN_CENTER_DD=${DT_DD}
export NL_FIRST_BIN_CENTER_HH=${DT_HH}
export NL_LAST_BIN_CENTER_YY=${DT_YYYY}
export NL_LAST_BIN_CENTER_MM=${DT_MM}
export NL_LAST_BIN_CENTER_DD=${DT_DD}
export NL_LAST_BIN_CENTER_HH=${DT_HH}
export NL_BIN_SEPERATION_YY=0
export NL_BIN_SEPERATION_MM=0
export NL_BIN_SEPERATION_DD=0
export NL_BIN_SEPERATION_HH=0
export NL_BIN_WIDTH_YY=0
export NL_BIN_WIDTH_MM=0
export NL_BIN_WIDTH_DD=0
export NL_BIN_WIDTH_HH=0
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
# &restart_file_tool_nml
export NL_INPUT_FILE_NAME="'assim_model_state_tp'"
export NL_OUTPUT_FILE_NAME="'assim_model_state_ic'"
export NL_OUTPUT_IS_MODEL_ADVANCE_FILE=.true.
export NL_OVERWRITE_ADVANCE_TIME=.true.
export NL_NEW_ADVANCE_DAYS=${NEXT_DAY_GREG}
export NL_NEW_ADVANCE_SECS=${NEXT_SEC_GREG}
#
# &preprocess_nml
export NL_INPUT_OBS_KIND_MOD_FILE="'"${DART_DIR}/obs_kind/DEFAULT_obs_kind_mod.F90"'"
export NL_OUTPUT_OBS_KIND_MOD_FILE="'"${DART_DIR}/obs_kind/obs_kind_mod.f90"'"
export NL_INPUT_OBS_DEF_MOD_FILE="'"${DART_DIR}/obs_kind/DEFAULT_obs_def_mod.F90"'"
export NL_OUTPUT_OBS_DEF_MOD_FILE="'"${DART_DIR}/obs_kind/obs_def_mod.f90"'"
export NL_INPUT_FILES="'${DART_DIR}/obs_def/obs_def_reanalysis_bufr_mod.f90',
                    '${DART_DIR}/obs_def/obs_def_radar_mod.f90',
                    '${DART_DIR}/obs_def/obs_def_metar_mod.f90',
                    '${DART_DIR}/obs_def/obs_def_dew_point_mod.f90',
                    '${DART_DIR}/obs_def/obs_def_altimeter_mod.f90',
                    '${DART_DIR}/obs_def/obs_def_gps_mod.f90',
                    '${DART_DIR}/obs_def/obs_def_gts_mod.f90',
                    '${DART_DIR}/obs_def/obs_def_vortex_mod.f90',
                    '${DART_DIR}/obs_def/obs_def_IASI_CO_mod.f90',
                    '${DART_DIR}/obs_def/obs_def_MOPITT_CO_mod.f90',
                    '${DART_DIR}/obs_def/obs_def_MODIS_AOD_mod.f90'"
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
                                   'MOPITT_CO_RETRIEVAL'"
#export NL_EVALUATE_THESE_OBS_TYPES="'MOPITT_CO_RETRIEVAL'"
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
# &location_nml
export NL_HORIZ_DIST_ONLY=.true.
export NL_VERTICAL_NORMALIZATION_HEIGHT=8000.0
#
# ASSIMILATION WINDOW PARAMETERS
export ASIM_DATE_MIN=$(${BUILD_DIR}/da_advance_time.exe ${DATE} -${ASIM_WINDOW} 2>/dev/null)
export ASIM_DATE_MAX=$(${BUILD_DIR}/da_advance_time.exe ${DATE} +${ASIM_WINDOW} 2>/dev/null)
export ASIM_MN_YYYY=$(echo $ASIM_DATE_MIN | cut -c1-4)
export ASIM_MN_MM=$(echo $ASIM_DATE_MIN | cut -c5-6)
export ASIM_MN_DD=$(echo $ASIM_DATE_MIN | cut -c7-8)
export ASIM_MN_HH=$(echo $ASIM_DATE_MIN | cut -c9-10)
#
export ASIM_MX_YYYY=$(echo $ASIM_DATE_MAX | cut -c1-4)
export ASIM_MX_MM=$(echo $ASIM_DATE_MAX | cut -c5-6)
export ASIM_MX_DD=$(echo $ASIM_DATE_MAX | cut -c7-8)
export ASIM_MX_HH=$(echo $ASIM_DATE_MAX | cut -c9-10)
#
# WRFCHEM FIRE PARAMETERS:
export FIRE_START_DATE=${YYYY}-${MM}-${DD}
export E_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} ${FCST_PERIOD} 2>/dev/null)
export E_YYYY=$(echo $E_DATE | cut -c1-4)
export E_MM=$(echo $E_DATE | cut -c5-6)
export E_DD=$(echo $E_DATE | cut -c7-8)
export E_HH=$(echo $E_DATE | cut -c9-10)
export FIRE_END_DATE=${E_YYYY}-${E_MM}-${E_DD}
#
#
#########################################################################
#
# CREATE RUN DIRECTORY
#
#########################################################################
#
if [[ ! -e ${RUN_DIR} ]]; then mkdir ${RUN_DIR}; fi
cd ${RUN_DIR}
#
#########################################################################
#
# RUN GEOGRID
#
#########################################################################
#
if [[ ${RUN_GEOGRID} = "true" ]] then
   mkdir -p ${RUN_DIR}/geogrid
   cd ${RUN_DIR}/geogrid
#
   cp ${WPS_DIR}/geogrid.exe ./.
   export NL_DX=${DX_CR}
   export NL_DY=${DX_CR}
#   export NL_START_DATE=${FILE_DATE}
#   export NL_END_DATE=${NEXT_FILE_DATE}
   ${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist_RT.ksh
#
   if [[ -f job.ksh ]]; then rm -rf job.ksh; fi
   touch job.ksh
   RANDOM=$$
   export JOBRND=geogrid_$RANDOM
   rm -rf *.jerr
   rm -rf *.jout
   cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -n 1                                          # number of total (MPI) tasks
#BSUB -R "span[ptile=${GEOGRID_TASKS_PER_NODE}]"    # mpi tasks per node
#BSUB -J ${JOBRND}                                  # job name
#BSUB -o ${JOBRND}.jout                             # output filename
#BSUB -e ${JOBRND}.jerr                             # error filename
#BSUB -W ${GEOGRID_TIME_LIMIT}                      # wallclock time (minutes)
#BSUB -q geyser 
#
mpirun.lsf ./geogrid.exe  > index.html 2>&1 
#
export RC=\$?     
if [[ -f SUCCESS ]]; then rm -rf SUCCESS; fi     
if [[ -f FAILED ]]; then rm -rf FAILED; fi          
if [[ \$RC = 0 ]]; then
   touch SUCCESS
else
   touch FAILED 
   exit
fi
EOF
   bsub -K < job.ksh 
fi
#
#########################################################################
#
# RUN UNGRIB
#
#########################################################################
#
if [[ ${RUN_UNGRIB} = "true" ]]; then 
   mkdir -p ${RUN_DIR}/${DATE}/ungrib
   cd ${RUN_DIR}/${DATE}/ungrib
   rm -rf GRIBFILE.*
#
   cp ${VTABLE_DIR}/Vtable.${VTABLE_TYPE} Vtable
   cp ${WPS_DIR}/ungrib.exe ./.
#
   export L_FCST_RANGE=${LBC_END}
   export L_START_DATE=${DATE}
   export L_END_DATE=$($BUILD_DIR/da_advance_time.exe ${L_START_DATE} ${L_FCST_RANGE} 2>/dev/null)
   export L_START_YEAR=$(echo $L_START_DATE | cut -c1-4)
   export L_START_MONTH=$(echo $L_START_DATE | cut -c5-6)
   export L_START_DAY=$(echo $L_START_DATE | cut -c7-8)
   export L_START_HOUR=$(echo $L_START_DATE | cut -c9-10)
   export L_END_YEAR=$(echo $L_END_DATE | cut -c1-4)
   export L_END_MONTH=$(echo $L_END_DATE | cut -c5-6)
   export L_END_DAY=$(echo $L_END_DATE | cut -c7-8)
   export L_END_HOUR=$(echo $L_END_DATE | cut -c9-10)
   export NL_START_YEAR=$(echo $L_START_DATE | cut -c1-4),$(echo $L_START_DATE | cut -c1-4)
   export NL_START_MONTH=$(echo $L_START_DATE | cut -c5-6),$(echo $L_START_DATE | cut -c5-6)
   export NL_START_DAY=$(echo $L_START_DATE | cut -c7-8),$(echo $L_START_DATE | cut -c7-8)
   export NL_START_HOUR=$(echo $L_START_DATE | cut -c9-10),$(echo $L_START_DATE | cut -c9-10)
   export NL_END_YEAR=$(echo $L_END_DATE | cut -c1-4),$(echo $L_END_DATE | cut -c1-4)
   export NL_END_MONTH=$(echo $L_END_DATE | cut -c5-6),$(echo $L_END_DATE | cut -c5-6)
   export NL_END_DAY=$(echo $L_END_DATE | cut -c7-8),$(echo $L_END_DATE | cut -c7-8)
   export NL_END_HOUR=$(echo $L_END_DATE | cut -c9-10),$(echo $L_END_DATE | cut -c9-10)
   export NL_START_DATE="'"${L_START_YEAR}-${L_START_MONTH}-${L_START_DAY}_${L_START_HOUR}:00:00"'","'"${L_START_YEAR}-${L_START_MONTH}-${L_START_DAY}_${L_START_HOUR}:00:00"'"
   export NL_END_DATE="'"${L_END_YEAR}-${L_END_MONTH}-${L_END_DAY}_${L_END_HOUR}:00:00"'","'"${L_END_YEAR}-${L_END_MONTH}-${L_END_DAY}_${L_END_HOUR}:00:00"'"
   ${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist_RT.ksh
#
# UNTAR THE PARENT FORECAST FILES
   FILES=''
   if [[ -e ${GRIB_DIR}/${DATE} ]]; then
      if [[ -e ${GRIB_DIR}/${DATE}/${GRIB_PART1}${DATE}${GRIB_PART2} ]]; then
         cd ${GRIB_DIR}/${DATE}
         tar -xf ${GRIB_DIR}/${DATE}/${GRIB_PART1}${DATE}${GRIB_PART2}
         cd ${RUN_DIR}/${DATE}/ungrib
      else
         echo 'APM: ERROR - No GRIB files in directory'
         exit
      fi
#  
      if [[ ${SINGLE_FILE} == false ]]; then
         export CCHH=${HH}00
         (( LBC_ITR=${LBC_START} ))
         while [[ ${LBC_ITR} -le ${LBC_END} ]]; do
            if [[ ${LBC_ITR} -lt 1000 ]]; then export CFTM=${LBC_ITR}; fi
            if [[ ${LBC_ITR} -lt 100  ]]; then export CFTM=0${LBC_ITR}; fi
            if [[ ${LBC_ITR} -lt 10   ]]; then export CFTM=00${LBC_ITR}; fi
            if [[ ${LBC_ITR} -eq 0    ]]; then export CFTM=000; fi
            export FILE=${GRIB_DIR}/${DATE}/${GRIB_PART1}${START_YEAR}${START_MONTH}${START_DAY}_${CCHH}_${CFTM}.grb2
            FILES="${FILES} ${FILE}"
            (( LBC_ITR=${LBC_ITR}+${LBC_FREQ} ))
         done
      else
         export FILE=${GRIB_DIR}/${DATE}/GFS_Global_0p5deg_20080612_1800.grib2
         FILES="${FILES} ${FILE}"
      fi
   fi
#
# LINK GRIB FILES
   ${WPS_DIR}/link_grib.csh $FILES
#
# RUN UNGRIB
   ./ungrib.exe
   RC=$?
   if [[ $RC != 0 ]]; then
      echo ungrib failed with error $RC
      exit $RC
   fi
#
# TAR THE PARENT FORECAST FILES
   if [[ -e ${GRIB_DIR}/${DATE}/${GRIB_PART1}${DATE}${GRIB_PART2} ]]; then
      rm -rf ${GRIB_DIR}/${DATE}/${GRIB_PART1}*.grb2
   else
      cd ${GRIB_DIR}
      tar -cf ${GRIB_PART1}${DATE}${GRIB_PART2} ${DATE}
      mv ${GRIB_PART1}${DATE}${GRIB_PART2} ${DATE}/.
      if [[ -e ${DATE}/${GRIB_PART1}${DATE}${GRIB_PART2} ]]; then
         rm -rf ${DATE}/${GRIB_PART1}*.grb2
      else
         echo 'APM: Failed to created tar file'
         exit
      fi
      cd ${RUN_DIR}/${DATE}/ungrib
   fi
fi
#
#########################################################################
#
# RUN METGRID
#
#########################################################################
#
if [[ ${RUN_METGRID} = "true" ]]; then 
   mkdir -p ${RUN_DIR}/${DATE}/metgrid
   cd ${RUN_DIR}/${DATE}/metgrid
#
   ln -fs ${GEOGRID_DIR}/geo_em.d${CR_DOMAIN}.nc ./.
   ln -fs ${GEOGRID_DIR}/geo_em.d${FR_DOMAIN}.nc ./.
   ln -fs ../ungrib/FILE:* ./.
   ln -fs ${WPS_DIR}/metgrid/METGRID.TBL.${METGRID_TABLE_TYPE} METGRID.TBL
   ln -fs ${WPS_DIR}/metgrid.exe .
#
   export L_FCST_RANGE=${LBC_END}
   export L_START_DATE=${DATE}
   export L_END_DATE=$($BUILD_DIR/da_advance_time.exe ${L_START_DATE} ${L_FCST_RANGE} 2>/dev/null)
   export L_START_YEAR=$(echo $L_START_DATE | cut -c1-4)
   export L_START_MONTH=$(echo $L_START_DATE | cut -c5-6)
   export L_START_DAY=$(echo $L_START_DATE | cut -c7-8)
   export L_START_HOUR=$(echo $L_START_DATE | cut -c9-10)
   export L_END_YEAR=$(echo $L_END_DATE | cut -c1-4)
   export L_END_MONTH=$(echo $L_END_DATE | cut -c5-6)
   export L_END_DAY=$(echo $L_END_DATE | cut -c7-8)
   export L_END_HOUR=$(echo $L_END_DATE | cut -c9-10)
   export NL_START_YEAR=$(echo $L_START_DATE | cut -c1-4),$(echo $L_START_DATE | cut -c1-4)
   export NL_START_MONTH=$(echo $L_START_DATE | cut -c5-6),$(echo $L_START_DATE | cut -c5-6)
   export NL_START_DAY=$(echo $L_START_DATE | cut -c7-8),$(echo $L_START_DATE | cut -c7-8)
   export NL_START_HOUR=$(echo $L_START_DATE | cut -c9-10),$(echo $L_START_DATE | cut -c9-10)
   export NL_END_YEAR=$(echo $L_END_DATE | cut -c1-4),$(echo $L_END_DATE | cut -c1-4)
   export NL_END_MONTH=$(echo $L_END_DATE | cut -c5-6),$(echo $L_END_DATE | cut -c5-6)
   export NL_END_DAY=$(echo $L_END_DATE | cut -c7-8),$(echo $L_END_DATE | cut -c7-8)
   export NL_END_HOUR=$(echo $L_END_DATE | cut -c9-10),$(echo $L_END_DATE | cut -c9-10)
   export NL_START_DATE="'"${L_START_YEAR}-${L_START_MONTH}-${L_START_DAY}_${L_START_HOUR}:00:00"'","'"${L_START_YEAR}-${L_START_MONTH}-${L_START_DAY}_${L_START_HOUR}:00:00"'"
   export NL_END_DATE="'"${L_END_YEAR}-${L_END_MONTH}-${L_END_DAY}_${L_END_HOUR}:00:00"'","'"${L_END_YEAR}-${L_END_MONTH}-${L_END_DAY}_${L_END_HOUR}:00:00"'"
   ${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist_RT.ksh
#
# JOB SCRIPT
   if [[ -f job.ksh ]]; then rm -rf job.ksh; fi
   rm -rf *.err
   rm -rf *.out
   touch job.ksh
   RANDOM=$$
   export JOBRND=metgrid_$RANDOM
   cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -n 1                                  # number of total (MPI) tasks
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W 00:10                              # wallclock time (minutes)
#BSUB -q geyser
#
# Run metgrid
mpirun.lsf ./metgrid.exe  > index.html 2>&1 
#
export RC=\$?     
if [[ -f SUCCESS ]]; then rm -rf SUCCESS; fi     
if [[ -f FAILED ]]; then rm -rf FAILED; fi          
if [[ \$RC = 0 ]]; then
   touch SUCCESS
else
   touch FAILED 
   exit
fi
EOF
#
# Submit script and wait until job completes
   bsub -K < job.ksh 
fi
#
#########################################################################
#
# RUN REAL
#
#########################################################################
#
if [[ ${RUN_REAL} = "true" ]]; then 
   mkdir -p ${RUN_DIR}/${DATE}/real
   cd ${RUN_DIR}/${DATE}/real
#
   cp ${WRF_DIR}/main/real.exe ./.
   cp ${DART_DIR}/models/wrf_chem/run_scripts/hist_io_flds ./.
#
# LINK IN THE METGRID FILES
   export P_DATE=${DATE}
   export P_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${LBC_END} 2>/dev/null)
   while [[ ${P_DATE} -le ${P_END_DATE} ]] ; do
      export P_YYYY=$(echo $P_DATE | cut -c1-4)
      export P_MM=$(echo $P_DATE | cut -c5-6)
      export P_DD=$(echo $P_DATE | cut -c7-8)
      export P_HH=$(echo $P_DATE | cut -c9-10)
      export P_FILE_DATE=${P_YYYY}-${P_MM}-${P_DD}_${P_HH}:00:00.nc
      ln -sf ${RUN_DIR}/${DATE}/metgrid/met_em.d${CR_DOMAIN}.${P_FILE_DATE} ./.
      ln -sf ${RUN_DIR}/${DATE}/metgrid/met_em.d${FR_DOMAIN}.${P_FILE_DATE} ./.
      export P_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${LBC_FREQ} 2>/dev/null) 
   done
#
# LOOP THROUGH BDY TENDENCY TIMES FOR PERTURB_BC
   export P_DATE=${DATE}
   export P_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${FCST_PERIOD} 2>/dev/null)
   while [[ ${P_DATE} -le ${P_END_DATE} ]] ; do      
#
# CREATE WRF NAMELIST
      export L_FCST_RANGE=${FCST_PERIOD}
      export NL_DX=${DX_CR},${DX_FR}
      export NL_DY=${DX_CR},${DX_FR}
      export L_START_DATE=${P_DATE}
      export L_END_DATE=$($BUILD_DIR/da_advance_time.exe ${L_START_DATE} ${L_FCST_RANGE} 2>/dev/null)
      export L_START_YEAR=$(echo $L_START_DATE | cut -c1-4)
      export L_START_MONTH=$(echo $L_START_DATE | cut -c5-6)
      export L_START_DAY=$(echo $L_START_DATE | cut -c7-8)
      export L_START_HOUR=$(echo $L_START_DATE | cut -c9-10)
      export L_END_YEAR=$(echo $L_END_DATE | cut -c1-4)
      export L_END_MONTH=$(echo $L_END_DATE | cut -c5-6)
      export L_END_DAY=$(echo $L_END_DATE | cut -c7-8)
      export L_END_HOUR=$(echo $L_END_DATE | cut -c9-10)
      export NL_START_YEAR=$(echo $L_START_DATE | cut -c1-4),$(echo $L_START_DATE | cut -c1-4)
      export NL_START_MONTH=$(echo $L_START_DATE | cut -c5-6),$(echo $L_START_DATE | cut -c5-6)
      export NL_START_DAY=$(echo $L_START_DATE | cut -c7-8),$(echo $L_START_DATE | cut -c7-8)
      export NL_START_HOUR=$(echo $L_START_DATE | cut -c9-10),$(echo $L_START_DATE | cut -c9-10)
      export NL_END_YEAR=$(echo $L_END_DATE | cut -c1-4),$(echo $L_END_DATE | cut -c1-4)
      export NL_END_MONTH=$(echo $L_END_DATE | cut -c5-6),$(echo $L_END_DATE | cut -c5-6)
      export NL_END_DAY=$(echo $L_END_DATE | cut -c7-8),$(echo $L_END_DATE | cut -c7-8)
      export NL_END_HOUR=$(echo $L_END_DATE | cut -c9-10),$(echo $L_END_DATE | cut -c9-10)
      export NL_START_DATE="'"${L_START_YEAR}-${L_START_MONTH}-${L_START_DAY}_${L_START_HOUR}:00:00"'","'"${L_START_YEAR}-${L_START_MONTH}-${L_START_DAY}_${L_START_HOUR}:00:00"'"
      export NL_END_DATE="'"${L_END_YEAR}-${L_END_MONTH}-${L_END_DAY}_${L_END_HOUR}:00:00"'","'"${L_END_YEAR}-${L_END_MONTH}-${L_END_DAY}_${L_END_HOUR}:00:00"'"
      ${HYBRID_SCRIPTS_DIR}/da_create_wrfchem_namelist_nested_RT.ksh
#
# JOB SCRIPT
      if [[ -f job.ksh ]]; then rm -rf job.ksh; fi
      rm -rf *.err
      rm -rf *.out
      touch job.ksh
      RANDOM=$$
      export JOBRND=real_$RANDOM
      cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -n 1                                  # number of total (MPI) tasks
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W 00:10                              # wallclock time (minutes)
#BSUB -q geyser
#
# Run real
mpirun.lsf ./real.exe  > index.html 2>&1 
#
export RC=\$?     
if [[ -f SUCCESS ]]; then rm -rf SUCCESS; fi     
if [[ -f FAILED ]]; then rm -rf FAILED; fi          
if [[ \$RC = 0 ]]; then
   touch SUCCESS
else
   touch FAILED 
   exit
fi
EOF
#
      bsub -K < job.ksh 
#
      mv wrfinput_d${CR_DOMAIN} wrfinput_d${CR_DOMAIN}_$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} 0 -W 2>/dev/null)
      mv wrfinput_d${FR_DOMAIN} wrfinput_d${FR_DOMAIN}_$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} 0 -W 2>/dev/null)
      mv wrfbdy_d${CR_DOMAIN} wrfbdy_d${CR_DOMAIN}_$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} 0 -W 2>/dev/null)
#      mv wrfbdy_d${FR_DOMAIN} wrfbdy_d${FR_DOMAIN}_$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} 0 -W 2>/dev/null)
      export P_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${LBC_FREQ} 2>/dev/null) 
   done   
fi
#
#########################################################################
#
# RUN PERT_WRFCHEM_MET_IC
#
#########################################################################
#
if [[ ${RUN_PERT_WRFCHEM_MET_IC} = "true" ]]; then 
   if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_met_ic ]]; then
      mkdir -p ${RUN_DIR}/${DATE}/wrfchem_met_ic
      cd ${RUN_DIR}/${DATE}/wrfchem_met_ic
   else
      cd ${RUN_DIR}/${DATE}/wrfchem_met_ic
   fi
#
   export NL_MAX_DOM=1
   export NL_OB_FORMAT=1
   export L_FCST_RANGE=${FCST_PERIOD}
   export L_START_DATE=${DATE}
   export L_END_DATE=$($BUILD_DIR/da_advance_time.exe ${L_START_DATE} ${L_FCST_RANGE} 2>/dev/null)
   export L_START_YEAR=$(echo $L_START_DATE | cut -c1-4)
   export L_START_MONTH=$(echo $L_START_DATE | cut -c5-6)
   export L_START_DAY=$(echo $L_START_DATE | cut -c7-8)
   export L_START_HOUR=$(echo $L_START_DATE | cut -c9-10)
   export L_END_YEAR=$(echo $L_END_DATE | cut -c1-4)
   export L_END_MONTH=$(echo $L_END_DATE | cut -c5-6)
   export L_END_DAY=$(echo $L_END_DATE | cut -c7-8)
   export L_END_HOUR=$(echo $L_END_DATE | cut -c9-10)
   export NL_START_YEAR=$(echo $L_START_DATE | cut -c1-4),$(echo $L_START_DATE | cut -c1-4)
   export NL_START_MONTH=$(echo $L_START_DATE | cut -c5-6),$(echo $L_START_DATE | cut -c5-6)
   export NL_START_DAY=$(echo $L_START_DATE | cut -c7-8),$(echo $L_START_DATE | cut -c7-8)
   export NL_START_HOUR=$(echo $L_START_DATE | cut -c9-10),$(echo $L_START_DATE | cut -c9-10)
   export NL_END_YEAR=$(echo $L_END_DATE | cut -c1-4),$(echo $L_END_DATE | cut -c1-4)
   export NL_END_MONTH=$(echo $L_END_DATE | cut -c5-6),$(echo $L_END_DATE | cut -c5-6)
   export NL_END_DAY=$(echo $L_END_DATE | cut -c7-8),$(echo $L_END_DATE | cut -c7-8)
   export NL_END_HOUR=$(echo $L_END_DATE | cut -c9-10),$(echo $L_END_DATE | cut -c9-10)
   export NL_START_DATE="'"${L_START_YEAR}-${L_START_MONTH}-${L_START_DAY}_${L_START_HOUR}:00:00"'","'"${L_START_YEAR}-${L_START_MONTH}-${L_START_DAY}_${L_START_HOUR}:00:00"'"
   export NL_END_DATE="'"${L_END_YEAR}-${L_END_MONTH}-${L_END_DAY}_${L_END_HOUR}:00:00"'","'"${L_END_YEAR}-${L_END_MONTH}-${L_END_DAY}_${L_END_HOUR}:00:00"'"
#
# LOOP THROUGH ALL BDY TENDENCY TIMES
   export P_DATE=${DATE}
   export P_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${FCST_PERIOD} 2>/dev/null)
   while [[ ${P_DATE} -le ${P_END_DATE} ]] ; do
#
# SET WRFDA PARAMETERS
      export ANALYSIS_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} 0 -W 2>/dev/null)
      export NL_ANALYSIS_DATE="'"${ANALYSIS_DATE}"'"
      export NL_TIME_WINDOW_MIN="'"$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} -${ASIM_WINDOW} -W 2>/dev/null)"'"
      export NL_TIME_WINDOW_MAX="'"$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} +${ASIM_WINDOW} -W 2>/dev/null)"'"
      export NL_ANALYSIS_TYPE="'"RANDOMCV"'"
      export NL_PUT_RAND_SEED=true
#
# LOOP THROUGH ALL MEMBERS IN THE ENSEMBLE
      export RAN_APM=${RANDOM}
      let MEM=1
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
#
# COARSE RESOLUTION GRID
         cd ${RUN_DIR}/${DATE}/wrfchem_met_ic
         export LCR_DIR=wrfda_cr_${MEM}
         if [[ ! -e ${LCR_DIR} ]]; then
            mkdir ${LCR_DIR}
            cd ${LCR_DIR}
         else
            cd ${LCR_DIR}
            rm -rf *
         fi
         export NL_E_WE=${NNXP_CR}
         export NL_E_SN=${NNYP_CR}
         export NL_DX=${DX_CR}
         export NL_DY=${DX_CR}
         export NL_GRID_ID=1
         export NL_PARENT_ID=1
         export NL_PARENT_GRID_RATIO=1
         export NL_I_PARENT_START=${ISTR_CR}
         export NL_J_PARENT_START=${JSTR_CR}
         export DA_INPUT_FILE=../../real/wrfinput_d${CR_DOMAIN}_${ANALYSIS_DATE}
         export NL_SEED_ARRAY1=$(${BUILD_DIR}/da_advance_time.exe ${DATE} 0 -f hhddmmyycc)
         export NL_SEED_ARRAY2=`echo ${MEM} \* 100000 | bc -l `
         ${HYBRID_SCRIPTS_DIR}/da_create_wrfda_namelist.ksh
         cp ${OB_PREPBUFR_DIR}/${DATE}/prepbufr.gdas.${DATE}.wo40.be ob.bufr
         cp ${DA_INPUT_FILE} fg
         cp ${BE_DIR}/be.dat.cv3 be.dat
         cp ${WRFVAR_DIR}/run/LANDUSE.TBL ./.
         cp ${WRFVAR_DIR}/var/da/da_wrfvar.exe ./.
#
# JOB SCRIPT 
         if [[ -f job.ksh ]]; then rm -rf job.ksh; fi
         rm -rf *.err
         rm -rf *.out
         touch job.ksh
         export JOBRND=wrfda_${RAN_APM}
         cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -x                                    # exclusive use of node (not_shared)
#BSUB -n ${WRFDA_NUM_TASKS}                       # number of total (MPI) tasks
#BSUB -R "span[ptile=${WRFDA_TASKS_PER_NODE}]"    # mpi tasks per node
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W ${WRFDA_TIME_LIMIT}               # wallclock time (minutes)
#BSUB -q ${WRFDA_JOB_CLASS}
#
mpirun.lsf ./da_wrfvar.exe  > index.html 2>&1 
#
export RC=\$?     
if [[ -f SUCCESS ]]; then rm -rf SUCCESS_*; fi     
if [[ -f FAILED ]]; then rm -rf FAILED_*; fi          
if [[ \$RC = 0 ]]; then
   touch SUCCESS_${RAN_APM}
else
   touch FAILED_${RAN_APM} 
   exit
fi
EOF
#
         bsub < job.ksh 
#
# FINE RESOLUTION GRID
         cd ${RUN_DIR}/${DATE}/wrfchem_met_ic
         export LFR_DIR=wrfda_fr_${MEM}
         if [[ ! -e ${LFR_DIR} ]]; then
            mkdir ${LFR_DIR}
            cd ${LFR_DIR}
         else
            cd ${LFR_DIR}
            rm -rf *
         fi
         export NL_E_WE=${NNXP_FR}
         export NL_E_SN=${NNYP_FR}
         export NL_DX=${DX_FR}
         export NL_DY=${DX_FR}
         export NL_GRID_ID=2
         export NL_PARENT_ID=1
         export NL_PARENT_GRID_RATIO=5
         export NL_I_PARENT_START=${ISTR_FR}
         export NL_J_PARENT_START=${JSTR_FR}
         export DA_INPUT_FILE=../../real/wrfinput_d${FR_DOMAIN}_${ANALYSIS_DATE}
         ${HYBRID_SCRIPTS_DIR}/da_create_wrfda_namelist.ksh
         cp ${OB_PREPBUFR_DIR}/${DATE}/prepbufr.gdas.${DATE}.wo40.be ob.bufr
         cp ${DA_INPUT_FILE} fg
         cp ${BE_DIR}/be.dat.cv3 be.dat
         cp ${WRFVAR_DIR}/run/LANDUSE.TBL ./.
         cp ${WRFVAR_DIR}/var/da/da_wrfvar.exe ./.
#
# JOB SCRIPT 
         if [[ -f job.ksh ]]; then rm -rf job.ksh; fi
         rm -rf *.err
         rm -rf *.out
         touch job.ksh
         export JOBRND=wrfda_${RAN_APM}
         cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -x                                    # exclusive use of node (not_shared)
#BSUB -n ${WRFDA_NUM_TASKS}                       # number of total (MPI) tasks
#BSUB -R "span[ptile=${WRFDA_TASKS_PER_NODE}]"    # mpi tasks per node
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W ${WRFDA_TIME_LIMIT}               # wallclock time (minutes)
#BSUB -q ${WRFDA_JOB_CLASS}
#
mpirun.lsf ./da_wrfvar.exe  > index.html 2>&1 
#
export RC=\$?     
if [[ -f SUCCESS ]]; then rm -rf SUCCESS_*; fi     
if [[ -f FAILED ]]; then rm -rf FAILED_*; fi          
if [[ \$RC = 0 ]]; then
   touch SUCCESS_${RAN_APM}
else
   touch FAILED_${RAN_APM}
   exit
fi
EOF
#
         bsub < job.ksh 
         let MEM=${MEM}+1
      done
#
# Wait for WRFDA to complete for all members
      cd ${RUN_DIR}/${DATE}/wrfchem_met_ic
      ${HYBRID_SCRIPTS_DIR}/da_run_hold.ksh ${RAN_APM}
#
      let MEM=1
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
         export LCR_DIR=wrfda_cr_${MEM}
         export LFR_DIR=wrfda_fr_${MEM}
         cp ${LCR_DIR}/wrfvar_output wrfinput_d${CR_DOMAIN}_${ANALYSIS_DATE}.${CMEM}
         cp ${LFR_DIR}/wrfvar_output wrfinput_d${FR_DOMAIN}_${ANALYSIS_DATE}.${CMEM}
         let MEM=${MEM}+1
      done
      export P_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${LBC_FREQ} 2>/dev/null) 
   done
   export NL_E_WE=${NNXP_CR},${NNXP_FR}
   export NL_E_SN=${NNYP_CR},${NNYP_FR}
   export NL_DX=${DX_CR},${DX_FR}
   export NL_DY=${DX_CR},${DX_FR}
   export NL_GRID_ID=1,2
   export NL_PARENT_ID=1,1
   export NL_PARENT_GRID_RATIO=1,5
   export NL_I_PARENT_START=${ISTR_CR},${ISTR_FR}
   export NL_J_PARENT_START=${JSTR_CR},${JSTR_FR}
fi
#
#########################################################################
#
# RUN PERT WRFCHEM MET BC
#
#########################################################################
#
if [[ ${RUN_PERT_WRFCHEM_MET_BC} = "true" ]]; then 
   if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_met_bc ]]; then
      mkdir -p ${RUN_DIR}/${DATE}/wrfchem_met_bc
      cd ${RUN_DIR}/${DATE}/wrfchem_met_bc
   else
      cd ${RUN_DIR}/${DATE}/wrfchem_met_bc
   fi
#
# LOOP THROUGH ALL MEMBERS IN THE ENSEMBLE
   let MEM=1
   while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
      export CMEM=e${MEM}
      if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
      if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
      export ANALYSIS_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} 0 -W 2>/dev/null)
      if [[ -f wrfbdy_this ]]; then
         rm -rf wrfbdy_this
         export DA_BDY_PATH=${RUN_DIR}/${DATE}/real
         export DA_BDY_FILE=${DA_BDY_PATH}/wrfbdy_d${CR_DOMAIN}_${ANALYSIS_DATE}
         cp ${DA_BDY_FILE} wrfbdy_this
      else
         export DA_BDY_PATH=${RUN_DIR}/${DATE}/real
         export DA_BDY_FILE=${DA_BDY_PATH}/wrfbdy_d${CR_DOMAIN}_${ANALYSIS_DATE}
         cp ${DA_BDY_FILE} wrfbdy_this
      fi
      rm -rf pert_wrf_bc
      cp ${TRUNK_DIR}/${DART_VER}/models/wrf_chem/work/pert_wrf_bc ./.
      rm -rf input.nml
      ${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_input.nml.ksh
#
# LOOP THROUGH ALL BDY TENDENCY TIMES FOR THIS MEMBER.
      export L_DATE=${DATE}
      export L_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${FCST_PERIOD} 2>/dev/null)
      while [[ ${L_DATE} -lt ${L_END_DATE} ]]; do
         export ANALYSIS_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 0 -W 2>/dev/null)
         export NEXT_L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${LBC_FREQ} 2>/dev/null)
         export NEXT_ANALYSIS_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${LBC_FREQ} -W 2>/dev/null)
         rm -rf wrfinput_this
         rm -rf wrfinput_next
         export DA_INPUT_PATH=${RUN_DIR}/${DATE}/wrfchem_met_ic
         ln -fs ${DA_INPUT_PATH}/wrfinput_d${CR_DOMAIN}_${ANALYSIS_DATE}.${CMEM} wrfinput_this 
         ln -fs ${DA_INPUT_PATH}/wrfinput_d${CR_DOMAIN}_${NEXT_ANALYSIS_DATE}.${CMEM} wrfinput_next 
#
# JOB SCRIPT 
         if [[ -f job.ksh ]]; then rm -rf job.ksh; fi
         rm -rf *.err
         rm -rf *.out
         touch job.ksh
         RANDOM=$$
         export JOBRND=pert_bc_$RANDOM
         cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -n 1                                  # number of total (MPI) tasks
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W 00:10                              # wallclock time (minutes)
#BSUB -q geyser
#
mpirun.lsf ./pert_wrf_bc  > index.html 2>&1 
#
export RC=\$?     
if [[ -f SUCCESS ]]; then rm -rf SUCCESS; fi     
if [[ -f FAILED ]]; then rm -rf FAILED; fi          
if [[ \$RC = 0 ]]; then
   touch SUCCESS
else
   touch FAILED 
   exit
fi
EOF
#
         bsub -K < job.ksh 
         RC=$?
         if [[ $RC != 0 ]]; then
            echo pert_wrf_bc failed with error $RC
            exit $RC
         fi
         export L_DATE=${NEXT_L_DATE}
      done
      export ANALYSIS_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} 0 -W 2>/dev/null)
      mv wrfbdy_this wrfbdy_d${CR_DOMAIN}_${ANALYSIS_DATE}.${CMEM}
      let MEM=${MEM}+1
   done
fi
#
#########################################################################
#
# RUN EXO_COLDENS
#
#########################################################################
#
if ${RUN_EXO_COLDENS}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/exo_coldens ]]; then
      mkdir ${RUN_DIR}/${DATE}/exo_coldens
      cd ${RUN_DIR}/${DATE}/exo_coldens
   else
      cd ${RUN_DIR}/${DATE}/exo_coldens
   fi
#
# LINK NEEDED FILES
   export FILE_CR=wrfinput_d${CR_DOMAIN}
   export FILE_FR=wrfinput_d${FR_DOMAIN}
   rm -rf ${FILE_CR}
   rm -rf ${FILE_FR}
   ln -sf ${REAL_DIR}/${FILE_CR}_${FILE_DATE} ${FILE_CR}   
   ln -sf ${REAL_DIR}/${FILE_FR}_${FILE_DATE} ${FILE_FR}   
   export FILE=exo_coldens.nc
   rm -rf ${FILE}
   ln -sf ${COLDENS_INPUT_DIR}/${FILE} ${FILE}
   export FILE=exo_coldens
   rm -rf ${FILE}
   ln -sf ${COLDENS_INPUT_DIR}/${FILE} ${FILE}
#
# CREATE INPUT FILE
   export FILE=exo_coldens.inp
   rm -rf ${FILE}
   cat << EOF > ${FILE}
&control
domains = 2,
/
EOF
#
# RUN exo_coldens
   ./exo_coldens < exo_coldens.inp
#
# TEST WHETHER OUTPUT EXISTS
   export FILE_CR=exo_coldens_d${CR_DOMAIN}
   export FILE_FR=exo_coldens_d${FR_DOMAIN}
   if [[ ! -e ${FILE_CR} || ! -e ${FILE_FR} ]]; then
      echo EXO_COLDENS FAILED
      exit
   else
      echo EXO_COLDENS SUCCESS
   fi
fi
#
#########################################################################
#
# RUN SEASONS_WES
#
#########################################################################
#
if ${RUN_SEASON_WES}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/seasons_wes ]]; then
      mkdir ${RUN_DIR}/${DATE}/seasons_wes
      cd ${RUN_DIR}/${DATE}/seasons_wes
   else
      cd ${RUN_DIR}/${DATE}/seasons_wes
   fi
#
# LINK NEEDED FILES
   export FILE_CR=wrfinput_d${CR_DOMAIN}
   export FILE_FR=wrfinput_d${FR_DOMAIN}
   rm -rf ${FILE_CR}
   rm -rf ${FILE_FR}
   ln -sf ${REAL_DIR}/${FILE_CR}_${FILE_DATE} ${FILE_CR}   
   ln -sf ${REAL_DIR}/${FILE_FR}_${FILE_DATE} ${FILE_FR}   
   export FILE=season_wes_usgs.nc
   rm -rf ${FILE}
   ln -sf ${COLDENS_INPUT_DIR}/${FILE} ${FILE}
   export FILE=wesely
   rm -rf ${FILE}
   ln -sf ${COLDENS_INPUT_DIR}/${FILE} ${FILE}
#
# CREATE INPUT FILE
   export FILE=wesely.inp
   rm -rf ${FILE}
   cat << EOF > ${FILE}
&control
domains = 2,
/
EOF
#
# RUN wesely
   ./wesely < wesely.inp
#
# TEST WHETHER OUTPUT EXISTS
   export FILE_CR=wrf_season_wes_usgs_d${CR_DOMAIN}.nc
   export FILE_FR=wrf_season_wes_usgs_d${FR_DOMAIN}.nc
   if [[ ! -e ${FILE_CR} || ! -e ${FILE_FR} ]]; then
      echo WESELY FAILED
      exit
   else
      echo WESELY SUCCESS
   fi
fi
#
#########################################################################
#
# RUN WRFCHEM_BIO
#
#########################################################################
#
if ${RUN_WRFCHEM_BIO}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_bio ]]; then
      mkdir ${RUN_DIR}/${DATE}/wrfchem_bio
      cd ${RUN_DIR}/${DATE}/wrfchem_bio
   else
      cd ${RUN_DIR}/${DATE}/wrfchem_bio
   fi
#
# LOOP THROUGHT CURRENT AND NEXT DATE
   export L_DATE=${DATE}
   export LE_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${FCST_PERIOD} 2>/dev/null)
   while [[ ${L_DATE} -le ${LE_DATE} ]]; do 
      export L_YYYY=$(echo $L_DATE | cut -c1-4)
      export L_MM=$(echo $L_DATE | cut -c5-6)
      export L_DD=$(echo $L_DATE | cut -c7-8)
      export L_HH=$(echo $L_DATE | cut -c9-10)
      export L_FILE_DATE=${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
#
# LINK NEEDED FILES
      export FILE_CR=wrfinput_d${CR_DOMAIN}
      export FILE_FR=wrfinput_d${FR_DOMAIN}
      rm -rf ${FILE_CR}
      rm -rf ${FILE_FR}
      cp ${REAL_DIR}/${FILE_CR}_${L_FILE_DATE} ${FILE_CR}   
      cp ${REAL_DIR}/${FILE_FR}_${L_FILE_DATE} ${FILE_FR}   
      export FILE_CR=wrfbiochemi_d${CR_DOMAIN}
      export FILE_FR=wrfbiochemi_d${FR_DOMAIN}
      if [[ ${L_DATE} -eq ${DATE} ]]; then
         rm -rf ${FILE_CR}
         rm -rf ${FILE_FR}
      fi
      rm -rf btr*.nc
      rm -rf DSW*.nc
      rm -rf hrb*.nc
      rm -rf iso*.nc
      rm -rf lai*.nc
      rm -rf ntr*.nc
      rm -rf shr*.nc
      rm -rf TAS*.nc
      cp ${MEGAN_INPUT_DIR}/*.nc ./.
      export FILE=megan_bio_emiss
      rm -rf ${FILE}
      cp ${MEGAN_BIO_DIR}/${FILE} ${FILE}
#
# CREATE INPUT FILE
      export FILE=megan_bio_emiss.inp
      rm -rf ${FILE}
      cat << EOF > ${FILE}
&control
domains = 2,
start_lai_mnth = 1,
end_lai_mnth = 12
/
EOF
#
# CREATE job.ksh
      rm -rf job.ksh
      rm -rf wrf_bio_*.err
      rm -rf wrf_bio_*.out
      rm -rf core.*
      touch job.ksh
      RANDOM=$$
      export JOBRND=wrf_bio_$RANDOM
      cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -n 1
#BSUB -R "span[ptile=8]"    
#BSUB -J ${JOBRND}
#BSUB -o ${JOBRND}.out
#BSUB -e ${JOBRND}.err
#BSUB -W 00:15        
#BSUB -q premium
#
# RUN megan_bio_emis
./megan_bio_emiss < megan_bio_emiss.inp > index_megan_bio 2>&1
# 
export RC=\$?     
rm -rf WRFCHEM_BIO_SUCCESS
rm -rf WRFCHEM_BIO_FAILED          
if [[ \$RC = 0 ]]; then
   touch WRFCHEM_BIO_SUCCESS
else
   touch WRFCHEM_BIO_FAILED 
   exit
fi
EOF
#
      bsub -K < job.ksh 
#
# TEST WHETHER OUTPUT EXISTS
      export FILE_CR=wrfbiochemi_d${CR_DOMAIN}
      export FILE_FR=wrfbiochemi_d${FR_DOMAIN}
      if [[ ! -e ${FILE_CR} || ! -e ${FILE_FR} ]]; then
         echo WRFCHEM_BIO FAILED
         exit
      else
         echo WRFCHEM_BIO SUCCESS
         mv ${FILE_CR} ${FILE_CR}_${L_FILE_DATE}
         mv ${FILE_FR} ${FILE_FR}_${L_FILE_DATE}
      fi
      export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 6 2>/dev/null)
   done
fi
#
#########################################################################
#
# RUN WRFCHEM_FIRE
#
#########################################################################
#
if ${RUN_WRFCHEM_FIRE}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_fire ]]; then
      mkdir ${RUN_DIR}/${DATE}/wrfchem_fire
      cd ${RUN_DIR}/${DATE}/wrfchem_fire
   else
      cd ${RUN_DIR}/${DATE}/wrfchem_fire
   fi
#
# LINK NEEDED FILES
   export FILE_CR=wrfinput_d${CR_DOMAIN}
   export FILE_FR=wrfinput_d${FR_DOMAIN}
   rm -rf ${FILE_CR}
   rm -rf ${FILE_FR}
   ln -sf ${REAL_DIR}/${FILE_CR}_${FILE_DATE} ${FILE_CR}   
   ln -sf ${REAL_DIR}/${FILE_FR}_${FILE_DATE} ${FILE_FR}   
   rm -rf GLOB_*.txt
   ln -sf ${FIRE_INPUT_DIR}/GLOB_*.txt ./.
   export FILE=fire_emis
   rm -rf ${FILE}
   ln -sf ${FIRE_INPUT_DIR}/src/${FILE} ${FILE}
   rm -rf grass_from_img.nc
   rm -rf shrub_from_img.nc
   rm -rf tempfor_from_img.nc
   rm -rf tropfor_from_img.nc
   ln -sf ${FIRE_INPUT_DIR}/grass_from_img.nc
   ln -sf ${FIRE_INPUT_DIR}/shrub_from_img.nc
   ln -sf ${FIRE_INPUT_DIR}/tempfor_from_img.nc
   ln -sf ${FIRE_INPUT_DIR}/tropfor_from_img.nc
#
# CREATE INPUT FILE
   export FILE=fire_emis.mozc.inp
   rm -rf ${FILE}
   cat << EOF > ${FILE}
&control
domains = 2,
fire_filename(1) = 'GLOB_MAY-SEPT2008_05062011_MOZ4.txt',
start_date = '${FIRE_START_DATE}', 
end_date = '${FIRE_END_DATE}',
fire_directory = './',
wrf_directory = './',
wrf2fire_map = 'co -> CO', 'no -> NO', 'so2 -> SO2', 'bigalk -> BIGALK',
               'bigene -> BIGENE', 'c2h4 -> C2H4', 'c2h5oh -> C2H5OH',
               'c2h6 -> C2H6', 'c3h8 -> C3H8','c3h6 -> C3H6','ch2o -> CH2O', 'ch3cho -> CH3CHO',
               'ch3coch3 -> CH3COCH3','ch3oh -> CH3OH','mek -> MEK','toluene -> TOLUENE',
               'nh3 -> NH3','no2 -> NO2','open -> BIGALD','c10h16 -> C10H16',
               'ch3cooh -> CH3COOH','cres -> CRESOL','glyald -> GLYALD','mgly -> CH3COCHO',
               'gly -> CH3COCHO','acetol -> HYAC','isop -> ISOP','macr -> MACR'
               'mvk -> MVK',
               'oc -> OC;aerosol','bc -> BC;aerosol'
/
EOF
#
# JOB SCRIPT
   if [[ -f job.ksh ]]; then rm -rf job.ksh; fi
   touch job.ksh
   RANDOM=$$
   export JOBRND=wrf_fire_$RANDOM
   cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -n 1
#BSUB -R "span[ptile=8]"    
#BSUB -J ${JOBRND}
#BSUB -o ${JOBRND}.out
#BSUB -e ${JOBRND}.err
#BSUB -W 00:20        
#BSUB -q small
#
./fire_emis < fire_emis.mozc.inp > index_fire_emis 2>&1
# 
export RC=\$?     
rm -rf WRFCHEM_FIRE_SUCCESS     
rm -rf WRFCHEM_FIRE_FAILED          
if [[ \$RC = 0 ]]; then
   touch WRFCHEM_FIRE_SUCCESS
else
   touch WRFCHEM_FIRE_FAILED 
   exit
fi
EOF
#
   bsub -K < job.ksh 
   export L_DATE=${DATE}
   while [[ ${L_DATE} -le ${END_DATE} ]]; do
      export L_YYYY=$(echo $L_DATE | cut -c1-4)
      export L_MM=$(echo $L_DATE | cut -c5-6)
      export L_DD=$(echo $L_DATE | cut -c7-8)
      export L_HH=$(echo $L_DATE | cut -c9-10)
      export L_FILE_DATE=${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
      export DD_DATE=${L_YYYY}${L_MM}${L_DD}
#
# TEST WHETHER OUTPUT EXISTS
      export FILE_CR=wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE}
      export FILE_FR=wrffirechemi_d${FR_DOMAIN}_${L_FILE_DATE}
      if [[ ! -e ${FILE_CR} || ! -e ${FILE_CR} ]]; then
         echo WRFFIRE FAILED
         exit
      else
         echo WRFFIRE SUCCESS
      fi
      export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 1 2>/dev/null)
   done
fi
#
#########################################################################
#
# RUN WRFCHEM_CHEMI
#
#########################################################################
#
if ${RUN_WRFCHEM_CHEMI}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_chemi ]]; then
      mkdir ${RUN_DIR}/${DATE}/wrfchem_chemi
      cd ${RUN_DIR}/${DATE}/wrfchem_chemi
   else
      cd ${RUN_DIR}/${DATE}/wrfchem_chemi
   fi
   export L_DATE=${DATE}
   export LE_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${FCST_PERIOD} 2>/dev/null)
#
   while [[ ${L_DATE} -le ${LE_DATE} ]]; do
      export L_YYYY=$(echo $L_DATE | cut -c1-4)
      export L_MM=$(echo $L_DATE | cut -c5-6)
      export L_DD=$(echo $L_DATE | cut -c7-8)
      export L_HH=$(echo $L_DATE | cut -c9-10)
#
      export FILE_PATH=${ACD_DIR}/AVE_TEST_DATA/chem_static_RT_test/${L_YYYY}${L_MM}${L_DD}
      cp ${FILE_PATH}/wrfchemi_d${CR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00 ./.
      cp ${FILE_PATH}/wrfchemi_d${FR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00 ./.
      export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 1 2>/dev/null)
   done
fi
#
#########################################################################
#
# RUN WRFCHEM PERTURB ICBC
#
#########################################################################
#
if ${RUN_PERT_WRFCHEM_CHEM_ICBC}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_chem_icbc ]]; then
      mkdir ${RUN_DIR}/${DATE}/wrfchem_chem_icbc
      cd ${RUN_DIR}/${DATE}/wrfchem_chem_icbc
   else
      cd ${RUN_DIR}/${DATE}/wrfchem_chem_icbc
   fi
#
# PERTURB CHEM ICBC
   cp ${PERT_CHEM_INPUT_DIR}/runICBC_parent_rt_CR.ksh ./.
   cp ${PERT_CHEM_INPUT_DIR}/runICBC_parent_rt_FR.ksh ./.
   cp ${PERT_CHEM_INPUT_DIR}/runICBC_setN_rt_CR.ksh ./.
   cp ${PERT_CHEM_INPUT_DIR}/runICBC_setN_rt_FR.ksh ./.
   cp ${PERT_CHEM_INPUT_DIR}/random.py ./.
   cp ${PERT_CHEM_INPUT_DIR}/run_mozbc_rt_CR.csh ./.
   cp ${PERT_CHEM_INPUT_DIR}/run_mozbc_rt_FR.csh ./.
   cp ${PERT_CHEM_INPUT_DIR}/mozbc-dart/mozbc ./.
   cp ${PERT_CHEM_INPUT_DIR}/set0 ./.
   cp ${PERT_CHEM_INPUT_DIR}/set00 ./.
#
# SELECT MOZART DATA FILE
  if [[ ${MM} -eq 06 ]]; then export MOZART_DATA=h0001.nc; fi
  if [[ ${MM} -eq 07 ]]; then export MOZART_DATA=h0002.nc; fi
#
# CREATE INPUT FILES COARSE DOMAIN
   rm -rf mozbc.both.inp
   cat << EOF > mozbc.both.inp
&control
do_bc     = .true.
do_ic     = .true.
domain    = 1
dir_wrf  = '${RUN_DIR}/${DATE}/wrfchem_chem_icbc/'
dir_moz = '${PERT_CHEM_INPUT_DIR}/MOZART/'
fn_moz  = '${MOZART_DATA}'
def_missing_var = .true.
met_file_prefix  = 'met_em'
met_file_suffix  = '.nc'
met_file_separator= '.'
spc_map = 'o3 -> O3', 'no -> NO',
          'no2 -> NO2', 'no3 -> NO3', 'nh3 -> NH3', 'hno3 -> HNO3', 'hno4 -> HO2NO2',
          'n2o5 -> N2O5','ho2 -> HO2', 'h2o2 -> H2O2',
          'ch4 -> CH4', 'co -> CO', 'ch3o2 -> CH3O2', 'ch3ooh -> CH3OOH',
          'hcho -> CH2O', 'ch3oh -> CH3OH', 'c2h4 -> C2H4',
          'ald -> CH3CHO', 'ch3cooh -> CH3COOH', 'acet -> CH3COCH3', 'mgly -> CH3COCHO',
          'pan -> PAN', 'mpan -> MPAN','macr -> MACR',
          'mvk -> MVK', 'c2h6 -> C2H6', 'c3h6 -> C3H6', 'c3h8 -> C3H8',
          'c2h5oh -> C2H5OH','c10h16 -> C10H16',
          'onit -> ONIT', 'onitr -> ONITR', 'isopr -> ISOP',
          'isopn -> ISOPNO3', 'acetol -> HYAC', 'glyald -> GLYALD',
          'hydrald -> HYDRALD', 'mek -> MEK',
          'bigene -> BIGENE', 'open -> BIGALD', 'bigalk -> BIGALK',
          'tol -> TOLUENE',
          'cres -> CRESOL', 'dms -> DMS', 'so2 -> SO2', 'sulf -> SO4',
          'BC1 -> .4143*CB1;1.e9', 'BC2 -> .4143*CB2;1.e9',
          'OC1 -> .4143*OC1;1.e9', 'OC2 -> .4143*OC2;1.e9',
          'SEAS_1 -> 2.*SA1;1.e9', 'SEAS_2 -> 2.*SA2;1.e9',
          'SEAS_3 -> 2.*SA3;1.e9', 'SEAS_4 -> 2.*SA4;1.e9'
          'DUST_1 -> 1.1738*[DUST1];1.e9', 'DUST_2 -> .939*[DUST2];1.e9',
          'DUST_3 -> .2348*[DUST2]+.939*[DUST3];1.e9',
          'DUST_4 -> .2348*[DUST3]+.5869*[DUST4];1.e9', 'DUST_5 -> .5869*[DUST4];1.e9'
/
EOF
   rm -rf mozbc.ic.inp
   cat << EOF > mozbc.ic.inp
&control
do_bc     = .false.
do_ic     = .true.
domain    = 1
dir_wrf  = '${RUN_DIR}/${DATE}/wrfchem_chem_icbc/'
dir_moz = '${PERT_CHEM_INPUT_DIR}/MOZART/'
fn_moz  = '${MOZART_DATA}'
def_missing_var = .true.
met_file_prefix  = 'met_em'
met_file_suffix  = '.nc'
met_file_separator= '.'
EOF
   rm -rf mozbc.bc.inp
   cat << EOF > mozbc.bc.inp
&control
do_bc     = .true.
do_ic     = .false.
domain    = 1
dir_wrf  = '${RUN_DIR}/${DATE}/wrfchem_chem_icbc/'
dir_moz = '${PERT_CHEM_INPUT_DIR}/MOZART/'
fn_moz  = '${MOZART_DATA}'
def_missing_var = .true.
met_file_prefix  = 'met_em'
met_file_suffix  = '.nc'
met_file_separator= '.'
EOF
   cp ${METGRID_DIR}/met_em.d${CR_DOMAIN}.*:00:00.nc ./.
   ./random.py ${MOZ_SPREAD} ${NUM_MEMBERS} ${PERT_CHEM_INPUT_DIR} ${RUN_DIR}/${DATE}/wrfchem_chem_icbc
   ./runICBC_parent_rt_CR.ksh
   ./runICBC_setN_rt_CR.ksh
#
# CREATE INPUT FILES FINE DOMAIN
   rm -rf mozbc.both.inp
   cat << EOF > mozbc.both.inp
&control
do_bc     = .true.
do_ic     = .true.
domain    = 2
dir_wrf  = '${RUN_DIR}/${DATE}/wrfchem_chem_icbc/'
dir_moz = '${PERT_CHEM_INPUT_DIR}/MOZART/'
fn_moz  = '${MOZART_DATA}'
def_missing_var = .true.
met_file_prefix  = 'met_em'
met_file_suffix  = '.nc'
met_file_separator= '.'
spc_map = 'o3 -> O3', 'no -> NO',
          'no2 -> NO2', 'no3 -> NO3', 'nh3 -> NH3', 'hno3 -> HNO3', 'hno4 -> HO2NO2',
          'n2o5 -> N2O5','ho2 -> HO2', 'h2o2 -> H2O2',
          'ch4 -> CH4', 'co -> CO', 'ch3o2 -> CH3O2', 'ch3ooh -> CH3OOH',
          'hcho -> CH2O', 'ch3oh -> CH3OH', 'c2h4 -> C2H4',
          'ald -> CH3CHO', 'ch3cooh -> CH3COOH', 'acet -> CH3COCH3', 'mgly -> CH3COCHO',
          'pan -> PAN', 'mpan -> MPAN','macr -> MACR',
          'mvk -> MVK', 'c2h6 -> C2H6', 'c3h6 -> C3H6', 'c3h8 -> C3H8',
          'c2h5oh -> C2H5OH','c10h16 -> C10H16',
          'onit -> ONIT', 'onitr -> ONITR', 'isopr -> ISOP',
          'isopn -> ISOPNO3', 'acetol -> HYAC', 'glyald -> GLYALD',
          'hydrald -> HYDRALD', 'mek -> MEK',
          'bigene -> BIGENE', 'open -> BIGALD', 'bigalk -> BIGALK',
          'tol -> TOLUENE',
          'cres -> CRESOL', 'dms -> DMS', 'so2 -> SO2', 'sulf -> SO4',
          'BC1 -> .4143*CB1;1.e9', 'BC2 -> .4143*CB2;1.e9',
          'OC1 -> .4143*OC1;1.e9', 'OC2 -> .4143*OC2;1.e9',
          'SEAS_1 -> 2.*SA1;1.e9', 'SEAS_2 -> 2.*SA2;1.e9',
          'SEAS_3 -> 2.*SA3;1.e9', 'SEAS_4 -> 2.*SA4;1.e9'
          'DUST_1 -> 1.1738*[DUST1];1.e9', 'DUST_2 -> .939*[DUST2];1.e9',
          'DUST_3 -> .2348*[DUST2]+.939*[DUST3];1.e9',
          'DUST_4 -> .2348*[DUST3]+.5869*[DUST4];1.e9', 'DUST_5 -> .5869*[DUST4];1.e9'
/
EOF
   rm -rf mozbc.ic.inp
   cat << EOF > mozbc.ic.inp
&control
do_bc     = .false.
do_ic     = .true.
domain    = 2
dir_wrf  = '${RUN_DIR}/${DATE}/wrfchem_chem_icbc/'
dir_moz = '${PERT_CHEM_INPUT_DIR}/MOZART/'
fn_moz  = '${MOZART_DATA}'
def_missing_var = .true.
met_file_prefix  = 'met_em'
met_file_suffix  = '.nc'
met_file_separator= '.'
EOF
   rm -rf mozbc.bc.inp
   cat << EOF > mozbc.bc.inp
&control
do_bc     = .true.
do_ic     = .false.
domain    = 2
dir_wrf  = '${RUN_DIR}/${DATE}/wrfchem_chem_icbc/'
dir_moz = '${PERT_CHEM_INPUT_DIR}/MOZART/'
fn_moz  = '${MOZART_DATA}'
def_missing_var = .true.
met_file_prefix  = 'met_em'
met_file_suffix  = '.nc'
met_file_separator= '.'
EOF
   cp ${METGRID_DIR}/met_em.d${FR_DOMAIN}.*:00:00.nc ./.
   ./random.py ${MOZ_SPREAD} ${NUM_MEMBERS} ${PERT_CHEM_INPUT_DIR} ${RUN_DIR}/${DATE}/wrfchem_chem_icbc
   ./runICBC_parent_rt_FR.ksh
   ./runICBC_setN_rt_FR.ksh
#
# COMBINE WRFCHEM WITH WRF CR PARENT FILES
   export WRFINPEN=wrfinput_d${CR_DOMAIN}_${YYYY}-${MM}-${DD}_${HH}:00:00
   export WRFBDYEN=wrfbdy_d${CR_DOMAIN}_${YYYY}-${MM}-${DD}_${HH}:00:00
   ncks -A ${REAL_DIR}/${WRFINPEN} ${WRFINPEN}
   ncks -A ${REAL_DIR}/${WRFBDYEN} ${WRFBDYEN}
#
# COMBINE WRFCHEM WITH WRF FR DOMAIN PARENT FILES
   export WRFINPEN=wrfinput_d${FR_DOMAIN}_${YYYY}-${MM}-${DD}_${HH}:00:00
   ncks -A ${REAL_DIR}/${WRFINPEN} ${WRFINPEN}
#
# LOOP THROUGH ALL MEMBERS IN THE ENSEMBLE
   let MEM=1
   while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
      export CMEM=e${MEM}
      if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
      if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
#
# COMBINE WRFCHEM WITH WRF CR DOMAIN
      export WRFINPEN=wrfinput_d${CR_DOMAIN}_${YYYY}-${MM}-${DD}_${HH}:00:00.${CMEM}
      export WRFBDYEN=wrfbdy_d${CR_DOMAIN}_${YYYY}-${MM}-${DD}_${HH}:00:00.${CMEM}
      ncks -A ${WRFCHEM_MET_IC_DIR}/${WRFINPEN} ${WRFINPEN}
      ncks -A ${WRFCHEM_MET_BC_DIR}/${WRFBDYEN} ${WRFBDYEN}
#
# COMBINE WRFCHEM WITH WRF FR DOMAIN
      export WRFINPEN=wrfinput_d${FR_DOMAIN}_${YYYY}-${MM}-${DD}_${HH}:00:00.${CMEM}
      ncks -A ${WRFCHEM_MET_IC_DIR}/${WRFINPEN} ${WRFINPEN}
      let MEM=MEM+1
   done
fi
#
#########################################################################
#
# RUN WRFCHEM PERTURB EMISSIONS
#
#########################################################################
#
if ${RUN_PERT_WRFCHEM_CHEM_EMISS}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_chem_emiss ]]; then
      mkdir ${RUN_DIR}/${DATE}/wrfchem_chem_emiss
      cd ${RUN_DIR}/${DATE}/wrfchem_chem_emiss
   else
      cd ${RUN_DIR}/${DATE}/wrfchem_chem_emiss
   fi
#
# COPY PERTURBATION CODE
   if [[ -e perturb_chem_emiss_RT_CR.exe ]]; then rm -rf perturb_chem_emiss_RT_CR.exe; fi
   if [[ -e perturb_chem_emiss_RT_FR.exe ]]; then rm -rf perturb_chem_emiss_RT_FR.exe; fi
   cp ${PERT_CHEM_EMISS_DIR}/perturb_chem_emiss_RT_CR.exe ./.
   cp ${PERT_CHEM_EMISS_DIR}/perturb_chem_emiss_RT_FR.exe ./.
#
   export L_DATE=${DATE}
   export LE_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${FCST_PERIOD} 2>/dev/null)
   while [[ ${L_DATE} -le ${LE_DATE} ]] ; do
      export L_YYYY=$(echo $L_DATE | cut -c1-4)
      export L_MM=$(echo $L_DATE | cut -c5-6)
      export L_DD=$(echo $L_DATE | cut -c7-8)
      export L_HH=$(echo $L_DATE | cut -c9-10)
#
# LOOP THROUGH ALL MEMBERS IN THE ENSEMBLE
      let MEM=1
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export NL_ENS_MEMBER=${MEM}
         export NL_PERT_CHEM=true
         export NL_PERT_FIRE=true
         export NL_PERT_BIO=false
         if [[ ${L_HH} -eq 00 || ${L_HH} -eq 06 || ${L_HH} -eq 12 || ${L_HH} -eq 18 ]]; then
            export NL_PERT_BIO=true
         fi
         export CMEM=e${MEM}
         if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
#
# GET CR EMISSON FILES FOR THIS MEMBER
         export WRFCHEMI=wrfchemi_d${CR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
         export WRFFIRECHEMI=wrffirechemi_d${CR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
         export WRFBIOCHEMI=wrfbiochemi_d${CR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
#
         if [[ ${NL_PERT_CHEM} == true ]]; then
            cp ${WRFCHEM_CHEMI_DIR}/${WRFCHEMI} ${WRFCHEMI}.${CMEM}
         fi
         if [[ ${NL_PERT_FIRE} == true ]]; then
            cp ${WRFCHEM_FIRE_DIR}/${WRFFIRECHEMI} ${WRFFIRECHEMI}.${CMEM}
         fi
         if [[ ${NL_PERT_BIO} == true ]]; then
            cp ${WRFCHEM_BIO_DIR}/${WRFBIOCHEMI} ${WRFBIOCHEMI}.${CMEM}
         fi
#
# CREATE NAMELIST
         rm -rf perturb_chem_emiss_nml.nl
         cat << EOF > perturb_chem_emiss_nml.nl
&perturb_chem_emiss_nml
idate=${L_DATE},
ens_member=${NL_ENS_MEMBER},
tr_mean=${NL_MEAN},
tr_stdev=${NL_SPREAD},
wrfchemi='${WRFCHEMI}.${CMEM}',
wrffirechemi='${WRFFIRECHEMI}.${CMEM}',
wrfbiochemi='${WRFBIOCHEMI}.${CMEM}',
pert_chem=${NL_PERT_CHEM},
pert_fire=${NL_PERT_FIRE},
pert_bio=${NL_PERT_BIO}
/
EOF
#
# RUN PERTURBATION CODE    
         ./perturb_chem_emiss_RT_CR.exe
#
# GET FR EMISSON FILES FOR THIS MEMBER
         export WRFCHEMI=wrfchemi_d${FR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
         export WRFFIRECHEMI=wrffirechemi_d${FR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
         export WRFBIOCHEMI=wrfbiochemi_d${FR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
#
         if [[ ${NL_PERT_CHEM} == true ]]; then
            cp ${WRFCHEM_CHEMI_DIR}/${WRFCHEMI} ${WRFCHEMI}.${CMEM}
         fi
         if [[ ${NL_PERT_FIRE} == true ]]; then
            cp ${WRFCHEM_FIRE_DIR}/${WRFFIRECHEMI} ${WRFFIRECHEMI}.${CMEM}
         fi
         if [[ ${NL_PERT_BIO} == true ]]; then
            cp ${WRFCHEM_BIO_DIR}/${WRFBIOCHEMI} ${WRFBIOCHEMI}.${CMEM}
         fi
#
# CREATE NAMELIST
         rm -rf perturb_chem_emiss_nml.nl
         cat << EOF > perturb_chem_emiss_nml.nl
&perturb_chem_emiss_nml
idate=${L_DATE},
ens_member=${NL_ENS_MEMBER},
tr_mean=${NL_MEAN},
tr_stdev=${NL_SPREAD},
wrfchemi='${WRFCHEMI}.${CMEM}',
wrffirechemi='${WRFFIRECHEMI}.${CMEM}',
wrfbiochemi='${WRFBIOCHEMI}.${CMEM}',
pert_chem=${NL_PERT_CHEM},
pert_fire=${NL_PERT_FIRE},
pert_bio=${NL_PERT_BIO}
/
EOF
#
# RUN PERTURBATION CODE    
         ./perturb_chem_emiss_RT_FR.exe
#
# GO TO NEXT MEMBER
         let MEM=${MEM}+1
      done
#
# ADVANCE TIME
      export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 1 2>/dev/null)
   done
fi
#
#########################################################################
#
# RUN MOPITT CO OBSERVATIONS
#
#########################################################################
#
if ${RUN_MOPITT_CO_OBS}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/mopitt_co_obs ]]; then
      mkdir ${RUN_DIR}/${DATE}/mopitt_co_obs
      cd ${RUN_DIR}/${DATE}/mopitt_co_obs
   else
      cd ${RUN_DIR}/${DATE}/mopitt_co_obs
   fi
#
# SET MOPITT PARAMETERS
   export MOPITT_FILE_PRE=MOP02J-
   export MOPITT_FILE_EXT=-L2V10.1.3.beta.hdf   
#
#  SET OBS WINDOW
   export BIN_BEG=${ASIM_MN_HH}
   export BIN_END=${ASIM_MX_HH}
   export FLG=0
   if [[ ${BIN_END} -eq 3 ]]; then
      export FLG=1
      export BIN_END=24
   fi
#
# SET MOPITT INPUT DATA DIR
   export MOPITT_INPUT_DIR=/glade/p/acd/mizzi/for_arthur/MOPITT/hdf_files
   export MOPITT_EXE_DIR=/glade/p/acd/mizzi/for_arthur/MOPITT/idl
   export MOPITT_IDL_DIR=/glade/p/acd/mizzi/for_arthur/MOPITT/idl
   export MOP_INFILE="'"${MOPITT_INPUT_DIR}/${MOPITT_FILE_PRE}${YYYY}${MM}${DD}${MOPITT_FILE_EXT}"'"
   export MOP_OUTFILE="'"MOPITT_CO_${D_DATE}'.dat'"'"
#
# COPY EXECUTABLE
   export FILE=mopitt_extract_no_transform_RT.pro
   cp ${MOPITT_IDL_DIR}/${FILE} ./.
   rm -rf job.ksh
   rm -rf idl_*.err
   rm -rf idl_*.out
   touch job.ksh
   RANDOM=$$
   export JOBRND=idl_${RANDOM}
   cat <<EOFF >job.ksh
#!/bin/csh -fx
#BSUB -P ${PROJ_NUMBER}
#BSUB -n 1                                  # number of total (MPI) tasks
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W 00:10                              # wallclock time (minutes)
#BSUB -q geyser
#
idl << EOF
.compile mopitt_extract_no_transform_RT.pro
mopitt_extract_no_transform_RT, ${MOP_INFILE}, ${MOP_OUTFILE}, ${BIN_BEG}, ${BIN_END}
exit
EOF
EOFF
   bsub -K < job.ksh
#
# GET ADDITIONAL DATA FOR DAY-TO-DAY CROSSOVER
   if [[ ${FLG} -eq 1 ]];  then
      export FLG=0
      export BIN_BEG=0
      export BIN_END=3
      export MOP_INFILE="'"${MOPITT_INPUT_DIR}/${MOPITT_FILE_PRE}${ASIM_MX_YYYY}${ASIM_MX_MM}${ASIM_MX_DD}${MOPITT_FILE_EXT}"'"
      rm -rf job.ksh
      rm -rf idl_*.err
      rm -rf idl_*.out
      touch job.ksh
      RANDOM=$$
      export JOBRND=idl_${RANDOM}
      cat <<EOFF >job.ksh
#!/bin/csh -fx
#BSUB -P ${PROJ_NUMBER}
#BSUB -n 1                                  # number of total (MPI) tasks
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W 00:10                              # wallclock time (minutes)
#BSUB -q geyser
#
idl << EOF
.compile mopitt_extract_no_transform_RT.pro
mopitt_extract_no_transform_RT, ${MOP_INFILE}, ${MOP_OUTFILE}, ${BIN_BEG}, ${BIN_END}
exit
EOF
EOFF
      bsub -K < job.ksh
   fi   
#
# SET NAMELIST TO CONVERT MOPITT ASCII TO OBS_SEQ 
   export NL_YEAR=${D_YYYY}
   export NL_MONTH=${D_MM}
   export NL_DAY=${D_DD}
   export NL_HOUR=${D_HH}
   if [[ ${D_HH} -eq 24 ]]; then
      export NL_BIN_BEG=21.01
      export NL_BIN_END=3.00
   elif [[ ${D_HH} -eq 6 ]]; then
      export NL_BIN_BEG=3.01
      export NL_BIN_END=9.00
   elif [[ ${D_HH} -eq 12 ]]; then
      export NL_BIN_BEG=9.01
      export NL_BIN_END=15.00
   elif [[ ${D_HH} -eq 18 ]]; then
      export NL_BIN_BEG=15.01
      export NL_BIN_END=21.00
   fi
   cp MOPITT_CO_${D_DATE}.dat ${D_DATE}.dat
   export NL_FILEDIR="'"./"'" 
   export NL_FILENAME=${D_DATE}.dat
#
# USE MOPITT DATA 
   rm -rf input.nml
   ${HYBRID_SCRIPTS_DIR}/da_create_dart_mopitt_input_nml.ksh
#
# GET EXECUTABLE
   cp ${DART_DIR}/observations/MOPITT_CO/work/mopitt_ascii_to_obs ./.
   ./mopitt_ascii_to_obs
#
# COPY OUTPUT TO ARCHIVE LOCATION
   export MOPITT_FILE=mopitt_obs_seq${D_DATE}
   cp ${MOPITT_FILE} obs_seq_mopitt_co_${DATE}.out
fi
#
#########################################################################
#
# RUN IASI CO OBSERVATIONS
#
#########################################################################
#
if ${RUN_IASI_CO_OBS}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/iasi_co_obs ]]; then
      mkdir ${RUN_DIR}/${DATE}/iasi_co_obs
      cd ${RUN_DIR}/${DATE}/iasi_co_obs
   else
      cd ${RUN_DIR}/${DATE}/iasi_co_obs
   fi
#
# NO CODE YET
fi
#
#########################################################################
#
# RUN IASI O3 OBSERVATIONS
#
#########################################################################
#
if ${RUN_IASI_O3_OBS}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/iasi_o3_obs ]]; then
      mkdir ${RUN_DIR}/${DATE}/iasio3_obs
      cd ${RUN_DIR}/${DATE}/iasio3_obs
   else
      cd ${RUN_DIR}/${DATE}/iasi_o3_obs
   fi
#
# INSERT CODE HERE (NO CODE YET)
fi
#
#########################################################################
#
# RUN PREPBUFR MET OBSERVATIONS
#
#########################################################################
#
# APM: This block needs to be revised so we can convert a single prepbufr
#      file in real time we can use only the obs that are on the current
#      prepbufr file.
#
if ${RUN_MET_OBS}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/prepbufr_met_obs ]]; then
      mkdir ${RUN_DIR}/${DATE}/prepbufr_met_obs
      cd ${RUN_DIR}/${DATE}/prepbufr_met_obs
   else
      cd ${RUN_DIR}/${DATE}/prepbufr_met_obs
   fi
#
# GET PREPBUFR FILES
# (THIS ASSUMES +/- 3HR OBS WINDOW) 
   if [[ ${D_HH} -eq 24 ]]; then
      export LL_HH=18
   elif [[ ${D_HH} -eq 6 ]]; then
      export LL_HH=00
   elif [[ ${D_HH} -eq 12 ]]; then
      export LL_HH=06
   elif [[ ${D_HH} -eq 18 ]]; then
      export LL_HH=12
   fi
#           
   export L_DATE=${D_YYYY}${D_MM}${D_DD}${LL_HH}
   export E_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} +6 2>/dev/null)
   while [[ ${L_DATE} -le ${E_DATE} ]]; do
      export L_YYYY=$(echo $L_DATE | cut -c1-4)
      export L_YY=$(echo $L_DATE | cut -c3-4)
      export L_MM=$(echo $L_DATE | cut -c5-6)
      export L_DD=$(echo $L_DATE | cut -c7-8)
      export L_HH=$(echo $L_DATE | cut -c9-10)
      cp ${PREPBUFR_DATA_DIR}/${L_YYYY}${L_MM}${L_DD}${L_HH}/prepbufr.gdas.${L_YYYY}${L_MM}${L_DD}${L_HH}.wo40.be prepqm${L_YY}${L_MM}${L_DD}${L_HH}
      export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} +6 2>/dev/null)
   done
#
# GET DART input.nml
   rm -rf input.nml
   cp ${DART_DIR}/observations/NCEP/prep_bufr/work/input.nml ./.
#
# RUN_PREPBUFR TO ASCII CONVERTER
   ${DART_DIR}/observations/NCEP/prep_bufr/work/prepbufr.csh_RT ${D_YYYY} ${DD_MM} ${DD_DD} ${DD_DD} ${DART_DIR}/observations/NCEP/prep_bufr/exe > index.file
#
# RUN ASCII TO OBS_SEQ CONVERTER
   ${HYBRID_SCRIPTS_DIR}/da_create_dart_ncep_ascii_to_obs_input_nml_RT.ksh
   ${DART_DIR}/observations/NCEP/ascii_to_obs/work/create_real_obs
#
   mv obs_seq${D_DATE} obs_seq_${DATE}.out
fi
#
#########################################################################
#
# RUN COMBINE OBSERVATIONS
#
#########################################################################
#
if ${RUN_COMBINE_OBS}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/combine_obs ]]; then
      mkdir ${RUN_DIR}/${DATE}/combine_obs
      cd ${RUN_DIR}/${DATE}/combine_obs
   else
      cd ${RUN_DIR}/${DATE}/combine_obs
   fi
#
# GET EXECUTABLES
   cp ${DART_DIR}/models/wrf_chem/work/obs_sequence_tool ./.
   cp ${DART_DIR}/models/wrf_chem/work/input.nml ./.
#
# GET OBS_SEQ FILES TO COMBINE
# MET OBS
   if [[ -e ${PREPBUFR_MET_OBS_DIR}/obs_seq_${DATE}.out && ${RUN_MET_OBS} ]]; then 
      cp ${PREPBUFR_MET_OBS_DIR}/obs_seq_${DATE}.out ./obs_seq_MET_${DATE}.out
   elif ${RUN_MET_OBS}; then
      echo APM: ERROR in COMBINE_OBS_DIR obs_seq_${DATE}.out does not exist    
      exit
   fi
#
# MOPITT CO
   if [[ -e ${MOPITT_CO_OBS_DIR}/obs_seq_mopitt_co_${DATE}.out && ${RUN_MOPITT_CO_OBS} ]]; then 
      cp ${MOPITT_CO_OBS_DIR}/obs_seq_mopitt_co_${DATE}.out ./obs_seq_MOP_CO_${DATE}.out
   elif ${RUN_MOPITT_CO_OBS}; then
      echo APM: ERROR in COMBINE_OBS_DIR obs_seq_mopitt_co_${DATE}.out does not exist    
      exit
   fi
#
# IASI CO
   if [[ -e ${IASI_CO_OBS_DIR}/obs_seq_iasi_co_${DATE}.out && ${RUN_IASI_CO_OBS} ]]; then 
      cp ${IASI_CO_OBS_DIR}/obs_seq_iasi_co_${DATE}.out ./obs_seq_IAS_CO_${DATE}.out
   elif ${RUN_IASI_CO_OBS}; then
      echo APM: ERROR in COMBINE_OBS_DIR obs_seq_iasi_co_${DATE}.out does not exist    
      exit
   fi
#
# IASI O3
   if [[ -e ${IASI_O3_OBS_DIR}/obs_seq_iasi_o3_${DATE}.out && ${RUN_IASI_O3_OBS} ]]; then 
      cp ${IASI_O3_OBS_DIR}/obs_seq_iasi_o3_${DATE}.out ./obs_seq_IAS_O3_${DATE}.out   
   elif ${RUN_IASI_O3_OBS}; then
      echo APM: ERROR in COMBINE_OBS_DIR obs_seq_iasi_o3_${DATE}.out does not exist    
      exit
   fi
#
# SETUP OBS_SEQUENCE_TOOL INPUT.NML
   export RUN_MET_OBS=true
   export RUN_MOPITT_CO_OBS=true
   export NUM_FILES=0
   if ${RUN_MET_OBS}; then let NUM_FILES=${NUM_FILES}+1; fi
   if ${RUN_MOPITT_CO_OBS}; then let NUM_FILES=${NUM_FILES}+1; fi
   if ${RUN_IASI_CO_OBS}; then let NUM_FILES=${NUM_FILES}+1; fi
   if ${RUN_IASI_O3_OBS}; then let NUM_FILES=${NUM_FILES}+1; fi
   export NL_NUM_INPUT_FILES=${NUM_FILES}
#
# APM: How to handle the following parameter definition?
   export NL_FILENAME_SEQ="'obs_seq_MET_${DATE}.out','obs_seq_MOP_CO_${DATE}.out'"
   export NL_FILENAME_OUT="'obs_seq.proc'"
   export NL_FIRST_OBS_DAYS=${ASIM_MIN_DAY_GREG}
   export NL_FIRST_OBS_SECONDS=${ASIM_MIN_SEC_GREG}
   export NL_LAST_OBS_DAYS=${ASIM_MAX_DAY_GREG}
   export NL_LAST_OBS_SECONDS=${ASIM_MAX_SEC_GREG}
   export NL_SYNONYMOUS_COPY_LIST="'NCEP BUFR observation','MOPITT CO observation','IASI CO observation','IASI O3 observation'"
   export NL_SYNONYMOUS_QC_LIST="'NCEP QC index','MOPITT CO QC index','IASI CO QC index','IASI O3 QC index'"
   rm -rf input.nml
   ${HYBRID_SCRIPTS_DIR}/da_create_dart_input_nml.ksh       
#
# Make obs_def_apm_nml for apm_scale to adcdjust observation error variance
   rm -rf obs_def_apm.nml
   cat <<EOF > obs_def_apm.nml
&obs_def_apm_nml
apm_scale=${NL_APM_SCALE}
apm_scale_sw=${NL_APM_SCALE_SW}
/
EOF
#
   ./obs_sequence_tool
   mv obs_seq.proc obs_seq_comb_${DATE}.out
fi
#
#########################################################################
#
# RUN PREPROCESS OBSERVATIONS
#
#########################################################################
#
if ${RUN_PREPROCESS_OBS}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/preprocess_obs ]]; then
      mkdir ${RUN_DIR}/${DATE}/preprocess_obs
      cd ${RUN_DIR}/${DATE}/preprocess_obs
   else
      cd ${RUN_DIR}/${DATE}/preprocess_obs
   fi
#
# GET WRFINPUT TEMPLATE
   cp ${WRFCHEM_MET_IC_DIR}/wrfinput_d${CR_DOMAIN}_${FILE_DATE}.e001 wrfinput_d${CR_DOMAIN}
#
# GET DART UTILITIES
   cp ${DART_DIR}/models/wrf_chem/work/wrf_dart_obs_preprocess ./.
   cp ${DART_DIR}/models/wrf_chem/WRF_DART_utilities/wrf_dart_obs_preprocess.nml ./.
   cp ${DART_DIR}/models/wrf_chem/work/input.nml ./.
#
# GET INPUT DATA
   cp ${COMBINE_OBS_DIR}/obs_seq_comb_${DATE}.out obs_seq.old
#
# MAKE DART OBS ERROR SCALING NAMELIST
   rm -rf obs_def_apm.nml
   cat <<EOF > obs_def_apm.nml
&obs_def_apm_nml
apm_scale=${NL_APM_SCALE}
apm_scale_sw=${NL_APM_SCALE_SW}
/
EOF
#
# CREATE JOB SCRIPT
   if [[ -f job.ksh ]]; then rm -rf job.ksh; fi
   touch job.ksh
   RANDOM=$$
   export JOBRND=pre_$RANDOM
   cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -n 1                                  # number of total (MPI) tas    758 ks
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W 00:10                              # wallclock time (minutes)
#BSUB -q geyser
#
# Run wrf_obs_preprocess
rm -rf pre_*.err
rm -rf pre_*.out
./wrf_dart_obs_preprocess ${DAY_GREG} ${SEC_GREG} > index_preprocess 2>&1 
#
export RC=\$?     
if [[ -f PRE_SUCCESS ]]; then rm -rf PRE_SUCCESS; fi     
if [[ -f PRE_FAILED ]]; then rm -rf PRE_FAILED; fi          
if [[ \$RC = 0 ]]; then
   touch PRE_SUCCESS
else
   touch PRE_FAILED 
   exit
fi
EOF
#
   bsub -K < job.ksh 
#
# SAVE OUTPUT
   mv obs_seq.new obs_seq_comb_filtered_${DATE}.out 
fi
#
#########################################################################
#
# RUN WRF-CHEM INITAL (NO CYCLING-BASED FIRST GUESS FOR DART)
#
#########################################################################
#
if ${RUN_WRFCHEM_INITIAL}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_initial ]]; then
      mkdir -p ${RUN_DIR}/${DATE}/wrfchem_initial
      cd ${RUN_DIR}/${DATE}/wrfchem_initial
   else
      cd ${RUN_DIR}/${DATE}/wrfchem_initial
   fi
#
# Run WRF-Chem for all ensemble members
   export RAN_APM=${RANDOM}
   let MEM=1
   while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
      export CMEM=e${MEM}
      export KMEM=${MEM}
      if [[ ${MEM} -lt 1000 ]]; then export KMEM=0${MEM}; fi
      if [[ ${MEM} -lt 100 ]]; then export KMEM=00${MEM}; export CMEM=e0${MEM}; fi
      if [[ ${MEM} -lt 10 ]]; then export KMEM=000${MEM}; export CMEM=e00${MEM}; fi
      export L_RUN_DIR=run_${CMEM}
      cd ${RUN_DIR}/${DATE}/wrfchem_initial
      if [[ ! -e ${L_RUN_DIR} ]]; then
         mkdir ${L_RUN_DIR}
         cd ${L_RUN_DIR}
      else
         cd ${L_RUN_DIR}
      fi
#
# Get WRF-Chem parameter files
      cp ${WRFCHEM_DIR}/test/em_real/wrf.exe ./.
      cp ${WRFCHEM_DIR}/test/em_real/aerosol.formatted ./.
      cp ${WRFCHEM_DIR}/test/em_real/aerosol_lat.formatted ./.
      cp ${WRFCHEM_DIR}/test/em_real/aerosol_lon.formatted ./.
      cp ${WRFCHEM_DIR}/test/em_real/aerosol_plev.formatted ./.
      cp ${WRFCHEM_DIR}/test/em_real/bulkdens.asc_s_0_03_0_9 ./.
      cp ${WRFCHEM_DIR}/test/em_real/bulkradii.asc_s_0_03_0_9 ./.
      cp ${WRFCHEM_DIR}/test/em_real/CAM_ABS_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/CAM_AEROPT_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/CAMtr_volume_mixing_ratio ./.
      cp ${WRFCHEM_DIR}/test/em_real/CAMtr_volume_mixing_ratio.A1B ./.
      cp ${WRFCHEM_DIR}/test/em_real/CAMtr_volume_mixing_ratio.A2 ./.
      cp ${WRFCHEM_DIR}/test/em_real/CAMtr_volume_mixing_ratio.RCP4.5 ./.
      cp ${WRFCHEM_DIR}/test/em_real/CAMtr_volume_mixing_ratio.RCP6 ./.
      cp ${WRFCHEM_DIR}/test/em_real/capacity.asc ./.
      cp ${WRFCHEM_DIR}/test/em_real/CCN_ACTIVATE.BIN ./.
      cp ${WRFCHEM_DIR}/test/em_real/CLM_ALB_ICE_DFS_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/CLM_ALB_ICE_DRC_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/CLM_ASM_ICE_DFS_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/CLM_ASM_ICE_DRC_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/CLM_DRDSDT0_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/CLM_EXT_ICE_DFS_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/CLM_EXT_ICE_DRC_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/CLM_KAPPA_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/CLM_TAU_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/coeff_p.asc ./.
      cp ${WRFCHEM_DIR}/test/em_real/coeff_q.asc ./.
      cp ${WRFCHEM_DIR}/test/em_real/constants.asc ./.
      cp ${WRFCHEM_DIR}/test/em_real/ETAMPNEW_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/ETAMPNEW_DATA.expanded_rain ./.
      cp ${WRFCHEM_DIR}/test/em_real/GENPARM.TBL ./.
      cp ${WRFCHEM_DIR}/test/em_real/grib2map.tbl ./.
      cp ${WRFCHEM_DIR}/test/em_real/gribmap.txt ./.
      cp ${WRFCHEM_DIR}/test/em_real/kernels.asc_s_0_03_0_9 ./.
      cp ${WRFCHEM_DIR}/test/em_real/kernels_z.asc ./.
      cp ${WRFCHEM_DIR}/test/em_real/LANDUSE.TBL ./.
      cp ${WRFCHEM_DIR}/test/em_real/masses.asc ./.
      cp ${WRFCHEM_DIR}/test/em_real/MPTABLE.TBL ./.
      cp ${WRFCHEM_DIR}/test/em_real/ozone.formatted ./.
      cp ${WRFCHEM_DIR}/test/em_real/ozone_lat.formatted ./.
      cp ${WRFCHEM_DIR}/test/em_real/ozone_plev.formatted ./.
      cp ${WRFCHEM_DIR}/test/em_real/RRTM_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/RRTMG_LW_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/RRTMG_SW_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/SOILPARM.TBL ./.
      cp ${WRFCHEM_DIR}/test/em_real/termvels.asc ./.
      cp ${WRFCHEM_DIR}/test/em_real/tr49t67 ./.
      cp ${WRFCHEM_DIR}/test/em_real/tr49t85 ./.
      cp ${WRFCHEM_DIR}/test/em_real/tr67t85 ./.
      cp ${WRFCHEM_DIR}/test/em_real/URBPARM.TBL ./.
      cp ${WRFCHEM_DIR}/test/em_real/VEGPARM.TBL ./.
      cp ${DART_DIR}/models/wrf_chem/run_scripts/hist_io_flds ./.
#
      cp ${REAL_TIME_DIR}/clim_p_trop.nc ./.
      cp ${REAL_TIME_DIR}/ubvals_b40.20th.track1_1996-2005.nc ./.
      cp ${EXO_COLDENS_DIR}/exo_coldens_d${CR_DOMAIN} ./.
      cp ${SEASONS_WES_DIR}/wrf_season_wes_usgs_d${CR_DOMAIN}.nc ./.
#
# Get WRF-Chem emissions files

      cp ${WRFCHEM_CHEM_EMISS_DIR}/wrfbiochemi_d${CR_DOMAIN}_${START_FILE_DATE}.${CMEM} wrfbiochemi_d${CR_DOMAIN}_${START_FILE_DATE}
      export L_DATE=${START_DATE}
      while [[ ${L_DATE} -le ${END_DATE} ]]; do
         export L_YY=`echo ${L_DATE} | cut -c1-4`
         export L_MM=`echo ${L_DATE} | cut -c5-6`
         export L_DD=`echo ${L_DATE} | cut -c7-8`
         export L_HH=`echo ${L_DATE} | cut -c9-10`
         export L_FILE_DATE=${L_YY}-${L_MM}-${L_DD}_${L_HH}:00:00
         cp ${WRFCHEM_CHEM_EMISS_DIR}/wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE}.${CMEM} wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE}
         cp ${WRFCHEM_CHEM_EMISS_DIR}/wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE}.${CMEM} wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE}
         export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} +1 2>/dev/null)
      done
#
# Get WR-Chem input and bdy files
      cp ${WRFCHEM_CHEM_ICBC_DIR}/wrfinput_d${CR_DOMAIN}_${START_FILE_DATE}.${CMEM} wrfinput_d${CR_DOMAIN}
      cp ${WRFCHEM_CHEM_ICBC_DIR}/wrfbdy_d${CR_DOMAIN}_${START_FILE_DATE}.${CMEM} wrfbdy_d${CR_DOMAIN}
#
# Create WRF-Chem namelist.input
      export NL_MAX_DOM=1
      rm -rf namelist.input
      ${HYBRID_SCRIPTS_DIR}/da_create_wrfchem_namelist_RT.ksh
#      rm -rf namelist.input
#      cp /glade/scratch/mizzi/real_time_wrf_chem/2008072000/wrfchem_initial/run_e001 ./.
#
# Create job script for this member and run it 
      rm -rf job.ksh
      touch job.ksh
      export JOBRND=advm_${RAN_APM}
      cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -x                                    # exclusive use of node (not_shared)
#BSUB -n ${WRFCHEM_NUM_TASKS}                       # number of total (MPI) tasks
#BSUB -R "span[ptile=${WRFCHEM_TASKS_PER_NODE}]"    # mpi tasks per node
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W ${WRFCHEM_TIME_LIMIT}               # wallclock time (minutes)
#BSUB -q ${WRFCHEM_JOB_CLASS}
#
mpirun.lsf ./wrf.exe > index_wrfchem_${KMEM} 2>&1 
export RC=\$?     
rm -rf WRFCHEM_SUCCESS_*; fi     
rm -rf WRFCHEM_FAILED_*; fi          
if [[ \$RC = 0 ]]; then
   touch WRFCHEM_SUCCESS_${RAN_APM}
else
   touch WRFCHEM_FAILED_${RAN_APM} 
   exit
fi
EOF
#
      bsub < job.ksh 
      let MEM=${MEM}+1
   done
#
# Wait for WRFCHEM to complete for each member
   ${HYBRID_SCRIPTS_DIR}/da_run_hold.ksh ${RAN_APM}
fi
#
#########################################################################
#
# RUN DART_FILTER
#
#########################################################################
#
if ${RUN_DART_FILTER}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/dart_filter ]]; then
      mkdir -p ${RUN_DIR}/${DATE}/dart_filter
      cd ${RUN_DIR}/${DATE}/dart_filter
   else
      cd ${RUN_DIR}/${DATE}/dart_filter
   fi
#
# Get DART files
   cp ${DART_DIR}/models/wrf_chem/work/filter      ./.
   cp ${DART_DIR}/system_simulation/final_full_precomputed_tables/final_full.${NUM_MEMBERS} ./.
#
# Get background forecasts
   if [[ ${DATE} -eq ${FIRST_FILTER_DATE} ]]; then
      export BACKGND_FCST_DIR=${WRFCHEM_INITIAL_DIR}
   else
      export BACKGND_FCST_DIR=${WRFCHEM_LAST_CYCLE_CR_DIR}
   fi
#
# Get observations
   if [[ ${PREPROCESS_OBS_DIR}/obs_seq_comb_filtered_${START_DATE}.out ]]; then      
      cp  ${PREPROCESS_OBS_DIR}/obs_seq_comb_filtered_${START_DATE}.out obs_seq.out
   else
      echo APM ERROR: NO DART OBSERVATIONS
      exit
   fi
#
# Run WRF_TO_DART
   export RAN_APM=${RANDOM}
   let MEM=1
   while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
      export CMEM=e${MEM}
      export KMEM=${MEM}
      if [[ ${MEM} -lt 1000 ]]; then export KMEM=0${MEM}; fi
      if [[ ${MEM} -lt 100 ]]; then export KMEM=00${MEM}; export CMEM=e0${MEM}; fi
      if [[ ${MEM} -lt 10 ]]; then export KMEM=000${MEM}; export CMEM=e00${MEM}; fi
#
      cd ${RUN_DIR}/${DATE}/dart_filter
      rm -rf dart_wrk_${CMEM}
      mkdir dart_wrk_${CMEM}
      cd dart_wrk_${CMEM}
#
# &wrf_to_dart_nml
      export NL_DART_RESTART_NAME="'../filter_ic_old.${KMEM}'"
      export NL_PRINT_DATA_RANGES=.false.
      ${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_input.nml.ksh
      cp ${DART_DIR}/models/wrf_chem/work/wrf_to_dart ./.
      cp ${BACKGND_FCST_DIR}/run_${CMEM}/wrfout_d${CR_DOMAIN}_${FILE_DATE} wrfinput_d${CR_DOMAIN}
#
# Create job script 
      rm -rf job.ksh
      touch job.ksh
      RANDOM=$$
      export JOBRND=wr2dt_${RAN_APM}
      cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -n 1                                  # number of total (MPI) tasks
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W 00:05                              # wallclock time (minutes)
#BSUB -q geyser
#
# Run wrf_to_dart
./wrf_to_dart > index_wrf_to_dart 2>&1 
#
export RC=\$?     
rm -rf WRF2DART_SUCCESS_*
rm -rf WRF2DART_FAILED_*
if [[ \$RC = 0 ]]; then
   touch WRF2DART_SUCCESS_${RAN_APM}
else
   touch WRF2DART_FAILED_${RAN_APM} 
   exit
fi
EOF
#
# Submit convert file script for each and wait until job completes
      bsub < job.ksh 
      let MEM=${MEM}+1
   done
#
# Wait for wrf_to_dart to complete for each member
   ${HYBRID_SCRIPTS_DIR}/da_run_hold.ksh ${RAN_APM}
#
   cd ${RUN_DIR}/${DATE}/dart_filter
   ${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_input.nml.ksh
   cp ${BACKGND_FCST_DIR}/run_e001/wrfout_d${CR_DOMAIN}_${FILE_DATE} wrfinput_d${CR_DOMAIN}
#
# Copy "out" inflation files from prior cycle to "in" inflation files for current cycle
   if ${USE_DART_INFL}; then
      if [[ ${DATE} -eq ${FIRST_FILTER_DATE} ]]; then
         export NL_INF_INITIAL_FROM_RESTART_PRIOR=.false.
         export NL_INF_SD_INITIAL_FROM_RESTART_PRIOR=.false.
         export NL_INF_INITIAL_FROM_RESTART_POST=.false.
         export NL_INF_SD_INITIAL_FROM_RESTART_POST=.false.
      else
         export NL_INF_INITIAL_FROM_RESTART_PRIOR=.true.
         export NL_INF_SD_INITIAL_FROM_RESTART_PRIOR=.true.
         export NL_INF_INITIAL_FROM_RESTART_POST=.true.
         export NL_INF_SD_INITIAL_FROM_RESTART_POST=.true.
      fi
      if [[ ${DATE} -ne ${FIRST_FILTER_DATE} ]]; then
         if [[ ${NL_INF_FLAVOR_PRIOR} != 0 ]]; then
            export INF_OUT_FILE_NAME_PRIOR=${RUN_DIR}/${PAST_DATE}/dart_filter/prior_inflate_ic_new
            cp ${INF_OUT_FILE_NAME_PRIOR} prior_inflate_ic_old
         fi
         if [[ ${NL_INF_FLAVOR_POST} != 0 ]]; then
            export INF_OUT_FILE_NAME_POST=${RUN_DIR}/${PAST_DATE}/dart_filter/post_inflate_ic_new
            cp ${NL_INF_OUT_FILE_NAME_POST} post_infalte_ic_old
         fi 
      fi
   fi
#
# Generate input.nml
   set -A temp `echo ${ASIM_MIN_DATE} 0 -g | ${DART_DIR}/models/wrf_chem/work/advance_time`
   (( temp[1]=${temp[1]}+1 ))
   export NL_FIRST_OBS_DAYS=${temp[0]}
   export NL_FIRST_OBS_SECONDS=${temp[1]}
   set -A temp `echo ${ASIM_MAX_DATE} 0 -g | ${DART_DIR}/models/wrf_chem/work/advance_time`
   export NL_LAST_OBS_DAYS=${temp[0]}
   export NL_LAST_OBS_SECONDS=${temp[1]}
   rm -rf input.nml
   ${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_input.nml.ksh
#
# Make filter_apm_nml for special_outlier_threshold
   rm -rf filter_apm.nml
   cat << EOF > filter_apm.nml
&filter_apm_nml
special_outlier_threshold=${NL_SPECIAL_OUTLIER_THRESHOLD}
/
EOF
#
# Make obs_def_apm_nml for apm_scale to adjust observation error variance
  rm -rf obs_def_apm.nml
  cat << EOF > obs_def_apm.nml
&obs_def_apm_nml
apm_scale=${NL_APM_SCALE}
apm_scale_sw=${NL_APM_SCALE_SW}
/
EOF
#
# Run DART_FILTER
# Create job script for this member and run it 
   rm -rf job.ksh
   touch job.ksh
   RANDOM=$$
   export JOBRND=filter_${RANDOM}
   cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -x                                    # exclusive use of node (not_shared)
#BSUB -n ${FILTER_NUM_TASKS}                       # number of total (MPI) tasks
#BSUB -R "span[ptile=${FILTER_TASKS_PER_NODE}]"    # mpi tasks per node
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W ${FILTER_TIME_LIMIT}               # wallclock time (minutes)
#BSUB -q ${FILTER_JOB_CLASS}
#
mpirun.lsf ./filter > index_filter 2>&1 
export RC=\$?     
rm -rf FILTER_SUCCESS     
rm -rf FILTER_FAILED          
if [[ \$RC = 0 ]]; then
   touch FILTER_SUCCESS
else
   touch FILTER_FAILED 
   exit
fi
EOF
#
   bsub -K < job.ksh
#
# Run DART_TO_WRF 
   export RAN_APM=${RANDOM}
   let MEM=1
   while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
      export CMEM=e${MEM}
      export KMEM=${MEM}
      if [[ ${MEM} -lt 1000 ]]; then export KMEM=0${MEM}; fi
      if [[ ${MEM} -lt 100 ]]; then export KMEM=00${MEM}; export CMEM=e0${MEM}; fi
      if [[ ${MEM} -lt 10 ]]; then export KMEM=000${MEM}; export CMEM=e00${MEM}; fi
#
      cd ${RUN_DIR}/${DATE}/dart_filter
      rm -rf dart_wrk_${CMEM}
      mkdir dart_wrk_${CMEM}
      cd dart_wrk_${CMEM}
#
# &dart_to_wrf_nml
      export NL_MODEL_ADVANCE_FILE=.false.
      export NL_DART_RESTART_NAME="'"../filter_ic_new.${KMEM}"'"
      rm -rf cd dinput.nml
      ${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_input.nml.ksh
      cp ${DART_DIR}/models/wrf_chem/work/dart_to_wrf ./.
      cp ${BACKGND_FCST_DIR}/run_${CMEM}/wrfout_d${CR_DOMAIN}_${FILE_DATE} wrfinput_d${CR_DOMAIN}
#
# Create job script 
      rm -rf job.ksh
      touch job.ksh
      RANDOM=$$
      export JOBRND=dt2wf_${RAN_APM}
      cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -n 1                                  # number of total (MPI) tasks
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W 00:05                              # wallclock time (minutes)
#BSUB -q geyser
#
# Run wrf_to_dart
./dart_to_wrf > index_dart_to_wrf 2>&1 
#
export RC=\$?     
rm -rf DART2WRF_SUCCESS_*
rm -rf DART2WRF_FAILED_*
if [[ \$RC = 0 ]]; then
   touch DART2WRF_SUCCESS_{RAN_APM}
else
   touch DART2WRF_FAILED_{RAN_APM}
   exit
fi
EOF
#
# Submit convert file script for each and wait until job completes
      bsub < job.ksh 
      cp wrfinput_d${CR_DOMAIN} ../wrfout_d${CR_DOMAIN}_${FILE_DATE}_filt.${CMEM} 
      let MEM=${MEM}+1
   done
#
# Wait for dart_to_wrf to complete for each member
   ${HYBRID_SCRIPTS_DIR}/da_run_hold.ksh ${RAN_APM}
fi
#
#########################################################################
#
# UPDATE COARSE RESOLUTION BOUNDARY CONDIIONS
#
#########################################################################
#
if ${RUN_UPDATE_BC}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/update_bc ]]; then
      mkdir -p ${RUN_DIR}/${DATE}/update_bc
      cd ${RUN_DIR}/${DATE}/update_bc
   else
      cd ${RUN_DIR}/${DATE}/update_bc
   fi
#
   let MEM=1
   while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
      export CMEM=e${MEM}
      export KMEM=${MEM}
      if [[ ${MEM} -lt 1000 ]]; then export KMEM=0${MEM}; fi
      if [[ ${MEM} -lt 100 ]]; then export KMEM=00${MEM}; export CMEM=e0${MEM}; fi
      if [[ ${MEM} -lt 10 ]]; then export KMEM=000${MEM}; export CMEM=e00${MEM}; fi
#
      export OPS_FORC_FILE=${WRFCHEM_CHEM_ICBC_DIR}/wrfinput_d${CR_DOMAIN}_${FILE_DATE}.${CMEM}
      export BDYCDN_IN=${WRFCHEM_CHEM_ICBC_DIR}/wrfbdy_d${CR_DOMAIN}_${FILE_DATE}.${CMEM}
      cp ${BDYCDN_IN} wrfbdy_d${CR_DOMAIN}_${FILE_DATE}_prior.${CMEM}
      export DA_OUTPUT_FILE=${DART_FILTER_DIR}/wrfout_d${CR_DOMAIN}_${FILE_DATE}_filt.${CMEM} 
      export BDYCDN_OUT=wrfbdy_d${CR_DOMAIN}_${FILE_DATE}_filt.${CMEM}    
      ${HYBRID_SCRIPTS_DIR}/da_run_update_bc.ksh > index_update_bc 2>&1
#
      let MEM=$MEM+1
   done
fi
#
#########################################################################
#
# RUN WRFCHEM_CYCLE_CR
#
#########################################################################
#
if ${RUN_WRFCHEM_CYCLE_CR}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_cycle_cr ]]; then
      mkdir -p ${RUN_DIR}/${DATE}/wrfchem_cycle_cr
      cd ${RUN_DIR}/${DATE}/wrfchem_cycle_cr
   else
      cd ${RUN_DIR}/${DATE}/wrfchem_cycle_cr
   fi
#
# Run WRF-Chem for all ensemble members
   export RAN_APM=${RANDOM}
   let MEM=1
   while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
      export CMEM=e${MEM}
      export KMEM=${MEM}
      if [[ ${MEM} -lt 1000 ]]; then export KMEM=0${MEM}; fi
      if [[ ${MEM} -lt 100 ]]; then export KMEM=00${MEM}; export CMEM=e0${MEM}; fi
      if [[ ${MEM} -lt 10 ]]; then export KMEM=000${MEM}; export CMEM=e00${MEM}; fi
      export L_RUN_DIR=run_${CMEM}
      cd ${RUN_DIR}/${DATE}/wrfchem_cycle_cr
      if [[ ! -e ${L_RUN_DIR} ]]; then
         mkdir ${L_RUN_DIR}
         cd ${L_RUN_DIR}
      else
         cd ${L_RUN_DIR}
      fi
#
# Get WRF-Chem parameter files
      cp ${WRFCHEM_DIR}/test/em_real/wrf.exe ./.
      cp ${WRFCHEM_DIR}/test/em_real/CAM_ABS_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/CAM_AEROPT_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/ETAMPNEW_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/GENPARM.TBL ./.
      cp ${WRFCHEM_DIR}/test/em_real/LANDUSE.TBL ./.
      cp ${WRFCHEM_DIR}/test/em_real/RRTMG_LW_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/RRTMG_SW_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/RRTM_DATA ./.
      cp ${WRFCHEM_DIR}/test/em_real/SOILPARM.TBL ./.
      cp ${WRFCHEM_DIR}/test/em_real/URBPARM.TBL ./.
      cp ${WRFCHEM_DIR}/test/em_real/VEGPARM.TBL ./.
      cp ${DART_DIR}/models/wrf_chem/run_scripts/hist_io_flds ./.
#
      cp ${REAL_TIME_DIR}/clim_p_trop.nc ./.
      cp ${REAL_TIME_DIR}/ubvals_b40.20th.track1_1996-2005.nc ./.
      cp ${EXO_COLDENS_DIR}/exo_coldens_d${CR_DOMAIN} ./.
      cp ${SEASONS_WES_DIR}/wrf_season_wes_usgs_d${CR_DOMAIN}.nc ./.
#
# Get WRF-Chem emissions files

      cp ${WRFCHEM_CHEM_EMISS_DIR}/wrfbiochemi_d${CR_DOMAIN}_${START_FILE_DATE}.${CMEM} wrfbiochemi_d${CR_DOMAIN}_${START_FILE_DATE}
      export L_DATE=${START_DATE}
      while [[ ${L_DATE} -le ${END_DATE} ]]; do
         export L_YY=`echo ${L_DATE} | cut -c1-4`
         export L_MM=`echo ${L_DATE} | cut -c5-6`
         export L_DD=`echo ${L_DATE} | cut -c7-8`
         export L_HH=`echo ${L_DATE} | cut -c9-10`
         export L_FILE_DATE=${L_YY}-${L_MM}-${L_DD}_${L_HH}:00:00
         cp ${WRFCHEM_CHEM_EMISS_DIR}/wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE}.${CMEM} wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE}
         cp ${WRFCHEM_CHEM_EMISS_DIR}/wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE}.${CMEM} wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE}
         export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} +1 2>/dev/null)
      done
#
# Get WR-Chem input and bdy files
      cp ${DART_FILTER_DIR}/wrfout_d${CR_DOMAIN}_${START_FILE_DATE}_filt.${CMEM} wrfinput_d${CR_DOMAIN}
      cp ${UPDATE_BC_DIR}/wrfbdy_d${CR_DOMAIN}_${START_FILE_DATE}_filt.${CMEM} wrfbdy_d${CR_DOMAIN}
#
# Create WRF-Chem namelist.input 
      export NL_MAX_DOM=1
      rm -rf namelist.input
      ${HYBRID_SCRIPTS_DIR}/da_create_wrfchem_namelist_RT.ksh
#
# Create job script for this member and run it 
      rm -rf job.ksh
      touch job.ksh
      export JOBRND=advm_${RAN_APM}
      cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -x                                    # exclusive use of node (not_shared)
#BSUB -n ${WRFCHEM_NUM_TASKS}                       # number of total (MPI) tasks
#BSUB -R "span[ptile=${WRFCHEM_TASKS_PER_NODE}]"    # mpi tasks per node
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W ${WRFCHEM_TIME_LIMIT}               # wallclock time (minutes)
#BSUB -q ${WRFCHEM_JOB_CLASS}
#
mpirun.lsf ./wrf.exe > index_wrfchem_${KMEM} 2>&1 
export RC=\$?     
rm -rf WRFCHEM_SUCCESS_*; fi     
rm -rf WRFCHEM_FAILED_*; fi          
if [[ \$RC = 0 ]]; then
   touch WRFCHEM_SUCCESS_${RAN_APM}     
else
   touch WRFCHEM_FAILED_${RAN_APM} 
   exit
fi
EOF
#
      bsub < job.ksh 
      let MEM=${MEM}+1
   done
#
# Wait for WRFCHEM to complete for each member
   ${HYBRID_SCRIPTS_DIR}/da_run_hold.ksh ${RAN_APM}
fi
#
#########################################################################
#
# CALCULATE ENSEMBLE MEAN
#
#########################################################################
#
#
#########################################################################
#
# FIND ENSEMBLE MEMBER CLOSEST TO ENSEMBLE MEAN
#
#########################################################################
#
   export CLOSE_MEM_ID=e001
#
#########################################################################
#
# INTERPOLATE CLOSEST MEMBER FROM COARS TO FINE GRID
#
#########################################################################
#
#
#########################################################################
#
# RUN WRFCHEM_CYCLE_FR
#
#########################################################################
#
if ${RUN_WRFCHEM_CYCLE_FR}; then
   if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_cycle_fr ]]; then
      mkdir -p ${RUN_DIR}/${DATE}/wrfchem_cycle_fr
      cd ${RUN_DIR}/${DATE}/wrfchem_cycle_fr
   else
      cd ${RUN_DIR}/${DATE}/wrfchem_cycle_fr
   fi
#
# Get WRF-Chem parameter files
   cp ${WRFCHEM_DIR}/test/em_real/wrf.exe ./.
   cp ${WRFCHEM_DIR}/test/em_real/CAM_ABS_DATA ./.
   cp ${WRFCHEM_DIR}/test/em_real/CAM_AEROPT_DATA ./.
   cp ${WRFCHEM_DIR}/test/em_real/ETAMPNEW_DATA ./.
   cp ${WRFCHEM_DIR}/test/em_real/GENPARM.TBL ./.
   cp ${WRFCHEM_DIR}/test/em_real/LANDUSE.TBL ./.
   cp ${WRFCHEM_DIR}/test/em_real/RRTMG_LW_DATA ./.
   cp ${WRFCHEM_DIR}/test/em_real/RRTMG_SW_DATA ./.
   cp ${WRFCHEM_DIR}/test/em_real/RRTM_DATA ./.
   cp ${WRFCHEM_DIR}/test/em_real/SOILPARM.TBL ./.
   cp ${WRFCHEM_DIR}/test/em_real/URBPARM.TBL ./.
   cp ${WRFCHEM_DIR}/test/em_real/VEGPARM.TBL ./.
   cp ${DART_DIR}/models/wrf_chem/run_scripts/hist_io_flds ./.
#
   cp ${REAL_TIME_DIR}/clim_p_trop.nc ./.
   cp ${REAL_TIME_DIR}/ubvals_b40.20th.track1_1996-2005.nc ./.
   cp ${EXO_COLDENS_DIR}/exo_coldens_d${CR_DOMAIN} ./.
   cp ${EXO_COLDENS_DIR}/exo_coldens_d${FR_DOMAIN} ./.
   cp ${SEASONS_WES_DIR}/wrf_season_wes_usgs_d${CR_DOMAIN}.nc ./.
   cp ${SEASONS_WES_DIR}/wrf_season_wes_usgs_d${FR_DOMAIN}.nc ./.
#
# Get WRF-Chem emissions files
#   cp ${WRFCHEM_BIO_DIR}/wrfbiochemi_d${CR_DOMAIN}_${START_FILE_DATE} wrfbiochemi_d${CR_DOMAIN}_${START_FILE_DATE}
#   cp ${WRFCHEM_BIO_DIR}/wrfbiochemi_d${FR_DOMAIN}_${START_FILE_DATE} wrfbiochemi_d${FR_DOMAIN}_${START_FILE_DATE}
#
   cp ${WRFCHEM_CHEM_EMISS_DIR}/wrfbiochemi_d${CR_DOMAIN}_${START_FILE_DATE}.${CLOSE_MEM_ID} wrfbiochemi_d${CR_DOMAIN}_${START_FILE_DATE}
   cp ${WRFCHEM_CHEM_EMISS_DIR}/wrfbiochemi_d${FR_DOMAIN}_${START_FILE_DATE}.${CLOSE_MEM_ID} wrfbiochemi_d${FR_DOMAIN}_${START_FILE_DATE}
#
   export L_DATE=${START_DATE}
   while [[ ${L_DATE} -le ${END_DATE} ]]; do
      export L_YY=`echo ${L_DATE} | cut -c1-4`
      export L_MM=`echo ${L_DATE} | cut -c5-6`
      export L_DD=`echo ${L_DATE} | cut -c7-8`
      export L_HH=`echo ${L_DATE} | cut -c9-10`
      export L_FILE_DATE=${L_YY}-${L_MM}-${L_DD}_${L_HH}:00:00
#      cp ${WRFCHEM_FIRE_DIR}/wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE} wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE}
#      cp ${WRFCHEM_CHEMI_DIR}/wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE} wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE}
#      cp ${WRFCHEM_FIRE_DIR}/wrffirechemi_d${FR_DOMAIN}_${L_FILE_DATE} wrffirechemi_d${FR_DOMAIN}_${L_FILE_DATE}
#      cp ${WRFCHEM_CHEMI_DIR}/wrfchemi_d${FR_DOMAIN}_${L_FILE_DATE} wrfchemi_d${FR_DOMAIN}_${L_FILE_DATE}
#
      cp ${WRFCHEM_CHEM_EMISS_DIR}/wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE}.${CLOSE_MEM_ID} wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE}
      cp ${WRFCHEM_CHEM_EMISS_DIR}/wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE}.${CLOSE_MEM_ID} wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE}
      cp ${WRFCHEM_CHEM_EMISS_DIR}/wrffirechemi_d${FR_DOMAIN}_${L_FILE_DATE}.${CLOSE_MEM_ID} wrffirechemi_d${FR_DOMAIN}_${L_FILE_DATE}
      cp ${WRFCHEM_CHEM_EMISS_DIR}/wrfchemi_d${FR_DOMAIN}_${L_FILE_DATE}.${CLOSE_MEM_ID} wrfchemi_d${FR_DOMAIN}_${L_FILE_DATE}
#
      export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} +1 2>/dev/null)
   done
#
# Get WR-Chem input and bdy files
#   cp ${REAL_DIR}/wrfout_d${CR_DOMAIN}_${START_FILE_DATE}_filt wrfinput_d${CR_DOMAIN}
#   cp ${REAL_DIR}/wrfbdy_d${CR_DOMAIN}_${START_FILE_DATE}_filt wrfbdy_d${CR_DOMAIN}
#   cp ${REAL_DIR}/wrfout_d${FR_DOMAIN}_${START_FILE_DATE}_filt wrfinput_d${FR_DOMAIN}
#   cp ${WRFCHEM_CHEM_ICBC_DIR}/wrfinput_d${CR_DOMAIN}_${START_FILE_DATE} wrfinput_d${CR_DOMAIN}
#   cp ${WRFCHEM_CHEM_ICBC_DIR}/wrfbdy_d${CR_DOMAIN}_${START_FILE_DATE} wrfbdy_d${CR_DOMAIN}
#   cp ${WRFCHEM_CHEM_ICBC_DIR}/wrfinput_d${FR_DOMAIN}_${START_FILE_DATE} wrfinput_d${FR_DOMAIN}
   cp ${DART_FILTER_DIR}/wrfout_d${CR_DOMAIN}_${START_FILE_DATE}_filt.${CLOSE_MEM_ID} wrfinput_d${CR_DOMAIN}
   cp ${UPDATE_BC_DIR}/wrfbdy_d${CR_DOMAIN}_${START_FILE_DATE}_filt.${CLOSE_MEM_ID} wrfbdy_d${CR_DOMAIN}
   cp ${REAL_DIR}/wrfinput_d${FR_DOMAIN}_${START_FILE_DATE} wrfinput_d${FR_DOMAIN}

#
# Create WRF-Chem namelist.input 
   export NL_MAX_DOM=2
   rm -rf namelist.input
   ${HYBRID_SCRIPTS_DIR}/da_create_wrfchem_namelist_nested_RT.ksh
#
# Create job script for this member and run it 
   rm -rf job.ksh
   touch job.ksh
   RANDOM=$$
   export JOBRND=advm_${RANDOM}
   cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -x                                    # exclusive use of node (not_shared)
#BSUB -n ${WRFCHEM_NUM_TASKS}                       # number of total (MPI) tasks
#BSUB -R "span[ptile=${WRFCHEM_TASKS_PER_NODE}]"    # mpi tasks per node
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W ${WRFCHEM_TIME_LIMIT}               # wallclock time (minutes)
#BSUB -q ${WRFCHEM_JOB_CLASS}
#
mpirun.lsf ./wrf.exe > index_wrfchem 2>&1 
export RC=\$?     
rm -rf WRFCHEM_SUCCESS; fi     
rm -rf WRFCHEM_FAILED; fi          
if [[ \$RC = 0 ]]; then
   touch WRFCHEM_SUCCESS
else
   touch WRFCHEM_FAILED 
   exit
fi
EOF
#
   bsub < job.ksh 
#
fi
#
