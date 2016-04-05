#!/bin/ksh -aeux
#########################################################################
#
# Purpose: Set global environment variables for real_time_wrf_chem
#
#########################################################################
#
# CYCLE DATE-TIME:
export DATE=2008072000
export YYYY=$(echo $DATE | cut -c1-4)
export YY=$(echo $DATE | cut -c3-4)
export MM=$(echo $DATE | cut -c5-6)
export DD=$(echo $DATE | cut -c7-8)
export HH=$(echo $DATE | cut -c9-10)
export FILE_DATE=${YYYY}-${MM}-${DD}_${HH}:00:00
#
export BUILD_DIR=/glade/p/work/mizzi/TRUNK/WRFDAv3.4_dmpar/var/da
export DART_DIR=/glade/p/work/mizzi/TRUNK/DART_CHEM
cp ${DART_DIR}/models/wrf_chem/work/advance_time ./.
cp ${DART_DIR}/models/wrf_chem/work/input.nml ./.
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
   export PAST_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} -1 2>/dev/null)
   export PAST_YYYY=$(echo $PAST_DATE | cut -c1-4)
   export PAST_YY=$(echo $PAST_DATE | cut -c3-4)
   export PAST_MM=$(echo $PAST_DATE | cut -c5-6)
   export PAST_DD=$(echo $PAST_DATE | cut -c7-8)
   export PAST_HH=$(echo $PAST_DATE | cut -c9-10)
   export D_YYYY=${PAST_YYYY}
   export D_YY=${PAST_YY}
   export D_MM=${PAST_MM}
   export D_DD=${PAST_DD}
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
set -A GREG_DATA `echo $DATE 0 -g | ./advance_time`
export DAY_GREG=${GREG_DATA[0]}
export SEC_GREG=${GREG_DATA[1]}
export ASIM_WINDOW=3
export ASIM_MIN_DATE=$($BUILD_DIR/da_advance_time.exe $DATE -$ASIM_WINDOW 2>/dev/null)
export ASIM_MAX_DATE=$($BUILD_DIR/da_advance_time.exe $DATE +$ASIM_WINDOW 2>/dev/null)
set -A temp `echo $ASIM_MIN_DATE 0 -g | ./advance_time`
export ASIM_MIN_DAY_GREG=${temp[0]}
export ASIM_MIN_SEC_GREG=${temp[1]}
set -A temp `echo $ASIM_MAX_DATE 0 -g | ./advance_time` 
export ASIM_MAX_DAY_GREG=${temp[0]}
export ASIM_MAX_SEC_GREG=${temp[1]}
#
# SELECT COMPONENT RUN OPTIONS:
export RUN_GEOGRID=false
export RUN_UNGRIB=false
export RUN_METGRID=false
export RUN_REAL=false
export RUN_PERT_WRFCHEM_MET_IC=false
export RUN_PERT_WRFCHEM_MET_BC=false
export RUN_EXO_COLDENS=false
export RUN_SEASON_WES=false
export RUN_WRFCHEM_BIO=false
export RUN_WRFCHEM_FIRE=false
export RUN_WRFCHEM_CHEMI=false
export RUN_PERT_WRFCHEM_CHEM_ICBC=false
export RUN_PERT_WRFCHEM_CHEM_EMISS=false
#
export RUN_MOPITT_CO_OBS=false
export RUN_IASI_CO_OBS=false
export RUN_IASI_O3_OBS=false
export RUN_MET_OBS=false
export RUN_COMBINE_OBS=false
export RUN_PREPROCESS_OBS=true
#
export RUN_WRFCHEM_IN=false
export RUN_DART=true
export RUN_WRFCHEM_CR=false
export RUN_WRFCHEM_FR=false
#
export NL_APM_SCALE=1.
export NL_APM_SCALE_SW=.FALSE.
#
# FORECAST PARAMETERS:
export FCST_PERIOD=24
export CYCLE_PERIOD=6
export NUM_MEMBERS=20
export NNXP=101
export NNYP=41
(( LBC_END=2*${FCST_PERIOD} ))
export LBC_FREQ=3
(( INTERVAL_SECONDS=${LBC_FREQ}*60*60 ))
export LBC_START=0
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
export GEOGRID_JOB_CLASS=small
#
# CODE VERSIONS:
export WPS_VER=WPSv3.4_dmpar
export WPS_GEOG_VER=WPSv3.4_GEOG_DATA
export WRFDA_VER=WRFDAv3.4_dmpar
export WRFDA_TOOLS_VER=WRFDA_TOOLSv3.4
export WRF_VER=WRFv3.4_dmpar
export WRFCHEM_VER=WRFCHEMv3.4_dmpar
export DART_VER=DART_CHEM
#
# ROOT DIRECTORIES:
export SCRATCH_DIR=/glade/scratch/mizzi
export WORK_DIR=/glade/p/work/mizzi
export ACD_DIR=/glade/p/acd/mizzi
#
# DEPENDENT DIRECTORIES:
export TRUNK_DIR=${WORK_DIR}/TRUNK
export WPS_DIR=${TRUNK_DIR}/${WPS_VER}
export WPS_GEOG_DIR=${TRUNK_DIR}/${WPS_GEOG_VER}/geog
export WRFVAR_DIR=${TRUNK_DIR}/${WRFDA_VER}
export DART_DIR=${TRUNK_DIR}/${DART_VER}
export BUILD_DIR=${WRFVAR_DIR}/var/da
export WRF_DIR=${TRUNK_DIR}/${WRF_VER}
export HYBRID_TRUNK_DIR=${WORK_DIR}/HYBRID_TRUNK
export HYBRID_SCRIPTS_DIR=${HYBRID_TRUNK_DIR}/hybrid_scripts
export SCRIPTS_DIR=${TRUNK_DIR}/${WRFDA_TOOLS_VER}/scripts
#
export INPUT_DATA_DIR=${ACD_DIR}/AVE_TEST_DATA
export GRIB_DIR=${INPUT_DATA_DIR}/GFS_DATA_NOMADS
export OBSPROC_DIR=${WRFVAR_DIR}/var/obsproc
export COLDENS_INPUT_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_INITIAL_EMISS/WES-COLDENS
export VTABLE_DIR=${WPS_DIR}/ungrib/Variable_Tables
export OB_PREPBUFR_DIR=${INPUT_DATA_DIR}/obs_MET
export BE_DIR=${WRFVAR_DIR}/var/run
export MET_EM_DIR=${ACD_DIR}/AVE_TEST_DATA/met_em_100km
export PREPBUFR_DATA_DIR=${INPUT_DATA_DIR}/obs_MET
export MEGAN_BIO_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_INITIAL_EMISS/MEGAN-BIO
export MEGAN_INPUT_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_INITIAL_EMISS/MEGAN-DATA
export FIRE_INPUT_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_INITIAL_EMISS/WRF-FIRE
export PERT_CHEM_INPUT_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_PERT_CHEM/ICBC_PERT
export PERT_CHEM_EMISS_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_PERT_CHEM/EMISS_PERT
export RUN_DIR=${SCRATCH_DIR}/real_time_wrf_chem
export GEOGRID_DIR=${RUN_DIR}/geogrid
export REAL_DIR=${RUN_DIR}/${DATE}/real
export WRFCHEM_MET_IC_DIR=${RUN_DIR}/${DATE}/wrfchem_met_ic
export WRFCHEM_MET_BC_DIR=${RUN_DIR}/${DATE}/wrfchem_met_bc
export WRFCHEM_BIO_DIR=${RUN_DIR}/${DATE}/wrfchem_bio
export WRFCHEM_FIRE_DIR=${RUN_DIR}/${DATE}/wrfchem_fire
export WRFCHEM_CHEMI_DIR=${RUN_DIR}/${DATE}/wrfchem_chemi
export PREPBUFR_MET_OBS_DIR=${RUN_DIR}/${DATE}/prepbufr_met_obs
export MOPITT_CO_OBS_DIR=${RUN_DIR}/${DATE}/mopitt_co_obs
export IASI_CO_OBS_DIR=${RUN_DIR}/${DATE}/iasi_co_obs
export IASI_O3_OBS_DIR=${RUN_DIR}/${DATE}/iasi_o3_obs
export COMBINE_OBS_DIR=${RUN_DIR}/${DATE}/combine_obs
export PREPROCESS_OBS_DIR=${RUN_DIR}/${DATE}/preprocess_obs
#
if [[ ! -e ${RUN_DIR} ]]; then mkdir ${RUN_DIR}; fi
cd ${RUN_DIR}
#
# WPS PARAMETERS:
export SINGLE_FILE=false
export HOR_SCALE=1500
export VTABLE_TYPE=GFS
export METGRID_TABLE_TYPE=ARW
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
export NL_MAX_DOM=1
export NL_IO_FORM_GEOGRID=2
export NL_OPT_OUTPUT_FROM_GEOGRID_PATH="'"${GEOGRID_DIR}"'"
export NL_DEBUG_LEVEL=0
#
# WPS GEOGRID NAMELIST:
export NL_PARENT_ID=1
export NL_PARENT_GRID_RATIO=1
export NL_I_PARENT_START=1
export NL_J_PARENT_START=1
export NL_S_WE=1
export NL_E_WE=${NNXP}
export NL_S_SN=1
export NL_E_SN=${NNYP}
export NL_GEOG_DATA_RES="'"30s"'"
export NL_DX=100000
export NL_DY=100000
export NL_MAP_PROJ="'"lambert"'"
export NL_REF_LAT=36.0
export NL_REF_LON=-115.0
export NL_TRUELAT1=10.0
export NL_TRUELAT2=50.0
export NL_STAND_LON=-111.0
export NL_GEOG_DATA_PATH="'"${WPS_GEOG_DIR}"'"
export NL_OPT_GEOGRID_TBL_PATH="'"${WPS_DIR}/geogrid"'"
#
# WPS UNGRIB NAMELIST:
export NL_OUT_FORMAT="'"WPS"'"
#
# WPS METGRID NAMELIST:
export NL_IO_FORM_METGRID=2
export NL_OPT_IGNORE_DOM_CENTER=.false.
export NL_OPT_OUTPUT_FROM_METGRID_PATH="'"${RUN_DIR}/rc"'"
#
# WRF NAMELIST:
# TIME CONTROL NAMELIST:
export NL_RUN_DAYS=0
export NL_RUN_HOURS=${FCST_PERIOD}
export NL_RUN_MINUTES=0
export NL_RUN_SECONDS=0
export NL_START_MINUTE=0
export NL_START_SECOND=0
export NL_END_MINUTE=0
export NL_END_SECOND=0
export NL_INTERVAL_SECONDS=${INTERVAL_SECONDS}
export NL_INPUT_FROM_FILE=".true."
export NL_HISTORY_INTERVAL=720
export NL_FRAMES_PER_OUTFILE=1
export NL_RESTART=".false."
export NL_RESTART_INTERVAL=2
export NL_IO_FORM_HISTORY=2
export NL_IO_FORM_RESTART=2
export NL_IO_FORM_INPUT=2
export NL_IO_FORM_BOUNDARY=2
export NL_WRITE_INPUT=".false."
export NL_INPUTOUT_INTERVAL=180
export NL_INPUT_OUTNAME="'"wrfout_d\<domain\>_\<date\>"'"
export NL_DEBUG_LEVEL=0
#
# DOMAINS NAMELIST:
export NL_TIME_STEP=360
export NL_TIME_STEP_FRACT_NUM=0
export NL_TIME_STEP_FRACT_DEN=1
export NL_MAX_DOM=1
export NL_S_WE=1
export NL_E_WE=${NNXP}
export NL_S_SN=1
export NL_E_SN=${NNYP}
#
# This is number of vertical levels in the wrfinput data
export NL_E_VERT=34
export NL_P_TOP_REQUESTED=1000
export NL_INTERP_TYPE=1
export NL_T_EXTRAP_TYPE=1
#
# This in number of vertical levels in the initial data
export NL_NUM_METGRID_LEVELS=27
export NL_NUM_METGRID_SOIL_LEVELS=4
export NL_DX=100000.0
export NL_DY=100000.0
export NL_GRID_ID=1
export NL_PARENT_ID=0
export NL_I_PARENT_START=0
export NL_J_PARENT_START=0
export NL_PARENT_GRID_RATIO=1
export NL_PARENT_TIME_STEP_RATIO=1
export NL_FEEDBACK=0
export NL_SMOOTH_OPTION=1
export NL_ETA_LEVELS=1.000,0.993,0.983,0.970,\
0.954,0.934,0.909,0.880,0.834,0.788,0.742,\
0.697,0.617,0.546,0.480,0.422,0.368,0.321,\
0.278,0.239,0.205,0.174,0.147,0.123,0.102,\
0.0832,0.0668,0.0526,0.0402,0.0295,0.0203,\
0.0124,0.00575,0.000
#
# PHYSICS NAMELIST:
export NL_MP_PHYSICS=2
export NL_RA_LW_PHYSICS=1
export NL_RA_SW_PHYSICS=2
export NL_RADT=40
export NL_SF_SFCLAY_PHYSICS=2
export NL_SF_SURFACE_PHYSICS=2
export NL_BL_PBL_PHYSICS=2
export NL_BLDT=0
export NL_CU_PHYSICS=3
export NL_CUDT=0
export NL_ISFFLX=1
export NL_IFSNOW=0
export NL_ICLOUD=0
export NL_SURFACE_INPUT_SOURCE=1
export NL_NUM_SOIL_LAYERS=4
export NL_SF_URBAN_PHYSICS=0
export NL_MAXIENS=1
export NL_MAXENS=3
export NL_MAXENS2=3
export NL_MAXENS3=16
export NL_ENSDIM=144
export NL_MP_ZERO_OUT=2
export NL_CU_RAD_FEEDBACK=.false.
export NL_CU_DIAG=1
export NL_PROGN=0
export NL_CUGD_AVEDX=1
#
# DYNAMICS NAMELIST:
export NL_USE_BASEPARAM_FR_NML=".true."
export NL_W_DAMPING=1
export NL_DIFF_OPT=1
export NL_KM_OPT=4
export NL_DIFF_6TH_OPT=0
export NL_DIFF_6TH_FACTOR=0.12
export NL_BASE_TEMP=290.
export NL_DAMP_OPT=3
export NL_ZDAMP=5000
export NL_DAMPCOEF=0.01
export NL_ISO_TEMP=0.
export NL_KHDIF=0
export NL_KVDIF=0
export NL_NON_HYDROSTATIC=.true.
export NL_TIME_STEP_SOUND=6
export NL_RK_ORD=3
export NL_MOIST_ADV_OPT=2
export NL_SCALAR_ADV_OPT=2
export NL_CHEM_ADV_OPT=2
export NL_TKE_ADV_OPT=2
#
# BDY_CONTROL NAMELIST:
export NL_SPEC_BDY_WIDTH=5
export NL_SPEC_ZONE=1
export NL_RELAX_ZONE=4
export NL_SPECIFIED=".true."
export NL_NESTED=".false."
export NL_REAL_DATA_INIT_TYPE=3
#
# QUILT NAMELIST:
export NL_NIO_TASKS_PER_GROUP=0
export NL_NIO_GROUPS=1
#
# FDDA NAMELIST:
export NL_GRID_FDDA=0
export NL_GFDDA_INNAME="'"wrffdda_d\<domain\>"'"
export NL_GFDDA_END_H=960
export NL_GFDDA_INTERVAL_M=360
export NL_FGDT=0
export NL_IF_NO_PBL_NUDGING_UV=0
export NL_IF_NO_PBL_NUDGING_T=0
export NL_IF_NO_PBL_NUDGING_q=0
export NL_IF_ZFAC_UV=0
export NL_K_ZFAC_UV=10
export NL_IF_ZFAC_T=0
export NL_K_ZFAC_T=10
export NL_IF_ZFAC_Q=0
export NL_K_ZFAC_Q=10
export NL_GUV=0.0012
export NL_GT=0.0012
export NL_GQ=0.0012
export NL_IF_RAMPING=1
export NL_DTRAMP_MIN=60.0
export NL_IO_FORM_GFDDA=2
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
#########################################################################
#
# RUN GEOGRID
#
#########################################################################
#
if [[ ${RUN_GEOGRID} = "true" ]]; then
   mkdir -p ${RUN_DIR}/geogrid
   cd ${RUN_DIR}/geogrid
#
   cp ${WPS_DIR}/geogrid.exe ./.
   ${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist.ksh
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
   export FCST_RANGE=${LBC_END}
   . ${SCRIPTS_DIR}/da_get_date_range.ksh
   export NL_START_DATE="'"${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:00:00"'"
   export NL_END_DATE="'"${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:00:00"'"
   ${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist.ksh
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
   ln -fs ${GEOGRID_DIR}/geo_em.d01.nc ./.
   ln -fs ../ungrib/FILE:* ./.
   ln -fs ${WPS_DIR}/metgrid/METGRID.TBL.${METGRID_TABLE_TYPE} METGRID.TBL
   ln -fs ${WPS_DIR}/metgrid.exe .
#
   export FCST_RANGE=${LBC_END}
   . ${SCRIPTS_DIR}/da_get_date_range.ksh
   export NL_START_DATE="'"${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:00:00"'"
   export NL_END_DATE="'"${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:00:00"'"
   ${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist.ksh
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
# Run real
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
#
# LINK IN THE METGRID FILES
   export L_DATE=${DATE}
   export L_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${LBC_END} 2>/dev/null)
   while [[ ${L_DATE} -le ${L_END_DATE} ]] ; do
      export L_YYYY=$(echo $L_DATE | cut -c1-4)
      export L_MM=$(echo $L_DATE | cut -c5-6)
      export L_DD=$(echo $L_DATE | cut -c7-8)
      export L_HH=$(echo $L_DATE | cut -c9-10)
      export L_FILE_DATE=${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00.nc
      ln -sf ${RUN_DIR}/${DATE}/metgrid/met_em.d01.${L_FILE_DATE} ./.
      export Q_DATE=${L_DATE}
      export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${Q_DATE} ${LBC_FREQ} 2>/dev/null) 
   done
#
# LOOP THROUGH BDY TENDENCY TIMES FOR PERTURB_BC
   export DATE_SAV=${DATE}
   export L_DATE=${DATE}
   export L_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${FCST_PERIOD} 2>/dev/null)
   while [[ ${L_DATE} -le ${L_END_DATE} ]] ; do      
#
# CREATE WRF NAMELIST
      export DATE=${L_DATE}
      . ${SCRIPTS_DIR}/da_get_date_range.ksh
      . ${HYBRID_SCRIPTS_DIR}/da_create_wrf_namelist.ksh
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
      mv wrfinput_d01 wrfinput_d01.$(${BUILD_DIR}/da_advance_time.exe ${DATE} 0 -W 2>/dev/null)
      mv wrfbdy_d01 wrfbdy_d01.$(${BUILD_DIR}/da_advance_time.exe ${DATE} 0 -W 2>/dev/null)
      export Q_DATE=${L_DATE}
      export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${Q_DATE} ${LBC_FREQ} 2>/dev/null) 
   done   
   export DATE=${DATE_SAV}
fi
#
#########################################################################
#
# RUN PERT WRFCHEM MET IC
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
   cp ${OB_PREPBUFR_DIR}/${DATE}/prepbufr.gdas.${DATE}.wo40.be ob.bufr
   export NL_OB_FORMAT=1
   export FCST_RANGE=${FCST_PERIOD}
   . ${SCRIPTS_DIR}/da_get_date_range.ksh
   export NL_START_DATE="'"${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:00:00"'"
   export NL_END_DATE="'"${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:00:00"'"
#
# LOOP THROUGH ALL BDY TENDENCY TIMES
   export L_DATE=${DATE}
   export L_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${FCST_PERIOD} 2>/dev/null)
   while [[ ${L_DATE} -le ${L_END_DATE} ]] ; do
#
# SET WRFDA PARAMETERS
      export ANALYSIS_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 0 -W 2>/dev/null)
      export NL_ANALYSIS_DATE="'"${ANALYSIS_DATE}"'"
      export NL_TIME_WINDOW_MIN="'"$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} -${ASIM_WINDOW} -W 2>/dev/null)"'"
      export NL_TIME_WINDOW_MAX="'"$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} +${ASIM_WINDOW} -W 2>/dev/null)"'"
      export DA_INPUT_FILE=../real/wrfinput_d01
      export DA_INPUT_FILE=${DA_INPUT_FILE}.${ANALYSIS_DATE}
      export NL_ANALYSIS_TYPE="'"RANDOMCV"'"
      export NL_PUT_RAND_SEED=true
      cp ${DA_INPUT_FILE} fg
      cp ${BE_DIR}/be.dat.cv3 be.dat
      cp ${WRFVAR_DIR}/run/LANDUSE.TBL ./.
      cp ${WRFVAR_DIR}/var/da/da_wrfvar.exe ./.
#
# LOOP THROUGH ALL MEMBERS IN THE ENSEMBLE
      let MEM=1
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
         export NL_SEED_ARRAY1=$(${BUILD_DIR}/da_advance_time.exe ${DATE} 0 -f hhddmmyycc)
         export NL_SEED_ARRAY2=`echo ${MEM} \* 100000 | bc -l `
         . ${HYBRID_SCRIPTS_DIR}/da_create_wrfda_namelist.ksh
#
# JOB SCRIPT 
         if [[ -f job.ksh ]]; then rm -rf job.ksh; fi
         rm -rf *.err
         rm -rf *.out
         touch job.ksh
         RANDOM=$$
         export JOBRND=wrfda_$RANDOM
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
mpirun.lsf ./da_wrfvar.exe  > index.html 2>&1 
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
         RC=$?
         if [[ $RC != 0 ]]; then
            echo wrfda-randomcv failed with error $RC
            exit $RC
         fi
         cp wrfvar_output wrfinput_d01.${ANALYSIS_DATE}.${CMEM}
         let MEM=${MEM}+1
      done
      export Q_DATE=${L_DATE}
      export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${Q_DATE} ${LBC_FREQ} 2>/dev/null) 
   done
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
         export DA_BDY_FILE=${DA_BDY_PATH}/wrfbdy_d01.${ANALYSIS_DATE}
         cp ${DA_BDY_FILE} wrfbdy_this
      else
         export DA_BDY_PATH=${RUN_DIR}/${DATE}/real
         export DA_BDY_FILE=${DA_BDY_PATH}/wrfbdy_d01.${ANALYSIS_DATE}
         cp ${DA_BDY_FILE} wrfbdy_this
      fi
      ln -fs ${TRUNK_DIR}/${DART_VER}/models/wrf/work/pert_wrf_bc ./.
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
         export DA_INPUT_PATH=${RUN_DIR}/${DATE}/wrfche_met_ic
         ln -fs ${DA_INPUT_PATH}/wrfinput_d01.${ANALYSIS_DATE}.${CMEM} wrfinput_this 
         ln -fs ${DA_INPUT_PATH}/wrfinput_d01.${NEXT_ANALYSIS_DATE}.${CMEM} wrfinput_next 
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
      mv wrfbdy_this wrfbdy_d01.${ANALYSIS_DATE}.${CMEM}
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
   export FILE=wrfinput_d01
   rm -rf ${FILE}
   ln -sf ${REAL_DIR}/${FILE}.${FILE_DATE} ${FILE}   
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
domains = 1,
/
EOF
#
# RUN exo_coldens
   ./exo_coldens < exo_coldens.inp
#
# TEST WHETHER OUTPUT EXISTS
   export FILE=exo_coldens_d01
   if [[ ! -e ${FILE} ]]; then
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
   export FILE=wrfinput_d01
   rm -rf ${FILE}
   ln -sf ${REAL_DIR}/${FILE}.${FILE_DATE} ${FILE}   
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
domains = 1,
/
EOF
#
# RUN wesely
   ./wesely < wesely.inp
#
# TEST WHETHER OUTPUT EXISTS
   export FILE=wrf_season_wes_usgs_d01.nc
   if [[ ! -e ${FILE} ]]; then
      echo WESELY FAILED
      exit
   else
      echo WESELY SUCCESS
   fi
fi
#
#########################################################################
#
# RUN WRFCHEM BIO
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
      export FILE=wrfinput_d01
      rm -rf ${FILE}
      ln -sf ${REAL_DIR}/${FILE}.${L_FILE_DATE} ${FILE}   
      export FILE=wrfbiochemi_d01
      if [[ ${L_DATE} -eq ${DATE} ]]; then
         rm -rf ${FILE}
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
domains = 1,
start_lai_mnth = ${MM},
end_lai_mnth = ${MM}
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
#BSUB -W 00:10        
#BSUB -q regular
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
# Run advance_model and wait until job completes
      bsub -K < job.ksh 
#
# TEST WHETHER OUTPUT EXISTS
      export FILE=wrfbiochemi_d01
      if [[ ! -e ${FILE} ]]; then
         echo WRFCHEM_BIO FAILED
         exit
      else
         echo WRFCHEM_BIO SUCCESS
         mv ${FILE} ${FILE}_${L_FILE_DATE}
      fi
      export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 6 2>/dev/null)
   done
fi
#
#########################################################################
#
# RUN WRFCHEM FIRE
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
   export FILE=wrfinput_d01
   rm -rf ${FILE}
   ln -sf ${REAL_DIR}/${FILE}.${FILE_DATE} ${FILE}   
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
domains = 1,
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
      export FILE=wrffirechemi_d01_${L_FILE_DATE}
      if [[ ! -e ${FILE} ]]; then
         echo WRFFIRE FAILED
         exit
      else
         echo WRFFIRE SUCCESS
      fi
      export P_DATE=${L_DATE}
      export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} 1 2>/dev/null)
   done
fi
#
#########################################################################
#
# RUN WRFCHEM CHEMI
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
      cp ${ACD_DIR}/AVE_TEST_DATA/chem_static_100km/${L_YYYY}${L_MM}${L_DD}/wrfchemi_d01_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00 ./.
      export P_DATE=${L_DATE}
      export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} 1 2>/dev/null)
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
   cp ${PERT_CHEM_INPUT_DIR}/runICBC_parent_rt.ksh ./.
   cp ${PERT_CHEM_INPUT_DIR}/runICBC_setN_rt.ksh ./.
   cp ${PERT_CHEM_INPUT_DIR}/random.py ./.
   cp ${PERT_CHEM_INPUT_DIR}/run_mozbc_rt.csh ./.
   cp ${PERT_CHEM_INPUT_DIR}/mozbc-dart/mozbc ./.
   cp ${PERT_CHEM_INPUT_DIR}/set0 ./.
   cp ${PERT_CHEM_INPUT_DIR}/set00 ./.
#
# SELECT MOZART DATA FILE
  if [[ ${MM} -eq 06 ]]; then export MOZART_DATA=h0001.nc; fi
  if [[ ${MM} -eq 07 ]]; then export MOZART_DATA=h0002.nc; fi
#
# CREATE INPUT FILES
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
   cp ${MET_EM_DIR}/${DATE}/met_em.d01.*:00:00.nc ./.
   ./random.py ${MOZ_SPREAD} ${NUM_MEMBERS} ${PERT_CHEM_INPUT_DIR} ${RUN_DIR}/${DATE}/wrfchem_chem_icbc
   ./runICBC_parent_rt.ksh
   ./runICBC_setN_rt.ksh
exit
#
# LOOP THROUGH ALL MEMBERS IN THE ENSEMBLE
   let MEM=${MEM_START}
   while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
      export CMEM=e${MEM}
      if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
      if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
#
# COMBINE WRFCHEM WITH WRF
      export WRFINPEN=wrfinput_d01_${YYYY}-${MM}-${DD}_${HH}:00:00.${CMEM}
      export WRFBDYEN=wrfbdy_d01_${YYYY}-${MM}-${DD}_${HH}:00:00.${CMEM}
      if [[ -e wrfinput_d01 ]]; then rm -rf wrfinput_d01; fi
      if [[ -e ${WRFINPEN} ]]; then rm -rf ${WRFINPEN}; fi
      if [[ -e wrfbdy_d01 ]]; then rm -rf wrfbdy_d01; fi
      if [[ -e ${WRFBDYEN} ]]; then rm -rf ${WRFBDYEN}; fi
      cp ${WRFCHEM_MET_IC_DIR}/${DATE}/${WRFINPEN} wrfinput.met
      cp ${WRFCHEM_MET_BC_DIR}/${DATE}/${WRFBDYEN} wrfbdy.met
      ncks -x -v MU,PSFC,Q2,T2,TH2,U10,V10,P,QVAPOR,T,U,V ${WRFINPEN} wrfinput_d01
      ncks -v MU,PSFC,Q2,T2,TH2,U10,V10,P,QVAPOR,T,U,V wrfinput.met ${WRFINPEN}
      ncks -A wrfinput_d01 ${WRFINPEN}
#
      ncks -x -v MU_BTXE,MU_BTXS,MU_BTYE,MU_BTYS,MU_BXE,MU_BXS,MU_BYE,MU_BYS,QVAPOR_BTXE,QVAPOR_BTXS,QVAPOR_BTYE,QVAPOR_BTYS,QVAPOR_BXE,QVAPOR_BXS,QVAPOR_BYE,QVAPOR_BYS,T_BTXE,T_BTXS,T_BTYE,T_BTYS,T_BXE,T_BXS,T_BYE,T_BYS,U_BTXE,U_BTXS,U_BTYE,U_BTYS,U_BXE,U_BXS,U_BYE,U_BYS,V_BTXE,V_BTXS,V_BTYE,V_BTYS,V_BXE,V_BXS,V_BYE,V_BYS,PH_BTXE,PH_BTXS,PH_BTYE,PH_BTYS,PH_BXE,PH_BXS,PH_BYE,PH_BYS ${WRFBDYEN} wrfbdy_d01
      ncks -v MU_BTXE,MU_BTXS,MU_BTYE,MU_BTYS,MU_BXE,MU_BXS,MU_BYE,MU_BYS,QVAPOR_BTXE,QVAPOR_BTXS,QVAPOR_BTYE,QVAPOR_BTYS,QVAPOR_BXE,QVAPOR_BXS,QVAPOR_BYE,QVAPOR_BYS,T_BTXE,T_BTXS,T_BTYE,T_BTYS,T_BXE,T_BXS,T_BYE,T_BYS,U_BTXE,U_BTXS,U_BTYE,U_BTYS,U_BXE,U_BXS,U_BYE,U_BYS,V_BTXE,V_BTXS,V_BTYE,V_BTYS,V_BXE,V_BXS,V_BYE,V_BYS,PH_BTXE,PH_BTXS,PH_BTYE,PH_BTYS,PH_BXE,PH_BXS,PH_BYE,PH_BYS wrfbdy.met ${WRFBDYEN}
      ncks -A wrfbdy_d01 ${WRFBDYEN}
      rm -rf wrfinput.met
      rm -rf wrfbdy.met
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
   if [[ -e perturb_chem_emiss.exe ]]; then rm -rf perturb_chem_emiss.exe; fi
   cp ${PERT_CHEM_EMISS_DIR}/perturb_chem_emiss.exe ./.
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
# GET EMISSON FILES FOR THIS MEMBER
         export WRFCHEMI=wrfchemi_d01_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
         export WRFFIRECHEMI=wrffirechemi_d01_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
         export WRFBIOCHEMI=wrfbiochemi_d01_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
#
         if [[ ${NL_PERT_CHEM} == true ]]; then
            cp ${WRFCHEM_CHEMI_DIR}/${WRFCHEMI} ./.
            cp ${WRFCHEMI} ${WRFCHEMI}.${CMEM}
         fi
         if [[ ${NL_PERT_FIRE} == true ]]; then
            cp ${WRFCHEM_FIRE_DIR}/${WRFFIRECHEMI} ./.
            cp ${WRFFIRECHEMI} ${WRFFIRECHEMI}.${CMEM}
         fi
         if [[ ${NL_PERT_BIO} == true ]]; then
            cp ${WRFCHEM_BIO_DIR}/${WRFBIOCHEMI} ./.
            cp ${WRFBIOCHEMI} ${WRFBIOCHEMI}.${CMEM}
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
         ./perturb_chem_emiss.exe
#
# GO TO NEXT MEMBER
         let MEM=${MEM}+1
      done
# ADVANCE TIME
      export P_DATE=${L_DATE}      
      export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} 1 2>/dev/null)
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
   export L_DATE=${D_YYYY}${D_MM}${D_DD}06
   export E_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} +24 2>/dev/null)
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
   cp ${WRFCHEM_MET_IC_DIR}/wrfinput_d01.${FILE_DATE}.e001 wrfinput_d01
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
# Run wrf_to_dart
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
