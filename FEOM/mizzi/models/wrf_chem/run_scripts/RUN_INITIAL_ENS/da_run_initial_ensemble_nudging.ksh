#!/bin/ksh -aeux
#########################################################################
# Script: APM_run_wps.ksh
#
# Purpose: Script to run WPS, REAL, WRFDA-RANDOMCV, and PERTURB_WRF_BC
#          for creating WRF ensemble input and bdy files.
#
#########################################################################
#
# SELECT RUN OPTIONS
export RUN_GEOGRID=false
export RUN_UNGRIB=false
export RUN_METGRID=false
export RUN_REAL=true
export RUN_RANDOMCV=false
export RUN_PERT_WRF_BC=false
#
# CODE VERSIONS
export WPS_VER=WPSv3.4_dmpar
export WRFDA_VER=WRFDAv3.4_dmpar
export WRFDA_TOOLS_VER=WRFDA_TOOLSv3.4
export WRF_VER=WRFv3.4_dmpar
export DART_VER=DART_DEVEL
#
# EXPERIMENT DETAILS:
export INITIAL_DATE=2012063000
export FINAL_DATE=2012070506
export NUM_MEMBERS=2
export NUM_PROCS=8
export MEM_START=1
export REGION=colo_fires_12km_nudging
# NAM forecast data
export FG_TYPE=NAM
export GRIB_PART1=nam_218_
export GRIB_PART2=.grb
# GFS forecast data
#export FG_TYPE=GFS
#export GRIB_PART1=gfs_4_
#export GRIB_PART2=.grb2
#
export USE_PREPBUFR=true
export HOR_SCALE=1500
#
# TIME DATA:
export DATE=${INITIAL_DATE}
export START_DATE=${INITIAL_DATE}
export END_DATE=${FINAL_DATE}
export YYYY=$(echo $DATE | cut -c1-4)
export MM=$(echo $DATE | cut -c5-6)
export DD=$(echo $DATE | cut -c7-8)
export HH=$(echo $DATE | cut -c9-10)
#
# CYCLING AND BOUNDARY CONDITION TIME DATA
export CYCLE_PERIOD=6
export FCST_RANGE=6
export ASIM_WINDOW=3
export LBC_FREQ=3
(( INTERVAL_SECONDS=${LBC_FREQ}*60*60 ))
export LBC_START=0
#
# DIRECTORIES:
export RUN_DIR=/Volumes/disk2/${REGION}
export TRUNK_DIR=/Volumes/data1/TRUNK
export HYBRID_TRUNK_DIR=${TRUNK_DIR}/HYBRID_TRUNK
export HYBRID_TRUNK_DATA_DIR=/Volumes/disk4/HYBRID_TRUNK_DATA
export WPS_GEOG_DIR=/Volumes/disk4/WPS_GEOG
export DAT_DIR=/Volumes/disk4/COLO_FIRE_DATA
#
export HYBRID_SCRIPTS_DIR=${TRUNK_DIR}/HYBRID_TRUNK/hybrid_scripts
export GRIB_DIR=${DAT_DIR}/nam_forecasts
if [[ ${FG_TYPE} == GFS ]]; then export GRIB_DIR=${DAT_DIR}/gfs_forecasts; fi
export SCRIPTS_DIR=${TRUNK_DIR}/${WRFDA_TOOLS_VER}/scripts
export WPS_DIR=${TRUNK_DIR}/${WPS_VER}
export WRFVAR_DIR=${TRUNK_DIR}/${WRFDA_VER}
export BUILD_DIR=${WRFVAR_DIR}/var/da
export WRF_DIR=${TRUNK_DIR}/${WRF_VER}
export OBSPROC_DIR=${WRFVAR_DIR}/var/obsproc
export VTABLE_DIR=${WPS_DIR}/ungrib/Variable_Tables
export OB_PREPBUFR_DIR=${DAT_DIR}/obs
export OB_ASCII_DIR=${DAT_DIR}/obs
export BE_DIR=${WRFVAR_DIR}/var/run
export FIXED_DATA_DIR=${HYBRID_TRUNK_DATA_DIR}/hybrid_data/fixed_data
export DART_INPUT_DIR=${HYBRID_TRUNK_DATA_DIR}/hybrid_data/dart_input
#
# SHARE NAMELIST PARAMETERS:
export NL_WRF_CORE="'"ARW"'"
export NL_MAX_DOM=1
export NL_IO_FORM_GEOGRID=2
export NL_OPT_OUTPUT_FROM_GEOGRID_PATH="'"${RUN_DIR}/rc"'"
export NL_DEBUG_LEVEL=0
export VTABLE_TYPE=NAM
if [[ ${FG_TYPE} == GFS ]]; then export VTABLE_TYPE=GFS; fi
export METGRID_TABLE_TYPE=ARW
#
# GEOGRID NAMELIST PARAMETERS:
export NL_PARENT_ID=1
export NL_PARENT_GRID_RATIO=1
export NL_I_PARENT_START=1
export NL_J_PARENT_START=1
export NL_S_WE=1
export NL_E_WE=249
export NL_S_SN=1
export NL_E_SN=195
export NL_GEOG_DATA_RES="'"30s"'"
export NL_DX=12000
export NL_DY=12000
export NL_MAP_PROJ="'"lambert"'"
export NL_REF_LAT=39.0
export NL_REF_LON=-111.0
export NL_TRUELAT1=30.0
export NL_TRUELAT2=60.0
export NL_STAND_LON=-111.0
export NL_GEOG_DATA_PATH="'"${WPS_GEOG_DIR}"'"
export NL_OPT_GEOGRID_TBL_PATH="'"${WPS_DIR}/geogrid"'"
#
# UNGRID NAMELIST PARAMETERS
export NL_OUT_FORMAT="'"WPS"'"
#
# METGRID NAMELIST PARAMETERS
export NL_IO_FORM_METGRID=2
export NL_OPT_IGNORE_DOM_CENTER=.false.
export NL_OPT_OUTPUT_FROM_METGRID_PATH="'"${RUN_DIR}/rc/2007081500"'"
#
# WRF NAMELIST PARAMETERS
#
# TIME CONTROL
export NL_RUN_DAYS=0
export NL_RUN_HOURS=${FCST_RANGE}
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
# DOMAINS
export NL_TIME_STEP=20
export NL_TIME_STEP_FRACT_NUM=0
export NL_TIME_STEP_FRACT_DEN=1
export NL_MAX_DOM=1
export NL_S_WE=1
export NL_E_WE=249
export NL_S_SN=1
export NL_E_SN=195
#
# This is number of vertical levels in the wrfinput data
export NL_E_VERT=28
export NL_P_TOP_REQUESTED=5000
export NL_INTERP_TYPE=1
export NL_T_EXTRAP_TYPE=1
#
# This in number of vertical levels in the initial data
export NL_NUM_METGRID_LEVELS=40
if [[ ${FG_TYPE} == GFS ]]; then export NL_NUM_METGRID_LEVELS=27; fi
export NL_NUM_METGRID_SOIL_LEVELS=4
export NL_DX=12000.0
export NL_DY=12000.0
export NL_GRID_ID=1
export NL_PARENT_ID=0
export NL_I_PARENT_START=0
export NL_J_PARENT_START=0
export NL_PARENT_GRID_RATIO=1
export NL_PARENT_TIME_STEP_RATIO=1
export NL_FEEDBACK=0
export NL_SMOOTH_OPTION=1
export NL_ETA_LEVELS=1.000,0.990,0.978,0.964,0.946,\
0.922,0.894,0.860,0.817,0.766,\
0.707,0.644,0.576,0.507,0.444,\
0.380,0.324,0.273,0.228,0.188,\
0.152,0.121,0.093,0.069,0.048,\
0.029,0.014,0.000
#
# PHYSICS
export NL_MP_PHYSICS=8
export NL_RA_LW_PHYSICS=1
export NL_RA_SW_PHYSICS=1
export NL_RADT=12
export NL_SF_SFCLAY_PHYSICS=1
export NL_SF_SURFACE_PHYSICS=2
export NL_BL_PBL_PHYSICS=1
export NL_BLDT=0
export NL_CU_PHYSICS=5
export NL_CUDT=0
export NL_ISFFLX=1
export NL_IFSNOW=1
export NL_ICLOUD=1
export NL_SURFACE_INPUT_SOURCE=1
export NL_NUM_SOIL_LAYERS=4
export NL_SF_URBAN_PHYSICS=1
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
# DYNAMICS
export NL_USE_BASEPARAM_FR_NML=".true."
export NL_W_DAMPING=0
export NL_DIFF_OPT=1
export NL_KM_OPT=4
export NL_DIFF_6TH_OPT=0
export NL_DIFF_6TH_FACTOR=0.12
export NL_BASE_TEMP=290.
export NL_DAMP_OPT=3
export NL_ZDAMP=5000
export NL_DAMPCOEF=0.2
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
# BDY_CONTROL
export NL_SPEC_BDY_WIDTH=5
export NL_SPEC_ZONE=1
export NL_RELAX_ZONE=4
export NL_SPECIFIED=".true."
export NL_NESTED=".false."
export NL_REAL_DATA_INIT_TYPE=3
#
# NAMELIST_QUILT
export NL_NIO_TASKS_PER_GROUP=0
export NL_NIO_GROUPS=1
#
# FDDA
export NL_GRID_FDDA=2
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
#
# WRFVAR1
export NL_PRINT_DETAIL_GRAD=false
export NL_VAR4D=false
export NL_MULTI_INC=0
#
# WRFVAR3
export NL_OB_FORMAT='set below'
export NL_NUM_FGAT_TIME=1
#
# WRFVAR4
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
# WRFVAR5
export NL_CHECK_MAX_IV=true
export NL_PUT_RAND_SEED='set below'
#
# WRFVAR6
export NL_NTMAX=100
#
# WRFVAR7
export NL_VAR_SCALING4=1.0
export NL_JE_FACTOR=1.0
export NL_CV_OPTIONS=3
#
# WRFVAR11
export NL_CV_OPTIONS_HUM=1
export NL_CHECK_RH=2
export NL_SEED_ARRAY1='set below'
export NL_SEED_ARRAY2='set=below'
export NL_CALCULATE_CG_COST_FN=true
export NL_LAT_STATS_OPTION=false
#
# WRFVAR15
export NL_NUM_PSEUDO=0
export NL_PSEUDO_X=0
export NL_PSEUDO_Y=0
export NL_PSEUDO_Z=0
export NL_PSEUDO_ERR=0.0
export NL_PSEUDO_VAL=0.0
#
# WRFVAR16
export NL_ALPHACV_METHOD=2
export NL_ENSDIM_ALPHA=0
export NL_ALPHA_CORR_TYPE=3
export NL_ALPHA_CORR_SCALE=${HOR_SCALE}
export NL_ALPHA_STD_DEV=1.0
export NL_ALPHA_VERTLOC=false
export NL_ALPHA_TRUNCATION=1
#
# WRFVAR17
export NL_ANALYSIS_TYPE='set below'
#
# WRFVAR18
export NL_ANALYSIS_DATE='set below'
#
# WRFVAR19
export NL_PSEUDO_VAR="'"t"'"
#
# WRFVAR21
export NL_TIME_WINDOW_MIN='set below'
#
# WRFVAR22
export NL_TIME_WINDOW_MAX='set below'
#
# WRFVAR23
export NL_JCDFI_USE=false
export NL_JCDFI_IO=false
#
#
# SET TIME DATA
export DATE=${INITIAL_DATE}
. ${SCRIPTS_DIR}/da_get_date_range.ksh
#
# RUN GEOGRID
if [[ ${RUN_GEOGRID} = "true" ]]; then
   mkdir -p ${RUN_DIR}/${DATE}/wps/geogrid
   cd ${RUN_DIR}/${DATE}/wps/geogrid
   export GEOGRID_PATH=${RUN_DIR}/${DATE}/wps/geogrid
#
# CREATE WPS NAMELIST
   ${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist.ksh
#
# RUN SCRIPT
   ln -sf ${WPS_DIR}/geogrid.exe ./.
   mpirun -np ${NUM_PROCS} ./geogrid.exe
   RC=$?
   if [[ $RC != 0 ]]; then
      echo geogrid failed with error $RC
      exit $RC
   fi
fi
#
# LOOP THROUGH DATES TO CREATE MET FILES
while [[ ${DATE} -le ${FINAL_DATE} ]] ; do 
export YYYY=$(echo $DATE | cut -c1-4)
export MM=$(echo $DATE | cut -c5-6)
export DD=$(echo $DATE | cut -c7-8)
export HH=$(echo $DATE | cut -c9-10)
#
# CREATE UNGRIB WORKING DIRECTORY
   if [[ ${RUN_UNGRIB} = "true" ]]; then 
      mkdir -p ${RUN_DIR}/${DATE}/wps/ungrib
      cd ${RUN_DIR}/${DATE}/wps/ungrib
#
# LINK VTABLE
      ln -fs ${VTABLE_DIR}/Vtable.${VTABLE_TYPE} Vtable
#
# CREATE WPS NAMELIST
      . ${SCRIPTS_DIR}/da_get_date_range.ksh
      export NL_START_DATE="'"${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:00:00"'"
      export NL_END_DATE="'"${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:00:00"'"
      ${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist.ksh
#
# CREATE UNGRIB FILE LIST
      (( LBC_END=${FCST_RANGE} ))
      export START_HOUR=${NL_START_HOUR}00
      integer ICNT=0
      FILES=''
      if [[ -e ${GRIB_DIR}/${DATE} ]]; then
         (( LBC_ITR=${LBC_START} ))
         while [[ ${LBC_ITR} -le ${LBC_END} ]]; do
               if [[ ${LBC_ITR} -lt 1000 ]]; then export CMEM=${LBC_ITR}; fi
               if [[ ${LBC_ITR} -lt 100  ]]; then export CMEM=0${LBC_ITR}; fi
               if [[ ${LBC_ITR} -lt 10   ]]; then export CMEM=00${LBC_ITR}; fi
               if [[ ${LBC_ITR} -eq 0    ]]; then export CMEM=000; fi
               export FILE=${GRIB_DIR}/${DATE}/${GRIB_PART1}${YYYY}${MM}${DD}_${START_HOUR}_${CMEM}${GRIB_PART2}
               FILES="${FILES} ${FILE}"
            ICNT=ICNT+1
            (( LBC_ITR=${LBC_ITR}+${LBC_FREQ} ))
         done
      fi
#
# LINK GRIB FILES
      ${WPS_DIR}/link_grib.csh $FILES
#
# RUN UNGRIB
#      ln -fs ${WPS_DIR}/ungrib.exe ./.
      cp ${WPS_DIR}/ungrib.exe ./.
      pwd
      ./ungrib.exe
      RC=$?
      if [[ $RC != 0 ]]; then
         echo ungrib failed with error $RC
         exit $RC
      fi
   fi
#
# RUN METGRID
   if [[ ${RUN_METGRID} = "true" ]]; then 
      mkdir -p ${RUN_DIR}/${DATE}/wps/metgrid
      cd ${RUN_DIR}/${DATE}/wps/metgrid
#      ln -fs ${FIXED_DATA_DIR}/geo_em.d01.nc ./.
      ln -fs ${GEOGRID_PATH}/geo_em.d01.nc ./.
      ln -fs ../ungrib/namelist.wps ./.
      ln -fs ../ungrib/FILE:* ./.
      ln -fs $WPS_DIR/metgrid/METGRID.TBL.${METGRID_TABLE_TYPE} METGRID.TBL
      ln -fs $WPS_DIR/metgrid.exe .
      mpirun -np ${NUM_PROCS} ./metgrid.exe
      RC=$?
      if [[ $RC != 0 ]]; then
         echo metgrid failed with error $RC
         exit $RC
      fi
   fi
#
# RUN REAL
   if [[ ${RUN_REAL} = "true" ]]; then 
#
# LOOP THROUGH ALL BDY TENDENCY TIMES FOR THIS FORECAST START DATE. 
      export FCST_RANGE_SAVE=${FCST_RANGE}
      export DATE_SAVE=${DATE}
      export L_DATE=${DATE}      
      export L_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} ${FCST_RANGE} 2>/dev/null)
      while [[ ${L_DATE} -le ${L_END_DATE} ]] ; do 
         mkdir -p ${RUN_DIR}/${DATE_SAVE}/real/${L_DATE}
         cd ${RUN_DIR}/${DATE_SAVE}/real/${L_DATE}
#
# CREATE WRF NAMELIST
         if [[ ${L_DATE} = ${DATE_SAVE} ]]; then
            export FCST_RANGE=${FCST_RANGE_SAVE}
         else
#            export FCST_RANGE=${LBC_FREQ}
            export FCST_RANGE=0
         fi        
         export DATE=${L_DATE}
         . ${SCRIPTS_DIR}/da_get_date_range.ksh
         . ${HYBRID_SCRIPTS_DIR}/da_create_wrf_namelist.ksh
#
# APM Begin special loop for PRE-EXISTING MET_EM DATA (NOTE: the 24 hr is ad hoc and depends on the forecast period)
         export LLL_DATE=${L_DATE}
         export LLL_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${LLL_DATE} 0 2>/dev/null)
         while [[ ${LLL_DATE} -le ${LLL_END_DATE} ]] ; do
            export YY=$(echo $LLL_DATE | cut -c1-4)
            export MM=$(echo $LLL_DATE | cut -c5-6)
            export DD=$(echo $LLL_DATE | cut -c7-8)
            export HH=$(echo $LLL_DATE | cut -c9-10)
            export LLL_FILE_DATE=${YY}-${MM}-${DD}_${HH}:00:00.nc
#            ln -sf ${DAT_DIR}/met_em_files/${LLL_DATE}/met_em.d01.${LLL_FILE_DATE} ./.
            ln -sf ${DAT_DIR}/met_em_files/${LLL_DATE}/met_em.d01.* ./.
            export LLQ_DATE=${LLL_DATE}
            export LLL_DATE=$(${BUILD_DIR}/da_advance_time.exe ${LLQ_DATE} ${LBC_FREQ} 2>/dev/null) 
         done
# APM End special loop for PRE-EXISTING MET_EM DATA
#
         ln -fs ${WRF_DIR}/main/real.exe ./.
         mpirun -np ${NUM_PROCS} ./real.exe
         RC=$?
         if [[ $RC != 0 ]]; then
            echo real failed with error $RC
            exit $RC
         fi
         if [[ ${L_DATE} = ${DATE_SAVE} ]]; then
            mv wrfinput_d01 wrfinput_d01.$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 0 -W 2>/dev/null)
            mv wrfbdy_d01 wrfbdy_d01.$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 0 -W 2>/dev/null)
            mv wrffdda_d01 wrffdda_d01.$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 0 -W 2>/dev/null)
         else
            mv wrfinput_d01 wrfinput_d01.$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 0 -W 2>/dev/null)
         fi        
         export NEXT_L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${LBC_FREQ} 2>/dev/null)   
         export L_DATE=${NEXT_L_DATE}
      done
      export DATE=${DATE_SAVE}
      export FCST_RANGE=${FCST_RANGE_SAVE}
   fi
#
# RUN RANDOMCV
   if [[ ${RUN_RANDOMCV} = "true" ]]; then 
#
# LOOP THROUGH ALL MEMBERS IN THE ENSEMBLE
      let MEM=${MEM_START}
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
#
# LOOP THROUGH ALL BDY TENDENCY TIMES FOR THIS MEMBER. 
         export L_DATE=${DATE}
         export L_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} ${FCST_RANGE} 2>/dev/null)
         while [[ ${L_DATE} -le ${L_END_DATE} ]]; do
            if [[ ! -d ${RUN_DIR}/${DATE}/wrfda/${L_DATE} ]]; then
               mkdir -p ${RUN_DIR}/${DATE}/wrfda/${L_DATE}
               cd ${RUN_DIR}/${DATE}/wrfda/${L_DATE}
            else
               cd ${RUN_DIR}/${DATE}/wrfda/${L_DATE}
            fi
#
# SET OBS PATH (DO I NEED OBS FOR RANDOMCV?  THESE OBS MAY BE AT THE WRONG TIME)
            if [[ ${USE_PREPBUFR} = "true" ]]; then
               cp ${OB_PREPBUFR_DIR}/${DATE}/ob.bufr_littleendian ob.bufr
               export NL_OB_FORMAT=1
            else
               cp ${OB_ASCII_DIR}/${DATE}/ob.ascii ob.ascii
               export NL_OB_FORMAT=2
            fi
#
# SET WRFDA PARAMETERS
            export ANALYSIS_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 0 -W 2>/dev/null)
            export NL_ANALYSIS_DATE="'"${ANALYSIS_DATE}"'"
            export NL_TIME_WINDOW_MIN="'"$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} -${ASIM_WINDOW} -W 2>/dev/null)"'"
            export NL_TIME_WINDOW_MAX="'"$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} +${ASIM_WINDOW} -W 2>/dev/null)"'"
            export NL_SEED_ARRAY1=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 0 -f hhddmmyycc)
            export NL_SEED_ARRAY2=`echo ${MEM} \* 100000 | bc -l `
            export DA_INPUT_FILE=../../real/${L_DATE}/wrfinput_d01
            export DA_INPUT_FILE=${DA_INPUT_FILE}.${ANALYSIS_DATE}
            export NL_ANALYSIS_TYPE="'"RANDOMCV"'"
            export NL_PUT_RAND_SEED=true
            cp ${DA_INPUT_FILE} fg
            ln -sf ${BE_DIR}/be.dat.cv3 be.dat
            ln -sf ${WRFVAR_DIR}/run/LANDUSE.TBL ./.
            ln -sf ${WRFVAR_DIR}/var/da/da_wrfvar.exe ./.
            . ${HYBRID_SCRIPTS_DIR}/da_create_wrfda_namelist.ksh
            mpirun -np ${NUM_PROCS} ./da_wrfvar.exe
            RC=$?
            if [[ $RC != 0 ]]; then
               echo wrfda-randomcv failed with error $RC
               exit $RC
            fi
            cp wrfvar_output wrfinput_d01.${ANALYSIS_DATE}.${CMEM}
            export NEXT_L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${LBC_FREQ} 2>/dev/null)
            export L_DATE=${NEXT_L_DATE}
         done
         let MEM=${MEM}+1
      done
   fi
#
# RUN PERTURB_WRF_BC
   if [[ ${RUN_PERT_WRF_BC} = "true" ]]; then 
         export APM_DUMMY=1
         cd ${RUN_DIR}/${DATE}/wrfda/${DATE}
#
# LOOP THROUGH ALL MEMBERS IN THE ENSEMBLE
      let MEM=${MEM_START}
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
#
# LOOP THROUGH ALL BDY TENDENCY TIMES FOR THIS MEMBER.
         export L_DATE=${DATE}
         export L_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} ${FCST_RANGE} 2>/dev/null)
         if [[ -f wrfbdy_this ]]; then
            rm wrfbdy_this
         fi
         while [[ ${L_DATE} -lt ${L_END_DATE} ]]; do
#
# SET OBS PATH (DO I NEED OBS FOR RANDOMCV?  THESE OBS MAY BE AT THE WRONG TIME)
            if [[ ${USE_PREPBUFR} = "true" ]]; then
               cp ${OB_PREPBUFR_DIR}/${DATE}/ob.bufr_littleendian ob.bufr
               export NL_OB_FORMAT=1
            else
               cp ${OB_ASCII_DIR}/${DATE}/ob.ascii ob.ascii
               export NL_OB_FORMAT=2
            fi
            export ANALYSIS_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 0 -W 2>/dev/null)
            export NL_ANALYSIS_DATE="'"${ANALYSIS_DATE}"'"
            export NL_TIME_WINDOW_MIN="'"$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} -${ASIM_WINDOW} -W 2>/dev/null)"'"
            export NL_TIME_WINDOW_MAX="'"$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} +${ASIM_WINDOW} -W 2>/dev/null)"'"
            export NL_SEED_ARRAY1=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 0 -f hhddmmyycc)
            export NL_SEED_ARRAY2=`echo ${MEM} \* 100000 | bc -l `
            export DA_INPUT_FILE=../../real/${L_DATE}/wrfinput_d01
            export DA_INPUT_FILE=${DA_INPUT_FILE}.${ANALYSIS_DATE}
            export NL_ANALYSIS_TYPE="'"RANDOMCV"'"
            export NL_PUT_RAND_SEED=true
            export NEXT_L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${LBC_FREQ} 2>/dev/null)
            export NEXT_ANALYSIS_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${LBC_FREQ} -W 2>/dev/null)
            export DA_INPUT_PATH=${RUN_DIR}/${DATE}/wrfda/${L_DATE}
            export DA_NEXT_INPUT_PATH=${RUN_DIR}/${DATE}/wrfda/${NEXT_L_DATE}
            if [[ ! -f wrfbdy_this ]]; then
               export DA_BDY_PATH=${RUN_DIR}/${DATE}/real/${DATE}
               export DA_BDY_FILE=${DA_BDY_PATH}/wrfbdy_d01.${ANALYSIS_DATE}
               cp ${DA_BDY_FILE} wrfbdy_this
            fi
            ln -fs ${DA_INPUT_PATH}/wrfinput_d01.${ANALYSIS_DATE}.${CMEM} wrfinput_this 
            ln -fs ${DA_NEXT_INPUT_PATH}/wrfinput_d01.${NEXT_ANALYSIS_DATE}.${CMEM} wrfinput_next 
            ln -fs ${DART_INPUT_DIR}/input.nml ./.
            ln -fs ${TRUNK_DIR}/${DART_VER}/models/wrf/work/pert_wrf_bc ./.
            . ${HYBRID_SCRIPTS_DIR}/da_create_wrfda_namelist.ksh
            ./pert_wrf_bc > pert_wrf_bc.out.${CMEM} 2>&1
            export L_DATE=${NEXT_L_DATE}
         done
         export ANALYSIS_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} 0 -W 2>/dev/null)
         mv wrfbdy_this wrfbdy_d01.${ANALYSIS_DATE}.${CMEM}
         let MEM=${MEM}+1
      done
# 
# COPY ENSEMBLE DATA TO WPB_RC
      export WPB_RC_DIR=${RUN_DIR}/wpb_rc/${DATE}
      export RC_DIR=${RUN_DIR}/rc/${DATE}
      mkdir -p ${WPB_RC_DIR}
      mkdir -p ${RC_DIR}
      export ANALYSIS_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} 0 -W 2>/dev/null)
      cp ${RUN_DIR}/${DATE}/real/${DATE}/wrfinput_d01.${ANALYSIS_DATE} ${RC_DIR}/wrfinput_d01_${ANALYSIS_DATE} 
      cp ${RUN_DIR}/${DATE}/real/${DATE}/wrfbdy_d01.${ANALYSIS_DATE} ${RC_DIR}/wrfbdy_d01_${ANALYSIS_DATE}
      let MEM=${MEM_START}
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
         cp wrfinput_d01.${ANALYSIS_DATE}.${CMEM} ${WPB_RC_DIR}/wrfinput_d01_${ANALYSIS_DATE}.${CMEM} 
         cp wrfbdy_d01.${ANALYSIS_DATE}.${CMEM} ${WPB_RC_DIR}/wrfbdy_d01_${ANALYSIS_DATE}.${CMEM} 
         let MEM=${MEM}+1
      done
   fi
   export NEXT_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} ${CYCLE_PERIOD} 2>/dev/null)
   export DATE=${NEXT_DATE}
done

