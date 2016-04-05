#!/bin/ksh -x
#########################################################################
# Script: APM_run_wps.ksh
#
# Purpose: Script to run WPS, REAL, WRFDA-RANDOMCV, and PERTURB_WRF_BC
#          for creating WRF ensemble input and bdy files.
#
#########################################################################
#
# SELECT RUN OPTIONS
export RUN_GEOGRID=true
export RUN_UNGRIB=true
export RUN_METGRID=true
export RUN_REAL=true
export RUN_RANDOMCV=true
export RUN_PERT_WRF_BC=true
#
# VERSIONS
export WRF_VER=WRFv3.3_dmpar
export WRFDA_VER=WRFDAv3.3_dmpar
export WRFDA_TOOLS_VER=WRFDA_TOOLSv3.3
export WPS_VER=WPSv3.3 
export DART_VER=DART
export MACHINE=MAC
#
# EXPERIMENT DETAILS:
export NUM_MEMBERS=80
export NUM_NODES=1
export NUM_PROCS=8
export MEM_START=1
export REGION=conus_200km
export EXPT=PAULA_INITIAL_DATA
export FG_TYPE=gfs
export GRIB_PART1=gfs.t
export GRIB_PART2=z.pgrb2f
export USE_PREPBUFR=false
export HOR_SCALE=1500
#
# TIME DATA:
export INITIAL_DATE=2010100100
export FINAL_DATE=2010100100
export DATE=${INITIAL_DATE}
export YYYY=$(echo $DATE | cut -c1-4)
export MM=$(echo $DATE | cut -c5-6)
export DD=$(echo $DATE | cut -c7-8)
export HH=$(echo $DATE | cut -c9-10)
#
# CYCLING AND BOUNDARY CONDITION TIME DATA
export CYCLE_PERIOD=6
export FCST_RANGE=120
export ASIM_WINDOW=3
export LBC_FREQ=3
(( INTERVAL_SECONDS=${LBC_FREQ}*60*60 ))
export LBC_START=0
#
# DIRECTORIES:
export RUN_DIR=/Volumes/cassia3/mizzi/${REGION}/${EXPT}
export TRUNK_DIR=/data2/da/mizzi/code/TRUNK
export WPS_GEOG_DIR=/gpfs01/home/fgao/DATA/WPS_GEOG
export GFS_DAT_DIR=/Volumes/cassia3/mizzi/DATA/PAULA_DATA
export PAULA_DATA_DIR=/Volumes/cassia3/mizzi/DATA/PAULA_DATA

export HYBRID_SCRIPTS_DIR=${TRUNK_DIR}/HYBRID_TRUNK/hybrid_scripts
export DAT_DIR=${TRUNK_DIR}/HYBRID_TRUNK/hybrid_data
export GRIB_DIR=${GFS_DAT_DIR}/gfs
export SCRIPTS_DIR=${TRUNK_DIR}/${WRFDA_TOOLS_VER}/scripts
export WPS_DIR=${TRUNK_DIR}/${WPS_VER}
export WRFVAR_DIR=${TRUNK_DIR}/${WRFDA_VER}
export WRF_DIR=${TRUNK_DIR}/${WRF_VER}
export DART_DIR=${TRUNK_DIR}/${DART_VER}
export OBSPROC_DIR=${WRFVAR_DIR}/var/obsproc
export VTABLE_DIR=${WPS_DIR}/ungrib/Variable_Tables
export OB_PREPBUFR_DIR=${PAULA_DATA_DIR}/ob
export OB_ASCII_DIR=${PAULA_DATA_DIR}/ob
export BE_DIR=${PAULA_DATA_DIR}/be
export FIXED_DATA_DIR=${DAT_DIR}/fixed_data
export DART_INPUT_DIR=${DAT_DIR}/dart_input
#
# SHARE NAMELIST PARAMETERS:
export NL_WRF_CORE="'"ARW"'"
export NL_MAX_DOM=1
export NL_IO_FORM_GEOGRID=2
export NL_OPT_OUTPUT_FROM_GEOGRID_PATH="'"${RUN_DIR}/rc"'"
export NL_DEBUG_LEVEL=0
export VTABLE_TYPE=GFS
export METGRID_TABLE_TYPE=ARW
#
# GEOGRID NAMELIST PARAMETERS:
export NL_PARENT_ID=0
export NL_PARENT_GRID_RATIO=1
export NL_I_PARENT_START=1
export NL_J_PARENT_START=1
export NL_S_WE=1
export NL_E_WE=718
export NL_S_SN=1
export NL_E_SN=373
export NL_GEOG_DATA_RES="'"30s"'"
export NL_DX=15000
export NL_DY=15000
export NL_MAP_PROJ="'"lambert"'"
export NL_REF_LAT=30.095
export NL_REF_LON=-77.610
export NL_TRUELAT1=30.095
export NL_TRUELAT2=0.
export NL_STAND_LON=-77.610
export NL_GEOG_DATA_PATH="'"${WPS_GEOG_DIR}"'"
export NL_OPT_GEOGRID_TBL_PATH="'"${WPS_DIR}/geogrid"'"
#
# UNGRID NAMELIST PARAMETERS
export NL_OUT_FORMAT="'"WPS"'"
#
# METGRID NAMELIST PARAMETERS
export NL_IO_FORM_METGRID=2
export NL_OPT_IGNORE_DOM_CENTER=.false.
export NL_OPT_OUTPUT_FROM_METGRID_PATH="'"${RUN_DIR}/rc/2010091500"'"
#
# WRF NAMELIST PARAMETERS
#
# TIME CONTROL
export NL_RUN_DAYS=0
export NL_RUN_HOURS=${FCST_RANGE}
export NL_RUN_MINUTES=0
export NL_RUN_SECONDS=0
export NL_START_YEAR=2010
export NL_START_MONTH=09
export NL_START_DAY=15
export NL_START_HOUR=00
export NL_START_MINUTE=0
export NL_START_SECOND=0
export NL_END_YEAR=2010
export NL_END_MONTH=09
export NL_END_DAY=15
export NL_END_HOUR=0
export NL_END_MINUTE=0
export NL_END_SECOND=0
export NL_INTERVAL_SECONDS=${INTERVAL_SECONDS}
export NL_INPUT_FROM_FILE=".true."
export NL_HISTORY_INTERVAL=60
export NL_FRAMES_PER_OUTFILE=1
export NL_RESTART=".false."
export NL_RESTART_INTERVAL=2
export NL_IO_FORM_HISTORY=2
export NL_IO_FORM_RESTART=2
export NL_IO_FORM_INPUT=2
export NL_IO_FORM_BOUNDARY=2
export NL_WRITE_INPUT=".false."
export NL_INPUTOUT_INTERVAL=60
export NL_INPUT_OUTNAME="'"wrfout_gsi_d\<domain\>_\<date\>"'"
export NL_DEBUG_LEVEL=0
#
# DOMAINS
export NL_TIME_STEP=60
export NL_TIME_STEP_FRACT_NUM=0
export NL_TIME_STEP_FRACT_DEN=1
export NL_MAX_DOM=1
export NL_S_WE=1
export NL_E_WE=718
export NL_S_SN=1
export NL_E_SN=373
export NL_E_VERT=43
export NL_P_TOP_REQUESTED=3000
export NL_INTERP_TYPE=1
export NL_T_EXTRAP_TYPE=1
export NL_NUM_METGRID_LEVELS=27
export NL_NUM_METGRID_SOIL_LEVELS=4
export NL_DX=15000.0
export NL_DY=15000.0
export NL_GRID_ID=1
export NL_PARENT_ID=0
export NL_I_PARENT_START=1
export NL_J_PARENT_START=1
export NL_PARENT_GRID_RATIO=1
export NL_PARENT_TIME_STEP_RATIO=1
export NL_FEEDBACK=0
export NL_SMOOTH_OPTION=1
export NL_ETA_LEVELS=1.000,.9920,.9827,.9722,.9601,.9463,.9306,.9129,.8931,.8709,.8462,.8190,.7893,.7571,.7225,.6856,.6469,.6066,.5652,.5230,.4808,.4389,.3978,.3580,.3200,.2840,.2503,.2190,.1903,.1641,.1404,.1191,.1000,.0832,.0682,.0551,.0436,.0336,.0248,.0172,.0106,.0049,.0000
#
# PHYSICS
export NL_MP_PHYSICS=7
export NL_RA_LW_PHYSICS=1
export NL_RA_SW_PHYSICS=2
export NL_RADT=12
export NL_SF_SFCLAY_PHYSICS=1
export NL_SF_SURFACE_PHYSICS=2
export NL_BL_PBL_PHYSICS=1
export NL_BLDT=0
export NL_CU_PHYSICS=1
export NL_CUDT=5
export NL_ISFFLX=1
export NL_IFSNOW=0
export NL_ICLOUD=1
export NL_SURFACE_INPUT_SOURCE=1
export NL_NUM_SOIL_LAYERS=4
export NL_SF_URBAN_PHYSICS=1
export NL_MAXIENS=1
export NL_MAXENS=3
export NL_MAXENS2=3
export NL_MAXENS3=16
export NL_ENSDIM=144
#
# DYNAMICS
export NL_USE_BASEPARAM_FR_NML=".true."
export NL_W_DAMPING=1
export NL_DIFF_OPT=1
export NL_KM_OPT=4
export NL_DIFF_6TH_OPT=1
export NL_DIFF_6TH_FACTOR=0.12
export NL_BASE_TEMP=290.
export NL_DAMP_OPT=0
export NL_ZDAMP=5000
export NL_DAMPCOEF=0.2
export NL_KHDIF=0
export NL_KVDIF=0
export NL_NON_HYDROSTATIC=.true.
export NL_MOIST_ADV_OPT=1
export NL_SCALAR_ADV_OPT=1
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
export NL_NIO_TASKS_PER_GROUP=${NUM_PROCS}
export NL_NIO_GROUPS=${NUM_NODES}
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
export NL_CV_OPTIONS=6
#
# WRFVAR11
export NL_CV_OPTIONS_HUM=1
export NL_CHECK_RH=2
export NL_SEED_ARRAY1='set below'
export NL_SEED_ARRAY2='set=below'
export NL_CALCULATE_CG_COST_FN=true
export NL_LAT_STATS_OPTION=false
#
# WRFVAR16
export NL_ALPHA_METHOD=2
export NL_ENSDIM_ALPHA=0
export NL_ALPHA_TRUNCATION=1
export NL_ALPHA_CORR_TYPE=3
export NL_ALPHA_CORR_SCALE=${HOR_SCALE}
export NL_ALPHA_STD_DEV=1.0
#
# WRFVAR17
export NL_ANALYSIS_TYPE='set below'
#
# WRFVAR18
export NL_ANALYSIS_DATE='set below'
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
# SET TIME DATA
export DATE=${INITIAL_DATE}
. ${SCRIPTS_DIR}/da_get_date_range.ksh
#
# RUN GEOGRID
if [[ ${RUN_GEOGRID} = "true" ]]; then
mkdir -p ${RUN_DIR}/${DATE}/wps/geogrid
cd ${RUN_DIR}/${DATE}/wps/geogrid
#
# CREATE WPS NAMELIST
${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist.ksh
#
# RUN SCRIPT
ln -sf ${WPS_DIR}/geogrid.exe ./.
if [[ ${MACHINE} = "MAC" ]]; then
mpirun -np ${NUM_PROCS} ./geogrid.exe
fi
if [[ ${MACHINE} = "IBM" ]]; then
./geogrid.exe
fi
RC=$?
if [[ $RC != 0 ]]; then
echo geogrid failed with error $RC
exit $RC
fi
fi
#
# LOOP THROUGH DATES TO CREATE MET FILES
while [[ ${DATE} -le ${FINAL_DATE} ]] ; do 
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
export DATEPT=$(echo $DATE | cut -c1-8)
export HHPT=$(echo $DATE | cut -c9-10)
${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist.ksh
#
# CREATE UNGRIB FILE LIST
(( LBC_END=${FCST_RANGE} ))
integer ICNT=0
FILES=''
if [[ -e ${GRIB_DIR} ]]; then
 (( LBC_ITR=${LBC_START} ))
 while [[ ${LBC_ITR} -le ${LBC_END} ]]; do
    if [[ ${LBC_ITR} -lt 100  ]]; then export CMEM=0${LBC_ITR}; fi
    if [[ ${LBC_ITR} -lt 10   ]]; then export CMEM=00${LBC_ITR}; fi
    if [[ ${LBC_ITR} -eq 0    ]]; then export CMEM=000; fi
    export FILE_ZIP=${GRIB_DIR}/${DATEPT}_i${HHPT}_f${CMEM}_GFS004.grb2.gz
    cp ${FILE_ZIP} ${DATEPT}_i${HHPT}_f${CMEM}_GFS004.grb2.gz
    gunzip ${DATEPT}_i${HHPT}_f${CMEM}_GFS004.grb2.gz
    export FILE_UNZIP=${DATEPT}_i${HHPT}_f${CMEM}_GFS004.grb2
    FILES="${FILES} ${FILE_UNZIP}"
    ICNT=ICNT+1
    (( LBC_ITR=${LBC_ITR}+6 ))
 done
#
# LINK GRIB FILES
 ${WPS_DIR}/link_grib.csh $FILES
#
# RUN UNGRIB
 ln -fs ${WPS_DIR}/ungrib.exe ./.
 if [[ ${MACHINE} = "MAC" ]]; then
    ./ungrib.exe
 fi
 if [[ ${MACHINE} = "IBM" ]]; then
    ./ungrib.exe
 fi
 RC=$?
 if [[ $RC != 0 ]]; then
    echo ungrib failed with error $RC
    exit $RC
 fi
#
# REMOVE LINKS TO GFS FILES
 while [[ ${LBC_ITR} -le ${LBC_END} ]]; do
    if [[ ${LBC_ITR} -lt 100 ]]; then export CMEM=0${LBC_ITR}; fi
    if [[ ${LBC_ITR} -lt 10  ]]; then export CMEM=00${LBC_ITR}; fi
    if [[ ${LBC_ITR} -eq 0   ]]; then export CMEM=000; fi
    export FILE_UNZIP=${DATEPT}_i${HHPT}_f${CMEM}_GFS004.grb2
    rm ${FILE_ZIP} 
    rm ${FILE_UNZIP}
    ICNT=ICNT+1
    (( LBC_ITR=${LBC_ITR}+${LBC_FREQ} ))
 done
else
  echo ungrib failed because ${GRIB_DIR} does not exit
  export RC
  exit $RC
fi
fi 
#
# RUN METGRID
if [[ ${RUN_METGRID} = "true" ]]; then 
mkdir -p ${RUN_DIR}/${DATE}/wps/metgrid
cd ${RUN_DIR}/${DATE}/wps/metgrid
ln -fs ${FIXED_DATA_DIR}/geo_em.d01.nc ./.
ln -fs ../ungrib/namelist.wps ./.
ln -fs ../ungrib/FILE:* ./.
ln -fs $WPS_DIR/metgrid/METGRID.TBL.${METGRID_TABLE_TYPE} METGRID.TBL
ln -fs $WPS_DIR/metgrid.exe .
if [[ ${MACHINE} = "MAC" ]]; then
 mpirun -np ${NUM_PROCS} ./metgrid.exe
fi
if [[ ${MACHINE} = "IBM" ]]; then
 ./metgrid.exe
fi
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
 ln -sf ../../wps/metgrid/met_em* ./.
 ln -fs ${WRF_DIR}/main/real.exe ./.
 if [[ ${MACHINE} = "MAC" ]]; then
    mpirun -np ${NUM_PROCS} ./real.exe
 fi
 if [[ ${MACHINE} = "IBM" ]]; then
    ./real.exe
 fi
 RC=$?
 if [[ $RC != 0 ]]; then
    echo real failed with error $RC
    exit $RC
 fi
 if [[ ${L_DATE} = ${DATE_SAVE} ]]; then
    mv wrfinput_d01 wrfinput_d01.$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 0 -W 2>/dev/null)
    mv wrfbdy_d01 wrfbdy_d01.$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 0 -W 2>/dev/null)
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
       cp ${OB_PREPBUFR_DIR}/${DATE}/ob.bufr ob.bufr
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
    ln -sf ${BE_DIR}/be.dat be.dat
    ln -sf ${WRFVAR_DIR}/run/LANDUSE.TBL ./.
    ln -sf ${WRFVAR_DIR}/var/da/da_wrfvar.exe ./.
    . ${HYBRID_SCRIPTS_DIR}/da_create_wrfda_namelist.ksh
    if [[ ${MACHINE} = "MAC" ]]; then
       mpirun -np ${NUM_PROCS} ./da_wrfvar.exe
    fi
    if [[ ${MACHINE} = "IBM" ]]; then
       ./da_wrfvar.exe
    fi
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
               cp ${OB_PREPBUFR_DIR}/${DATE}/ob.bufr ob.bufr
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
            ln -fs ${DART_DIR}/models/wrf/work/pert_wrf_bc ./.
            . ${HYBRID_SCRIPTS_DIR}/da_create_wrfda_namelist.ksh
            if [[ ${MACHINE} = "MAC" ]]; then
               ./pert_wrf_bc
            fi
            if [[ ${MACHINE} = "IBM" ]]; then
               ./pert_wrf_bc
            fi
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
exit
