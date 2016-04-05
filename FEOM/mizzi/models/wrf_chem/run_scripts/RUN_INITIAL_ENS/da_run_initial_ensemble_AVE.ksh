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
export RUN_GEOGRID=true
export RUN_UNGRIB=true
export RUN_METGRID=true
export RUN_REAL=true
export RUN_RANDOMCV=true
export RUN_PERT_WRF_BC=true
export RUN_ARCHIVE=true
#
# CODE VERSIONS
export WPS_VER=WPSv3.4_dmpar
export WPS_GEOG_VER=WPSv3.4_GEOG_DATA
export WRFDA_VER=WRFDAv3.4_dmpar
export WRFDA_TOOLS_VER=WRFDA_TOOLSv3.4
export WRF_VER=WRFv3.4_dmpar
export DART_VER=DART_CHEM
#
# Set job submission parameters
#export PROJ_NUMBER=P19010000
export PROJ_NUMBER=NACD0002
export TIME_LIMIT=0:10
export NUM_TASKS=32
export TASKS_PER_NODE=8
export JOB_CLASS=small
export SINGLE_FILE=false
#
# EXPERIMENT DETAILS:
#export INITIAL_DATE=2008061800
export INITIAL_DATE=2008072806
export FINAL_DATE=2008073118
export NUM_MEMBERS=20
export NUM_PROCS=8
export MEM_START=1
export REGION=dart_test
export NNXP=101
export NNYP=41
# NAM forecast data
#export FG_TYPE=NAM
#export GRIB_PART1=nam_218_
#export GRIB_PART2=.grb2
# GFS forecast data
export FG_TYPE=GFS

export GRIB_PART1=gfs_4_
export GRIB_PART2=.g2.tar
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
export FCST_RANGE=24
export FCST_RANGE=6
# APM: due to problem with the 06-12-18 GFS forecast, FCST_RANGE=6 for that date only
#export FCST_RANGE=6
(( LBC_END=2*${FCST_RANGE} ))
export ASIM_WINDOW=3
export LBC_FREQ=3
(( INTERVAL_SECONDS=${LBC_FREQ}*60*60 ))
export LBC_START=0
#
# DIRECTORIES:
export RUN_DIR=/glade/scratch/mizzi/${REGION}
export TRUNK_DIR=/glade/p/work/mizzi/TRUNK
export HYBRID_TRUNK_DIR=/glade/p/work/mizzi/HYBRID_TRUNK
export DAT_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA
#
export HYBRID_SCRIPTS_DIR=${HYBRID_TRUNK_DIR}/hybrid_scripts
export GRIB_DIR=${DAT_DIR}/NAM_DATA_NOMADS
if [[ ${FG_TYPE} == GFS ]]; then export GRIB_DIR=${DAT_DIR}/GFS_DATA_NOMADS; fi
export SCRIPTS_DIR=${TRUNK_DIR}/${WRFDA_TOOLS_VER}/scripts
export WPS_DIR=${TRUNK_DIR}/${WPS_VER}
export WPS_GEOG_DIR=${TRUNK_DIR}/${WPS_GEOG_VER}/geog
export WRFVAR_DIR=${TRUNK_DIR}/${WRFDA_VER}
export DART_DIR=${TRUNK_DIR}/${DART_VER}
export BUILD_DIR=${WRFVAR_DIR}/var/da
export WRF_DIR=${TRUNK_DIR}/${WRF_VER}
export OBSPROC_DIR=${WRFVAR_DIR}/var/obsproc
export VTABLE_DIR=${WPS_DIR}/ungrib/Variable_Tables
export OB_PREPBUFR_DIR=${DAT_DIR}/obs_MET
export OB_ASCII_DIR=${DAT_DIR}/obs_MET
export BE_DIR=${WRFVAR_DIR}/var/run
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
# APM: fix this path setting
export NL_GEOG_DATA_PATH="'"${WPS_GEOG_DIR}"'"
export NL_OPT_GEOGRID_TBL_PATH="'"${WPS_DIR}/geogrid"'"
#
# UNGRIB NAMELIST PARAMETERS
export NL_OUT_FORMAT="'"WPS"'"
#
# METGRID NAMELIST PARAMETERS
export NL_IO_FORM_METGRID=2
export NL_OPT_IGNORE_DOM_CENTER=.false.
# APM: problem with this setting
export NL_OPT_OUTPUT_FROM_METGRID_PATH="'"${RUN_DIR}/rc"'"
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
export NL_NUM_METGRID_LEVELS=40
if [[ ${FG_TYPE} == GFS ]]; then export NL_NUM_METGRID_LEVELS=27; fi
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
# PHYSICS
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
# DYNAMICS
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
export GEOGRID_PATH=${RUN_DIR}/geogrid
#
# RUN GEOGRID
if [[ ${RUN_GEOGRID} = "true" ]]; then
   mkdir -p ${RUN_DIR}/geogrid
   cd ${RUN_DIR}/geogrid
#
# CREATE WPS NAMELIST
   ${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist.ksh
#
# RUN SCRIPT
   ln -sf ${WPS_DIR}/geogrid.exe ./.
#
# APM: Create job script 
   if [[ -f job.ksh ]]; then rm -rf job.ksh; fi
   touch job.ksh
   RANDOM=$$
   export JOBRND=geogrid_$RANDOM
   rm -rf *.jerr
   rm -rf *.jout
   cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -n 1                                  # number of total (MPI) tasks
#BSUB -R "span[ptile=${TASKS_PER_NODE}]"    # mpi tasks per node
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.jout                      # output filename
#BSUB -e ${JOBRND}.jerr                      # error filename
#BSUB -W ${TIME_LIMIT}                      # wallclock time (minutes)
#BSUB -q geyser 
#
# Run geogrid
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
#
# Submit script and wait until job completes
   bsub -K < job.ksh 
#   mpirun -np ${NUM_PROCS} ./geogrid.exe
fi
#
# LOOP THROUGH DATES TO CREATE INPUT FILES
while [[ ${DATE} -le ${FINAL_DATE} ]] ; do 
   export YYYY=$(echo $DATE | cut -c1-4)
   export MM=$(echo $DATE | cut -c5-6)
   export DD=$(echo $DATE | cut -c7-8)
   export HH=$(echo $DATE | cut -c9-10)
#
# CREATE UNGRIB WORKING DIRECTORY
   if [[ ${RUN_UNGRIB} = "true" ]]; then 
      mkdir -p ${RUN_DIR}/${DATE}/ungrib
      cd ${RUN_DIR}/${DATE}/ungrib
      rm -rf GRIBFILE.*
      export FCST_RANGE_SAVE=${FCST_RANGE}
      export FCST_RANGE=${LBC_END}
      . ${SCRIPTS_DIR}/da_get_date_range.ksh
      export NL_START_DATE="'"${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:00:00"'"
      export NL_END_DATE="'"${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:00:00"'"
      ${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist.ksh
#
# LINK VTABLE
      ln -fs ${VTABLE_DIR}/Vtable.${VTABLE_TYPE} Vtable
#
# CREATE UNGRIB FILE LIST
      FILES=''
      if [[ -e ${GRIB_DIR}/${DATE} ]]; then
         if [[ -e ${GRIB_DIR}/${DATE}/${GRIB_PART1}${DATE}${GRIB_PART2} ]]; then
            cd ${GRIB_DIR}/${DATE}
            tar -xf ${GRIB_DIR}/${DATE}/${GRIB_PART1}${DATE}${GRIB_PART2}
            cd ${RUN_DIR}/${DATE}/ungrib
         else
            echo 'APM: No metfiles in directory'
            exit
         fi
# 
# Forecasts in separate GRIB files 
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
      cp ${WPS_DIR}/ungrib.exe ./.
      ./ungrib.exe
      RC=$?
      if [[ $RC != 0 ]]; then
         echo ungrib failed with error $RC
         exit $RC
      fi
#
# TAR THE PARENT FORECAST FILES (special fix for 2008061218 becaseu the 033 and 048 forecasts are missing)
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
      export FCST_RANGE=${FCST_RANGE_SAVE}
      . ${SCRIPTS_DIR}/da_get_date_range.ksh
      export NL_START_DATE="'"${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:00:00"'"
      export NL_END_DATE="'"${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:00:00"'"
   fi
#
# RUN METGRID
   if [[ ${RUN_METGRID} = "true" ]]; then 
      mkdir -p ${RUN_DIR}/${DATE}/metgrid
      cd ${RUN_DIR}/${DATE}/metgrid
      ln -fs ${GEOGRID_PATH}/geo_em.d01.nc ./.
#      ln -fs ../ungrib/namelist.wps ./.
      ln -fs ../ungrib/FILE:* ./.
      ln -fs $WPS_DIR/metgrid/METGRID.TBL.${METGRID_TABLE_TYPE} METGRID.TBL
      ln -fs $WPS_DIR/metgrid.exe .
      export FCST_RANGE_SAVE=${FCST_RANGE}
      export FCST_RANGE=${LBC_END}
      . ${SCRIPTS_DIR}/da_get_date_range.ksh
      export NL_START_DATE="'"${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:00:00"'"
      export NL_END_DATE="'"${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:00:00"'"
      ${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist.ksh
#
# APM: Create job script 
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
#      mpirun -np ${NUM_PROCS} ./metgrid.exe
      export FCST_RANGE=${FCST_RANGE_SAVE}
      . ${SCRIPTS_DIR}/da_get_date_range.ksh
      export NL_START_DATE="'"${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:00:00"'"
      export NL_END_DATE="'"${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:00:00"'"
   fi
#
# RUN REAL
   if [[ ${RUN_REAL} = "true" ]]; then 
      mkdir -p ${RUN_DIR}/${DATE}/real
      cd ${RUN_DIR}/${DATE}/real
      export DATE_SAV=${DATE}
#
# LINK IN THE METGRID FILES
      export L_DATE=${DATE}
      export L_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} ${LBC_END} 2>/dev/null)
      while [[ ${L_DATE} -le ${L_END_DATE} ]] ; do
         export YY=$(echo $L_DATE | cut -c1-4)
         export MM=$(echo $L_DATE | cut -c5-6)
         export DD=$(echo $L_DATE | cut -c7-8)
         export HH=$(echo $L_DATE | cut -c9-10)
         export L_FILE_DATE=${YY}-${MM}-${DD}_${HH}:00:00.nc
         ln -sf ${RUN_DIR}/${DATE}/metgrid/met_em.d01.${L_FILE_DATE} ./.
         export Q_DATE=${L_DATE}
         export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${Q_DATE} ${LBC_FREQ} 2>/dev/null) 
      done
      ln -fs ${WRF_DIR}/main/real.exe ./.
#
# LOOP THROUGH BDY TENDENCY TIMES FOR PERTURB_BC
      export L_DATE=${DATE}
      export L_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} ${FCST_RANGE} 2>/dev/null)
      while [[ ${L_DATE} -le ${L_END_DATE} ]] ; do      
#
# CREATE WRF NAMELIST
         export DATE=${L_DATE}
         . ${SCRIPTS_DIR}/da_get_date_range.ksh
         . ${HYBRID_SCRIPTS_DIR}/da_create_wrf_namelist.ksh
#
# APM: Create job script 
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
# Submit convert file script for each and wait until job completes
         bsub -K < job.ksh 
         mv wrfinput_d01 wrfinput_d01.$(${BUILD_DIR}/da_advance_time.exe ${DATE} 0 -W 2>/dev/null)
         mv wrfbdy_d01 wrfbdy_d01.$(${BUILD_DIR}/da_advance_time.exe ${DATE} 0 -W 2>/dev/null)
         export Q_DATE=${L_DATE}
         export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${Q_DATE} ${LBC_FREQ} 2>/dev/null) 
      done   
      export DATE=${DATE_SAV}
   fi
#
# RUN RANDOMCV
   if [[ ${RUN_RANDOMCV} = "true" ]]; then 
#
# GO TO RUN DIRECTORY
      if [[ ! -d ${RUN_DIR}/${DATE}/wrfda_ic ]]; then
         mkdir -p ${RUN_DIR}/${DATE}/wrfda_ic
         cd ${RUN_DIR}/${DATE}/wrfda_ic
      else
         cd ${RUN_DIR}/${DATE}/wrfda_ic
      fi
#
# SET OBS PATH (DO I NEED OBS FOR RANDOMCV?  THESE OBS MAY BE AT THE WRONG TIME)
      if [[ ${USE_PREPBUFR} = "true" ]]; then
         cp ${OB_PREPBUFR_DIR}/${DATE}/prepbufr.gdas.${DATE}.wo40.be ob.bufr
         export NL_OB_FORMAT=1
      else
         cp ${OB_ASCII_DIR}/${DATE}/ob.ascii ob.ascii
         export NL_OB_FORMAT=2
      fi
#
# LOOP THROUGH ALL BDY TENDENCY TIMES
      export L_DATE=${DATE}
      export L_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} ${FCST_RANGE} 2>/dev/null)
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
         ln -sf ${DA_INPUT_FILE} fg
         ln -sf ${BE_DIR}/be.dat.cv3 be.dat
         ln -sf ${WRFVAR_DIR}/run/LANDUSE.TBL ./.
         ln -sf ${WRFVAR_DIR}/var/da/da_wrfvar.exe ./.
#
# LOOP THROUGH ALL MEMBERS IN THE ENSEMBLE
         let MEM=${MEM_START}
         while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
            export CMEM=e${MEM}
            if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
            if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
            export NL_SEED_ARRAY1=$(${BUILD_DIR}/da_advance_time.exe ${DATE} 0 -f hhddmmyycc)
            export NL_SEED_ARRAY2=`echo ${MEM} \* 100000 | bc -l `
            . ${HYBRID_SCRIPTS_DIR}/da_create_wrfda_namelist.ksh
#
# APM: Create job script 
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
# Run real
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
# Submit convert file script for each and wait until job completes
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
# RUN PERTURB_WRF_BC
   if [[ ${RUN_PERT_WRF_BC} = "true" ]]; then 
      if [[ ! -d ${RUN_DIR}/${DATE}/wrfda_bc ]]; then
         mkdir -p ${RUN_DIR}/${DATE}/wrfda_bc
         cd ${RUN_DIR}/${DATE}/wrfda_bc
      else
         cd ${RUN_DIR}/${DATE}/wrfda_bc
      fi
#
# LOOP THROUGH ALL MEMBERS IN THE ENSEMBLE
      let MEM=${MEM_START}
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
         export L_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} ${FCST_RANGE} 2>/dev/null)
         while [[ ${L_DATE} -lt ${L_END_DATE} ]]; do
            export ANALYSIS_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 0 -W 2>/dev/null)
            export NEXT_L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${LBC_FREQ} 2>/dev/null)
            export NEXT_ANALYSIS_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${LBC_FREQ} -W 2>/dev/null)
            rm -rf wrfinput_this
            rm -rf wrfinput_next
            export DA_INPUT_PATH=${RUN_DIR}/${DATE}/wrfda_ic
            ln -fs ${DA_INPUT_PATH}/wrfinput_d01.${ANALYSIS_DATE}.${CMEM} wrfinput_this 
            ln -fs ${DA_INPUT_PATH}/wrfinput_d01.${NEXT_ANALYSIS_DATE}.${CMEM} wrfinput_next 
#
# APM: Create job script 
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
# Run real
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
# Submit convert file script for each and wait until job completes
            bsub -K < job.ksh 
#
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
# RUN ARCHIVE
   if [[ ${RUN_ARCHIVE} = "true" ]]; then 
# 
# COPY ENSEMBLE DATA TO WPB_RC
      export WPB_RC_DIR=${RUN_DIR}/wpb_rc/${DATE}
      export RC_DIR=${RUN_DIR}/rc/${DATE}
      export MET_DIR=${RUN_DIR}/met_em/${DATE}
      mkdir -p ${WPB_RC_DIR}
      mkdir -p ${RC_DIR}
      mkdir -p ${MET_DIR}
      export ANALYSIS_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} 0 -W 2>/dev/null)
      cp ${RUN_DIR}/${DATE}/real/wrfinput_d01.${ANALYSIS_DATE} ${RC_DIR}/wrfinput_d01_${ANALYSIS_DATE} 
      cp ${RUN_DIR}/${DATE}/real/wrfbdy_d01.${ANALYSIS_DATE} ${RC_DIR}/wrfbdy_d01_${ANALYSIS_DATE}
      cp ${RUN_DIR}/${DATE}/metgrid/met_em* ${MET_DIR}
      let MEM=${MEM_START}
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
         cp ${RUN_DIR}/${DATE}/wrfda_ic/wrfinput_d01.${ANALYSIS_DATE}.${CMEM} ${WPB_RC_DIR}/wrfinput_d01_${ANALYSIS_DATE}.${CMEM} 
         cp ${RUN_DIR}/${DATE}/wrfda_bc/wrfbdy_d01.${ANALYSIS_DATE}.${CMEM} ${WPB_RC_DIR}/wrfbdy_d01_${ANALYSIS_DATE}.${CMEM} 
         
         let MEM=${MEM}+1
      done
#
# REMOVE THE CURRENT RUN DIRECTORY
      cd ${RUN_DIR}
      rm -rf ${DATE} 
   fi
   export NEXT_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} ${CYCLE_PERIOD} 2>/dev/null)
   export DATE=${NEXT_DATE}
done

