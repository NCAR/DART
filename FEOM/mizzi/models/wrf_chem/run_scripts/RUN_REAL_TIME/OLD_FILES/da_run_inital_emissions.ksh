#!/bin/ksh -aeux
#########################################################################
#
# Purpose: PERTURB WRFCHEM ICBC
#
#########################################################################
#
# CODE VERSIONS
export WPS_VER=WPSv3.4_dmpar
export WRFDA_VER=WRFDAv3.4_dmpar
export WRFDA_TOOLS_VER=WRFDA_TOOLSv3.4
export WRF_VER=WRFv3.4_dmpar
export WRFCHEM_VER=WRFCHEMv3.4_dmpar
export DART_VER=DART_CHEM
#
# EXPERIMENT DETAILS:
#export START_DATE=2008070100
export START_DATE=2008072400
export END_DATE=2008072723
export FCST_RANGE=6
export PROJ_NUMBER=P19010000
export PROJ_NUMBER=NACD0002
#
# EMISSIONS TO RUN
export RUN_EXO_COLDENS=true
export RUN_SEASON_WES=true
export RUN_WRF_BIO=true
export RUN_WRF_FIRE=true
#
# DIRECTORIES:
export SCRATCH_DIR=/glade/scratch/mizzi
export HOME_DIR=/glade/p/work/mizzi
export ACD_DIR=/glade/p/acd/mizzi
#
export TRUNK_DIR=${HOME_DIR}/TRUNK
export DART_DIR=${TRUNK_DIR}/${DART_VER}
export WRFVAR_DIR=${TRUNK_DIR}/${WRFDA_VER}
export BUILD_DIR=${WRFVAR_DIR}/var/da
export COLDENS_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_INITIAL_EMISS/WES-COLDENS
export MEGAN_BIO_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_INITIAL_EMISS/MEGAN-BIO
export MEGAN_DATA_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_INITIAL_EMISS/MEGAN-DATA
export WRF_FIRE_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_INITIAL_EMISS/WRF-FIRE
export RC_DIR=${ACD_DIR}/AVE_TEST_DATA/rc_100km
export WPB_RC_DIR=${ACD_DIR}/AVE_TEST_DATA/wpb_rc_100km
export RUN_DIR=${SCRATCH_DIR}/INITIAL_EMISSIONS
export ARCHIVE_DIR=${ACD_DIR}/AVE_TEST_DATA/chem_static_100km
mkdir -p ${ARCHIVE_DIR}
#
# SET DATE PARAMETERS
export DATE=${END_DATE}
export YY=$(echo $DATE | cut -c1-4)
export MM=$(echo $DATE | cut -c5-6)
export DD=$(echo $DATE | cut -c7-8)
export HH=$(echo $DATE | cut -c9-10)
export FIRE_END_DATE=${YY}-${MM}-${DD}
#
export DATE=${START_DATE}
export YY=$(echo $DATE | cut -c1-4)
export MM=$(echo $DATE | cut -c5-6)
export DD=$(echo $DATE | cut -c7-8)
export HH=$(echo $DATE | cut -c9-10)
export FILE_DATE=${YY}-${MM}-${DD}_${HH}:00:00
export DD_DATE=${YY}${MM}${DD}
export FIRE_START_DATE=${YY}-${MM}-${DD}
#
if [[ ! -d ${RUN_DIR} ]]; then
   mkdir ${RUN_DIR}
   cd ${RUN_DIR}
else
   cd ${RUN_DIR}
fi
#
# CREATE TIME INDEPENDENT EMISSION FILES
#
# BEGIN RUN EXO_COLDENS
if ${RUN_EXO_COLDENS}; then
   if [[ ! -d ${RUN_DIR}/EXO_COLDENS ]]; then
      mkdir ${RUN_DIR}/EXO_COLDENS
      cd ${RUN_DIR}/EXO_COLDENS
   else
      cd ${RUN_DIR}/EXO_COLDENS
   fi
#
# LINK NEEDED FILES
   export FILE=wrfinput_d01
   rm -rf ${FILE}
   ln -sf ${RC_DIR}/${DATE}/${FILE}_${FILE_DATE} ${FILE}   
   export FILE=exo_coldens.nc
   rm -rf ${FILE}
   ln -sf ${COLDENS_DIR}/${FILE} ${FILE}
   export FILE=exo_coldens
   rm -rf ${FILE}
   ln -sf ${COLDENS_DIR}/${FILE} ${FILE}
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
      cp ${FILE} ${ARCHIVE_DIR}/.
   fi
fi
#
# END RUN EXO_COLDENS
#
# BEGIN RUN SEASON_WES
if ${RUN_SEASON_WES}; then
   if [[ ! -d ${RUN_DIR}/SEASON_WES ]]; then
      mkdir ${RUN_DIR}/SEASON_WES
      cd ${RUN_DIR}/SEASON_WES
   else
      cd ${RUN_DIR}/SEASON_WES
   fi
#
# LINK NEEDED FILES
   export FILE=wrfinput_d01
   rm -rf ${FILE}
   ln -sf ${RC_DIR}/${DATE}/${FILE}_${FILE_DATE} ${FILE}   
   export FILE=season_wes_usgs.nc
   rm -rf ${FILE}
   ln -sf ${COLDENS_DIR}/${FILE} ${FILE}
   export FILE=wesely
   rm -rf ${FILE}
   ln -sf ${COLDENS_DIR}/${FILE} ${FILE}
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
      cp ${FILE} ${ARCHIVE_DIR}/.
   fi
fi
#
# END RUN SEASON_WES
#
# BEGIN RUN WRF_BIO
if ${RUN_WRF_BIO}; then
   if [[ ! -d ${RUN_DIR}/WRFBIO ]]; then
      mkdir ${RUN_DIR}/WRFBIO
      cd ${RUN_DIR}/WRFBIO
   else
      cd ${RUN_DIR}/WRFBIO
   fi
   export L_DATE=${START_DATE}
   while [[ ${L_DATE} -le ${END_DATE} ]]; do
      export YY=$(echo $L_DATE | cut -c1-4)
      export MM=$(echo $L_DATE | cut -c5-6)
      export DD=$(echo $L_DATE | cut -c7-8)
      export HH=$(echo $L_DATE | cut -c9-10)
      export L_FILE_DATE=${YY}-${MM}-${DD}_${HH}:00:00
#
# LINK NEEDED FILES
      export FILE=wrfinput_d01
      rm -rf ${FILE}
      ln -sf ${RC_DIR}/${L_DATE}/${FILE}_${L_FILE_DATE} ${FILE}   
      export FILE=wrfbiochemi_d01
      rm -rf ${FILE}
      rm -rf btr*.nc
      rm -rf DSW*.nc
      rm -rf hrb*.nc
      rm -rf iso*.nc
      rm -rf lai*.nc
      rm -rf ntr*.nc
      rm -rf shr*.nc
      rm -rf TAS*.nc
      ln -sf ${MEGAN_DATA_DIR}/*.nc ./.
      export FILE=megan_bio_emiss
      rm -rf ${FILE}
      ln -sf ${MEGAN_BIO_DIR}/${FILE} ${FILE}
#
# CREATE INPUT FILE
      export FILE=megan_bio_emiss.inp
      rm -rf ${FILE}
      cat << EOF > ${FILE}
&control
domains = 1,
start_lai_mnth = 5, 
end_lai_mnth = 7
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
#BSUB -q small
#
#
# RUN megan_bio_emis
./megan_bio_emiss < megan_bio_emiss.inp > index_megan_bio 2>&1
# 
export RC=\$?     
rm -rf WRF_BIO_SUCCESS
rm -rf WRF_BIO_FAILED          
if [[ \$RC = 0 ]]; then
   touch WRF_BIO_SUCCESS
else
   touch WRF_BIO_FAILED 
   exit
fi
EOF
#
# Run advance_model and wait until job completes
      bsub -K < job.ksh 
#### ./megan_bio_emiss < megan_bio_emiss.inp > index_megan_bio 2>&1
#
# TEST WHETHER OUTPUT EXISTS
      export FILE=wrfbiochemi_d01
      if [[ ! -e ${FILE} ]]; then
         echo WRFBIO FAILED
         exit
      else
         echo WRFBIO SUCCESS
         mv ${FILE} ${FILE}_${L_FILE_DATE}
         mkdir -p ${ARCHIVE_DIR}/${YY}${MM}${DD}
         cp ${FILE}_${L_FILE_DATE} ${ARCHIVE_DIR}/${YY}${MM}${DD}/.
      fi
      export P_DATE=${L_DATE}
      export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} 6 2>/dev/null)
   done
fi
#
# END RUN WRF_BIO
#
# BEGIN RUN WRF_FIRE
if ${RUN_WRF_FIRE}; then
   if [[ ! -d ${RUN_DIR}/WRFFIRE ]]; then
      mkdir ${RUN_DIR}/WRFFIRE
      cd ${RUN_DIR}/WRFFIRE
   else
      cd ${RUN_DIR}/WRFFIRE
   fi
#
# LINK NEEDED FILES
   export FILE=wrfinput_d01
   rm -rf ${FILE}
   ln -sf ${RC_DIR}/${DATE}/${FILE}_${FILE_DATE} ${FILE}   
   rm -rf GLOB_*.txt
   ln -sf ${WRF_FIRE_DIR}/GLOB_*.txt ./.
   export FILE=fire_emis
   rm -rf ${FILE}
   ln -sf ${WRF_FIRE_DIR}/src/${FILE} ${FILE}
   rm -rf grass_from_img.nc
   rm -rf shrub_from_img.nc
   rm -rf tempfor_from_img.nc
   rm -rf tropfor_from_img.nc
   ln -sf ${WRF_FIRE_DIR}/grass_from_img.nc
   ln -sf ${WRF_FIRE_DIR}/shrub_from_img.nc
   ln -sf ${WRF_FIRE_DIR}/tempfor_from_img.nc
   ln -sf ${WRF_FIRE_DIR}/tropfor_from_img.nc
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
# CREATE job.ksh
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
#
# RUN fire_emis
./fire_emis < fire_emis.mozc.inp > index_fire_emis 2>&1
# 
export RC=\$?     
rm -rf WRF_FIRE_SUCCESS     
rm -rf WRF_FIRE_FAILED          
if [[ \$RC = 0 ]]; then
   touch WRF_FIRE_SUCCESS
else
   touch WRF_FIRE_FAILED 
   exit
fi
EOF
#
# Run advance_model and wait until job completes
   bsub -K < job.ksh 
###./fire_emis < fire_emis.mozc.inp > index_fire_emis 2>&1
   export L_DATE=${START_DATE}
   while [[ ${L_DATE} -le ${END_DATE} ]]; do
      export YY=$(echo $L_DATE | cut -c1-4)
      export MM=$(echo $L_DATE | cut -c5-6)
      export DD=$(echo $L_DATE | cut -c7-8)
      export HH=$(echo $L_DATE | cut -c9-10)
      export L_FILE_DATE=${YY}-${MM}-${DD}_${HH}:00:00
      export DD_DATE=${YY}${MM}${DD}
#
# TEST WHETHER OUTPUT EXISTS
      export FILE=wrffirechemi_d01_${L_FILE_DATE}
      if [[ ! -e ${FILE} ]]; then
         echo WRFFIRE FAILED
         exit
      else
         echo WRFFIRE SUCCESS
         mkdir -p ${ARCHIVE_DIR}/${DD_DATE}
         cp ${FILE} ${ARCHIVE_DIR}/${DD_DATE}/.
      fi
      export P_DATE=${L_DATE}
      export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} 1 2>/dev/null)
   done
fi
#
# END RUN WRF_FIRE
