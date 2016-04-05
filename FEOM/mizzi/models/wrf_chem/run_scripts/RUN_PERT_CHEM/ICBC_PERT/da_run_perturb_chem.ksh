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
export START_DATE=2008063006
export START_DATE=2008070100
export START_DATE=2008071100
export END_DATE=2008073118
export END_DATE=2008072718
export NUM_MEMBERS=20
export MEM_START=1
export FCST_RANGE=6
export USE_HSI=false
#
#export MOZ_SPREAD=0.05
export MOZ_SPREAD=0.10
#export MOZ_SPREAD=0.20
#export MOZ_SPREAD=0.30
#export MOZ_SPREAD=0.35
#export MOZ_SPREAD=0.40
#export CSPREAD=p05
export CSPREAD=p10
#export CSPREAD=p20
#export CSPREAD=p30
#export CSPREAD=p35
#export CSPREAD=p40
#
# DIRECTORIES:
export SCRATCH_DIR=/glade/scratch/mizzi
export HOME_DIR=/glade/p/work/mizzi
export ACD_DIR=/glade/p/acd/mizzi
#
export TRUNK_DIR=${HOME_DIR}/TRUNK
export DART_DIR=${TRUNK_DIR}/${DART_VER}
export DART_PERT_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_PERT_CHEM/ICBC_PERT
export WRFVAR_DIR=${TRUNK_DIR}/${WRFDA_VER}
export BUILD_DIR=${WRFVAR_DIR}/var/da
export MOZBC_DIR=${HOME_DIR}/mozbc-ens
export RUN_DIR=${SCRATCH_DIR}/PERTURB_ICBC_${CSPREAD}
export RC_CHEM_DIR=${ACD_DIR}/AVE_TEST_DATA/rc_chem_100km_${CSPREAD}
export HSI_RC_CHEM_DIR=/MIZZI/AVE_TEST_DATA/rc_chem_100km_${CSPREAD}
export WPB_RC_CHEM_DIR=${ACD_DIR}/AVE_TEST_DATA/wpb_rc_chem_100km_${CSPREAD}
export HSI_WPB_RC_CHEM_DIR=/MIZZI/AVE_TEST_DATA/wpb_rc_chem_100km_${CSPREAD}
#
export MET_EM_DIR=${ACD_DIR}/AVE_TEST_DATA/met_em_100km
export HSI_MET_EM_DIR=/MIZZI/AVE_TEST_DATA/met_em_100km
#
export RC_DIR=${ACD_DIR}/AVE_TEST_DATA/rc_100km
export HSI_RC_DIR=/MIZZI/AVE_TEST_DATA/rc_100km
export WPB_RC_DIR=${ACD_DIR}/AVE_TEST_DATA/wpb_rc_100km
export HSI_WPB_RC_DIR=/MIZZI/AVE_TEST_DATA/wpb_rc_100km
#
# LOOP THROUGH DATES
export DATE=${START_DATE}
export YYYY=$(echo $DATE | cut -c1-4)
export MM=$(echo $DATE | cut -c5-6)
export DD=$(echo $DATE | cut -c7-8)
export HH=$(echo $DATE | cut -c9-10)
export NEXT_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} ${FCST_RANGE} 2>/dev/null)
export NEXT_YYYY=$(echo $NEXT_DATE | cut -c1-4)
export NEXT_MM=$(echo $NEXT_DATE | cut -c5-6)
export NEXT_DD=$(echo $NEXT_DATE | cut -c7-8)
export NEXT_HH=$(echo $NEXT_DATE | cut -c9-10)
while [[ ${DATE} -le ${END_DATE} ]] ; do
   if [[ ! -d ${RUN_DIR}/${DATE} ]]; then
      mkdir -p ${RUN_DIR}/${DATE}
      cd ${RUN_DIR}/${DATE}
   else
      cd ${RUN_DIR}/${DATE}
   fi
   mkdir -p ${RC_CHEM_DIR}/${DATE} 
   mkdir -p ${WPB_RC_CHEM_DIR}/${DATE}  
   if ${USE_HSI}; then
      hsi mkdir -p ${HSI_RC_CHEM_DIR}/${DATE} 
      hsi mkdir -p ${HSI_WPB_RC_CHEM_DIR}/${DATE}  
   fi
#
# PERTURB CHEM ICBC
   cp ${DART_PERT_DIR}/runICBC_parent.ksh ./.
   cp ${DART_PERT_DIR}/runICBC_setN.ksh ./.
#   cp ${DART_PERT_DIR}/random_log_normal.py ./.
   cp ${DART_PERT_DIR}/random.py ./.
   cp ${DART_PERT_DIR}/run_mozbc.csh ./.
   cp ${DART_PERT_DIR}/mozbc-dart/mozbc ./.
#   cp ${DART_PERT_DIR}/ln_set0 .
   cp ${DART_PERT_DIR}/set0 ./.
   cp ${DART_PERT_DIR}/set00 ./.
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
dir_wrf  = '${RUN_DIR}/${DATE}/'
dir_moz = '${DART_PERT_DIR}/MOZART/'
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
dir_wrf  = '${RUN_DIR}/${DATE}/'
dir_moz = '${DART_PERT_DIR}/MOZART/'
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
dir_wrf  = '${RUN_DIR}/${DATE}/'
dir_moz = '${DART_PERT_DIR}/MOZART/'
fn_moz  = '${MOZART_DATA}'
def_missing_var = .true.
met_file_prefix  = 'met_em'
met_file_suffix  = '.nc'
met_file_separator= '.'
EOF
#
# APM: This needs to be fixed.  The number of met_em files depends on the forecast period.  The hsi get does not 
#      work.  The cp copies all met_em files in the current directory (assuming there are sufficent times there). 
   if ${USE_HSI}; then
      hsi get ${HSI_MET_EM_DIR}/${DATE}/met_em.d01.${YYYY}-${MM}-${DD}_${HH}:00:00.nc
      hsi get ${HSI_MET_EM_DIR}/${NEXT_DATE}/met_em.d01.${NEXT_YYYY}-${NEXT_MM}-${NEXT_DD}_${NEXT_HH}:00:00.nc
   else
#      cp ${MET_EM_DIR}/${DATE}/met_em.d01.${YYYY}-${MM}-${DD}_${HH}:00:00.nc ./.
      cp ${MET_EM_DIR}/${DATE}/met_em.d01.*:00:00.nc ./.
#      cp ${MET_EM_DIR}/${NEXT_DATE}/met_em.d01.${NEXT_YYYY}-${NEXT_MM}-${NEXT_DD}_${NEXT_HH}:00:00.nc ./.
   fi
#
#   ./random_log_normal.py ${MOZ_SPREAD} ${NUM_MEMBERS} ${DART_PERT_DIR} ${RUN_DIR}/${DATE}
   ./random.py ${MOZ_SPREAD} ${NUM_MEMBERS} ${DART_PERT_DIR} ${RUN_DIR}/${DATE}
   ./runICBC_parent.ksh
   ./runICBC_setN.ksh
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
         if ${USE_HSI}; then
            hsi get wrfinput.met : ${HSI_WPB_RC_DIR}/${DATE}/${WRFINPEN}
            hsi get wrfbdy.met : ${HSI_WPB_RC_DIR}/${DATE}/${WRFBDYEN}
         else
            cp ${WPB_RC_DIR}/${DATE}/${WRFINPEN} wrfinput.met
            cp ${WPB_RC_DIR}/${DATE}/${WRFBDYEN} wrfbdy.met
         fi
         ncks -x -v MU,PSFC,Q2,T2,TH2,U10,V10,P,QVAPOR,T,U,V ${WPB_RC_CHEM_DIR}/${DATE}/${WRFINPEN} wrfinput_d01
         ncks -v MU,PSFC,Q2,T2,TH2,U10,V10,P,QVAPOR,T,U,V wrfinput.met ${WRFINPEN}
         ncks -A wrfinput_d01 ${WRFINPEN}
#
         ncks -x -v MU_BTXE,MU_BTXS,MU_BTYE,MU_BTYS,MU_BXE,MU_BXS,MU_BYE,MU_BYS,QVAPOR_BTXE,QVAPOR_BTXS,QVAPOR_BTYE,QVAPOR_BTYS,QVAPOR_BXE,QVAPOR_BXS,QVAPOR_BYE,QVAPOR_BYS,T_BTXE,T_BTXS,T_BTYE,T_BTYS,T_BXE,T_BXS,T_BYE,T_BYS,U_BTXE,U_BTXS,U_BTYE,U_BTYS,U_BXE,U_BXS,U_BYE,U_BYS,V_BTXE,V_BTXS,V_BTYE,V_BTYS,V_BXE,V_BXS,V_BYE,V_BYS,PH_BTXE,PH_BTXS,PH_BTYE,PH_BTYS,PH_BXE,PH_BXS,PH_BYE,PH_BYS ${WPB_RC_CHEM_DIR}/${DATE}/${WRFBDYEN} wrfbdy_d01
         ncks -v MU_BTXE,MU_BTXS,MU_BTYE,MU_BTYS,MU_BXE,MU_BXS,MU_BYE,MU_BYS,QVAPOR_BTXE,QVAPOR_BTXS,QVAPOR_BTYE,QVAPOR_BTYS,QVAPOR_BXE,QVAPOR_BXS,QVAPOR_BYE,QVAPOR_BYS,T_BTXE,T_BTXS,T_BTYE,T_BTYS,T_BXE,T_BXS,T_BYE,T_BYS,U_BTXE,U_BTXS,U_BTYE,U_BTYS,U_BXE,U_BXS,U_BYE,U_BYS,V_BTXE,V_BTXS,V_BTYE,V_BTYS,V_BXE,V_BXS,V_BYE,V_BYS,PH_BTXE,PH_BTXS,PH_BTYE,PH_BTYS,PH_BXE,PH_BXS,PH_BYE,PH_BYS wrfbdy.met ${WRFBDYEN}
         ncks -A wrfbdy_d01 ${WRFBDYEN}
         rm -rf wrfinput.met
         rm -rf wrfbdy.met
#
# ARCHIVE MET/CHEM FILES
         export INP_FILE=${WPB_RC_CHEM_DIR}/${DATE}/${WRFINPEN}
         export BDY_FILE=${WPB_RC_CHEM_DIR}/${DATE}/${WRFBDYEN}
         if [[ -e ${INP_FILE} ]]; then rm -rf ${INP_FILE}; fi
         if [[ -e ${BDY_FILE} ]]; then rm -rf ${BDY_FILE}; fi
         cp ${WRFINPEN} ${INP_FILE}
         cp ${WRFBDYEN} ${BDY_FILE}
         let MEM=MEM+1
      done
#
# ADVANCE TIME
      export DATE=${NEXT_DATE}      
      export YYYY=$(echo $DATE | cut -c1-4)
      export MM=$(echo $DATE | cut -c5-6)
      export DD=$(echo $DATE | cut -c7-8)
      export HH=$(echo $DATE | cut -c9-10)
      export NEXT_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} ${FCST_RANGE} 2>/dev/null)
      export NEXT_YYYY=$(echo $NEXT_DATE | cut -c1-4)
      export NEXT_MM=$(echo $NEXT_DATE | cut -c5-6)
      export NEXT_DD=$(echo $NEXT_DATE | cut -c7-8)
      export NEXT_HH=$(echo $NEXT_DATE | cut -c9-10)
done
