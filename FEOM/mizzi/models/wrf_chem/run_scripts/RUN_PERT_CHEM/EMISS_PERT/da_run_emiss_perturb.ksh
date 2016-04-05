#!/bin/ksh -aeux
#########################################################################
#
# Purpose: EMISS PERTURB
#
#########################################################################
#
# CODE VERSIONS
export WRFDA_VER=WRFDAv3.4_dmpar
export WRFDA_TOOLS_VER=WRFDA_TOOLSv3.4
export WRF_VER=WRFv3.6.1_dmpar
export WRFCHEM_VER=WRFCHEMv3.6.1_dmpar
export DART_VER=DART_CHEM_MY_BRANCH
#
# EXPERIMENT DETAILS:
export NUM_MEMBERS=20
export MEM_START=1
export USE_HSI=false
#
# TIME DATA:
export START_DATE=2008060100
export START_DATE=2008062406
export END_DATE=2008062400
export END_DATE=2008063018
export DATE=${START_DATE}
export YYYY=$(echo $DATE | cut -c1-4)
export MM=$(echo $DATE | cut -c5-6)
export DD=$(echo $DATE | cut -c7-8)
export HH=$(echo $DATE | cut -c9-10)
#
# CYCLING AND BOUNDARY CONDITION TIME DATA
export EMISS_FREQ=1
export EMISS_FREQ=6
export NL_MEAN=1.0
#export NL_SPREAD=0.05
#export NL_SPREAD=0.20
#export NL_SPREAD=0.35
#export NL_SPREAD=0.40
#
export NL_SPREAD=0.00
#export NL_SPREAD=0.30
export CSPREAD=p00
#export CSPREAD=p30
#
# DIRECTORIES:
export SCRATCH_DIR=/glade/scratch/mizzi
export HOME_DIR=/glade/p/work/mizzi
export ACD_DIR=/glade/p/acd/mizzi
#
export RUN_DIR=${SCRATCH_DIR}/PERTURB_EMISS_${CSPREAD}
export SAVE_DIR=${ACD_DIR}/AVE_TEST_DATA/chem_static_100km_${CSPREAD}
export TRUNK_DIR=${HOME_DIR}/TRUNK
export WRFVAR_DIR=${TRUNK_DIR}/${WRFDA_VER}
export BUILD_DIR=${WRFVAR_DIR}/var/da
export WRF_DIR=${TRUNK_DIR}${WRF_VER}
export WRFCHEM_DIR=${TRUNK_DIR}/${WRFCHEM_VER}
export DATA_DIR=${ACD_DIR}/AVE_TEST_DATA/chem_static_100km
export HSI_DATA_DIR=/MIZZI/AVE_TEST_DATA/chem_static_100km
export DART_DIR=${TRUNK_DIR}/${DART_VER}/models/wrf_chem/run_scripts/RUN_PERT_CHEM/EMISS_PERT
#
# SAVE FILES TO ${SAVE_DIR}
if [[ ! -e ${SAVE_DIR} ]]; then
   mkdir -p ${SAVE_DIR}
fi
cp ${DATA_DIR}/clim_p_trop.nc ${SAVE_DIR}/.
cp ${DATA_DIR}/exo_coldens_d01 ${SAVE_DIR}/.
cp ${DATA_DIR}/ubvals_b40.20th.track1_1996-2005.nc ${SAVE_DIR}/.
cp ${DATA_DIR}/wrf_season_wes_usgs_d01.nc ${SAVE_DIR}/.
#
# LOOP THROUGH DATES
while [[ ${DATE} -le ${END_DATE} ]] ; do
   if [[ ! -d ${RUN_DIR}/${YYYY}${MM}${DD} ]]; then
      mkdir -p ${RUN_DIR}/${YYYY}${MM}${DD}
      cd ${RUN_DIR}/${YYYY}${MM}${DD}
   else
      cd ${RUN_DIR}/${YYYY}${MM}${DD}
   fi
   if [[ ! -e ${SAVE_DIR}/${YYYY}${MM}${DD} ]]; then
      mkdir -p ${SAVE_DIR}/${YYYY}${MM}${DD}
   fi
#
# COPY PERTURBATION CODE
   if [[ -e perturb_chem_emiss.exe ]]; then rm -rf perturb_chem_emiss.exe; fi
   cp ${DART_DIR}/perturb_chem_emiss.exe ./.
#
# LOOP THROUGH ALL MEMBERS IN THE ENSEMBLE
   let MEM=${MEM_START}
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export NL_ENS_MEMBER=${MEM}
         export NL_PERT_CHEM=false
         export NL_PERT_FIRE=false
         export NL_PERT_BIO=true
         if [[ ${HH} -eq 00 || ${HH} -eq 06 || ${HH} -eq 12 || ${HH} -eq 18 ]]; then
            export NL_PERT_BIO=true
         fi
         export CMEM=e${MEM}
         if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
#
# GET EMISSON FILES FOR THIS MEMBER
         export WRFCHEMI=wrfchemi_d01_${YYYY}-${MM}-${DD}_${HH}:00:00
         export WRFFIRECHEMI=wrffirechemi_d01_${YYYY}-${MM}-${DD}_${HH}:00:00
         export WRFBIOCHEMI=wrfbiochemi_d01_${YYYY}-${MM}-${DD}_${HH}:00:00
#
         if [[ ${NL_PERT_CHEM} == true ]]; then
            if [[ ! -e ${WRFCHEMI} ]]; then
               if [[ ${USE_HSI} == true ]]; then 
                  hsi -d4 get ${HSI_DATA_DIR}/${YYYY}${MM}${DD}/${WRFCHEMI}
               else
                  cp ${DATA_DIR}/${YYYY}${MM}${DD}/${WRFCHEMI} ./.
               fi
               cp ${WRFCHEMI} ${WRFCHEMI}.${CMEM}
            else
               cp ${WRFCHEMI} ${WRFCHEMI}.${CMEM}
            fi
            mv ${WRFCHEMI} ${SAVE_DIR}/${YYYY}${MM}${DD}/
         fi
         if [[ ${NL_PERT_FIRE} == true ]]; then
            if [[ ! -e ${WRFFIRECHEMI} ]]; then
               if [[ ${USE_HSI} == true ]]; then 
                  hsi -d4 get ${HSI_DATA_DIR}/${YYYY}${MM}${DD}/${WRFFIRECHEMI}
               else
                  cp ${DATA_DIR}/${YYYY}${MM}${DD}/${WRFFIRECHEMI} ./.
               fi
               cp ${WRFFIRECHEMI} ${WRFFIRECHEMI}.${CMEM}
            else
               cp ${WRFFIRECHEMI} ${WRFFIRECHEMI}.${CMEM}
            fi
            mv ${WRFFIRECHEMI} ${SAVE_DIR}/${YYYY}${MM}${DD}/
         fi
         if [[ ${NL_PERT_BIO} == true ]]; then
            if [[ ! -e  ${WRFBIOCHEMI} ]]; then
               if [[ ${USE_HSI} == true ]]; then 
                  hsi -d4 get ${HSI_DATA_DIR}/${YYYY}${MM}${DD}/${WRFBIOCHEMI}
               else
                  cp ${DATA_DIR}/${YYYY}${MM}${DD}/${WRFBIOCHEMI} ./.
               fi
               cp ${WRFBIOCHEMI} ${WRFBIOCHEMI}.${CMEM}
            else
               cp ${WRFBIOCHEMI} ${WRFBIOCHEMI}.${CMEM}
            fi
            mv ${WRFBIOCHEMI} ${SAVE_DIR}/${YYYY}${MM}${DD}/
         fi
#
# CREATE NAMELIST
         rm -rf perturb_chem_emiss_nml.nl
         cat << EOF > perturb_chem_emiss_nml.nl
&perturb_chem_emiss_nml
idate=${DATE},
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
# RUN PERTRUBATION CODE    
         ./perturb_chem_emiss.exe
#
# SAVE PERTURBED EMISSION FILES
         if [[ ! -e ${SAVE_DIR}/${YYYY}${MM}${DD} ]]; then
            mkdir -p ${DATA_DIR}/${YYYY}${MM}${DD}
         fi
         if [[ ${NL_PERT_CHEM} == true ]]; then
            mv ${WRFCHEMI}.${CMEM} ${SAVE_DIR}/${YYYY}${MM}${DD}/
         fi
         if [[ ${NL_PERT_FIRE} == true ]]; then
            mv ${WRFFIRECHEMI}.${CMEM} ${SAVE_DIR}/${YYYY}${MM}${DD}/
         fi
         if [[ ${NL_PERT_BIO} == true ]]; then
            mv ${WRFBIOCHEMI}.${CMEM} ${SAVE_DIR}/${YYYY}${MM}${DD}/
         fi
#
# GO TO NEXT MEMBER
         let MEM=${MEM}+1
      done
# ADVANCE TIME      
      export DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} ${EMISS_FREQ} 2>/dev/null)
      export YYYY=$(echo $DATE | cut -c1-4)
      export MM=$(echo $DATE | cut -c5-6)
      export DD=$(echo $DATE | cut -c7-8)
      export HH=$(echo $DATE | cut -c9-10)
done
