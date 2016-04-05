#!/bin/ksh -aeux
#
# Set experiment parameters
export DATE_START=2008061006
export DATE_END=2008062000
export DATE_START=2008060112
export DATE_END=2008062900
#
#export DATE_START=2008060112
#export DATE_END=2008060200
#
export TIME_INC=6
export DELETE_FLG=false
#
# Script to get files from HSI:
export WRFDA_VERSION=WRFDAv3.4_dmpar
export BUILD_DIR=/glade/p/work/mizzi/TRUNK/${WRFDA_VERSION}/var/build
#
# Experiment directories
#export EXP_DIR=MOPCOMB_Exp_2_RtDA_40M_p30p30_sp4
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_40M_p30p30_sp4
#export EXP_DIR=MOPCOMB_Exp_3_MgDA_40M_p30p30_sp4
#
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_p10p30
#export EXP_DIR=MOPCOMB_Exp_3_MgDA_20M_p10p00
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_p10p30_loc
#
#export EXP_DIR=MOPCOMB_Exp_3_MgDA_20M_100km_p10p00
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_100km_p10p30
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_100km_loc_p10p30
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_100km_bar_p10p30
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_100km_NoRot_p10p30
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_100km_NoRot_loc_p10p30
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_100km_NoRot_bloc_p10p30
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_100km_TEST_p75
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_100km_TEST_p50
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_100km_DBL_bloc_p10p30
export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_100km_DBL_NV_p10p30
#
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_p10p00
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_p10p30
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_p20p00
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_p30p00
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_p30p30
#
export HSI_DIR=/MIZZI/DART_TEST_AVE/${EXP_DIR}
export SCRATCH_DIR=/glade/scratch/mizzi/DART_TEST_AVE/${EXP_DIR}
export ACD_DIR=/glade/p/acd/mizzi/DART_TEST_AVE/${EXP_DIR}
cd ${ACD_DIR}
export L_DATE=${DATE_START}
#
# Copy file into ${L_DATE} subdirectory
while [[ ${L_DATE} -le ${DATE_END} ]] ; do
   echo 'copy '${L_DATE}
   if [[ ! -d ${ACD_DIR}/${L_DATE}/dart_filter ]]; then
      mkdir -p ${ACD_DIR}/${L_DATE}/dart_filter
   fi
   cd ${ACD_DIR}/${L_DATE}/dart_filter
   if [[ -f Prior_Diag.nc && ${DELETE_FLG} == true ]]; then rm -rf Prior_Diag.nc; fi
   if [[ -f Posterior_Diag.nc && ${DELETE_FLG} == true ]]; then rm -rf Posterior_Diag.nc; fi
   if [[ ! -f Prior_Diag.nc ]]; then
      if [[ -e ${SCRATCH_DIR}/${L_DATE}/dart_filter/Prior_Diag.nc ]]; then
         cp ${SCRATCH_DIR}/${L_DATE}/dart_filter/Prior_Diag.nc ./.
      else  
         hsi get ${ACD_DIR}/${L_DATE}/dart_filter/Prior_Diag.nc : ${HSI_DIR}/${L_DATE}/dart_filter/Prior_Diag.nc
      fi
   fi
   if [[ ! -f Posterior_Diag.nc ]]; then 
      if [[ -e ${SCRATCH_DIR}/${L_DATE}/dart_filter/Posterior_Diag.nc ]]; then
         cp ${SCRATCH_DIR}/${L_DATE}/dart_filter/Posterior_Diag.nc ./.
      else  
         hsi get ${ACD_DIR}/${L_DATE}/dart_filter/Posterior_Diag.nc : ${HSI_DIR}/${L_DATE}/dart_filter/Posterior_Diag.nc
      fi
   fi
   export P_DATE=${L_DATE}
   export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${TIME_INC} 2>/dev/null)
done
#
# Concatenate files
cd ${ACD_DIR}
export PRIOR_CAT_FILE=${ACD_DIR}/Cat_Prior_Diag.nc
export POST_CAT_FILE=${ACD_DIR}/Cat_Posterior_Diag.nc
if [[ -f ${PRIOR_CAT_FILE} ]]; then
  rm ${PRIOR_CAT_FILE}
fi
if [[ -f ${POST_CAT_FILE} ]]; then
  rm ${POST_CAT_FILE}
fi
#
# Copy first date
export L_DATE=${DATE_START}
cp ${ACD_DIR}/${L_DATE}/dart_filter/Prior_Diag.nc old_prior_file.nc 
cp ${ACD_DIR}/${L_DATE}/dart_filter/Posterior_Diag.nc old_post_file.nc 
#
# Advance date
export P_DATE=${L_DATE}
export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${TIME_INC} 2>/dev/null)
#
# Concatenate remaining dates
while [[ ${L_DATE} -le ${DATE_END} ]] ; do
   echo 'copy '${L_DATE}
   export NEW_PRIOR_FILE=${ACD_DIR}/${L_DATE}/dart_filter/Prior_Diag.nc
   export NEW_POST_FILE=${ACD_DIR}/${L_DATE}/dart_filter/Posterior_Diag.nc
   ncrcat old_prior_file.nc ${NEW_PRIOR_FILE} new_prior_file.nc
   ncrcat old_post_file.nc ${NEW_POST_FILE} new_post_file.nc
   rm old_prior_file.nc
   rm old_post_file.nc
   mv new_prior_file.nc old_prior_file.nc 
   mv new_post_file.nc old_post_file.nc 
   export P_DATE=${L_DATE}
   export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${TIME_INC} 2>/dev/null)
done
#
# Save archive files
mv old_prior_file.nc ${PRIOR_CAT_FILE}
mv old_post_file.nc ${POST_CAT_FILE}
exit
