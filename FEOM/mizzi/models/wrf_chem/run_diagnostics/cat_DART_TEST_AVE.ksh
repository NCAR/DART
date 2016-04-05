#!/bin/ksh -aeux
#
# Script to get files from HSI:

export WRFDA_VERSION=WRFDAv3.4_dmpar
export BUILD_DIR=/glade/p/work/mizzi/TRUNK/${WRFDA_VERSION}/var/build
#
#export EXP_DIR=MOPCOMB_Exp_2_RtDA_40M_p30p30_sp4
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_40M_p30p30_sp4
#export EXP_DIR=MOPCOMB_Exp_3_MgDA_40M_p30p30_sp4
#
#export EXP_DIR=MOPCOMB_Exp_3_MgDA_20M_p10p00
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_p10p30
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_p20p00
#export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_p30p00
export EXP_DIR=MOPCOMB_Exp_2_MgDA_20M_p30p30
#
export TARGET_PATH=/glade/p/acd/mizzi/DART_TEST_AVE/${EXP_DIR}
#
export DATE_START=2008061006
export DATE_END=2008062000
export TIME_INC=6
#
export L_DATE=${DATE_START}
cd ${TARGET_PATH}
export PRIOR_CAT_FILE=${TARGET_PATH}/Cat_Prior_Diag.nc
export POST_CAT_FILE=${TARGET_PATH}/Cat_Posterior_Diag.nc
if [[ -e ${PRIOR_CAT_FILE} ]]; then
  rm ${PRIOR_CAT_FILE}
fi
if [[ -e ${POST_CAT_FILE} ]]; then
  rm ${POST_CAT_FILE}
fi
cp ${TARGET_PATH}/${L_DATE}/dart_filter/Prior_Diag.nc old_prior_file.nc 
cp ${TARGET_PATH}/${L_DATE}/dart_filter/Posterior_Diag.nc old_post_file.nc 
export P_DATE=${L_DATE}
export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${TIME_INC} 2>/dev/null)
#
while [[ ${L_DATE} -le ${DATE_END} ]] ; do
   echo 'copy '${L_DATE}
   export NEW_PRIOR_FILE=${TARGET_PATH}/${L_DATE}/dart_filter/Prior_Diag.nc
   export NEW_POST_FILE=${TARGET_PATH}/${L_DATE}/dart_filter/Posterior_Diag.nc
   ncrcat old_prior_file.nc ${NEW_PRIOR_FILE} new_prior_file.nc
   ncrcat old_post_file.nc ${NEW_POST_FILE} new_post_file.nc
   rm old_prior_file.nc
   rm old_post_file.nc
   mv new_prior_file.nc old_prior_file.nc 
   mv new_post_file.nc old_post_file.nc 
   export P_DATE=${L_DATE}
   export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${TIME_INC} 2>/dev/null)
done
mv old_prior_file.nc ${PRIOR_CAT_FILE}
mv old_post_file.nc ${POST_CAT_FILE}
exit

