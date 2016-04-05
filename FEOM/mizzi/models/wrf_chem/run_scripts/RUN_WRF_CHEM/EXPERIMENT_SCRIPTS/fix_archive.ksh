#!/bin/ksh -x
#
# Save WRFCHEM forecasts to archive directory
export RUN_DIR=/glade/scratch/mizzi/DART_TEST_AVE/MOPCOMB_Exp_2_MgDA_20M_p10p30_loc
export L_DATE=2008061206
export NUM_MEMBERS=20
export DOMAIN=01
export L_FILE_DATE=2008-06-12_06:00:00
export NEXT_FILE_DATE=2008-06-12_12:00:00
export CENTRALDIR=/glade/scratch/mizzi/DART_TEST_AVE/MOPCOMB_Exp_2_MgDA_20M_p10p30_loc/DART_CENTRALDIR
#
         if [[ ! -d ${RUN_DIR}/${L_DATE}/wrfchem_forecast ]]; then
            mkdir -p ${RUN_DIR}/${L_DATE}/wrfchem_forecast
            cd ${RUN_DIR}/${L_DATE}/wrfchem_forecast
         else
            cd ${RUN_DIR}/${L_DATE}/wrfchem_forecast
         fi
         let IMEM=1
         while [[ ${IMEM} -le ${NUM_MEMBERS} ]]; do
            export KMEM=${IMEM}
            if [[ ${IMEM} -lt 1000 ]]; then export KMEM=0${IMEM}; fi
            if [[ ${IMEM} -lt 100 ]]; then export KMEM=00${IMEM}; fi
            if [[ ${IMEM} -lt 10 ]]; then export KMEM=000${IMEM}; fi
            export WRFOUT_FILE_ANL=wrfout_d${DOMAIN}_${L_FILE_DATE} 
            export WRFOUT_FILE_FOR=wrfout_d${DOMAIN}_${NEXT_FILE_DATE} 
            export WRFOUT_FILE_APM=wrfapm_d${DOMAIN}_${NEXT_FILE_DATE} 
            cp ${CENTRALDIR}/advance_temp_${KMEM}/wrfinput_d${DOMAIN} wrfinput_d${DOMAIN}_${KMEM}
            cp ${CENTRALDIR}/advance_temp_${KMEM}/${WRFOUT_FILE_ANL} ${WRFOUT_FILE_ANL}_${KMEM}
            cp ${CENTRALDIR}/advance_temp_${KMEM}/${WRFOUT_FILE_FOR} ${WRFOUT_FILE_FOR}_${KMEM}
            cp ${CENTRALDIR}/advance_temp_${KMEM}/${WRFOUT_FILE_APM} ${WRFOUT_FILE_APM}_${KMEM}
            let IMEM=${IMEM}+1
         done
         cp ${CENTRALDIR}/advance_temp_0001/namelist.input ./.
exit
