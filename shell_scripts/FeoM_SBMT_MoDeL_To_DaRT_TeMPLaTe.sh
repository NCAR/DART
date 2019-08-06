#!/bin/bash
#BSUB -J ENSM2DART[1-ENSEMBLEMEMBERNO]         # Name of the job.
#BSUB -o LOG/TSSM2D_%J_%I.out  # Appends std output to file %J.out.
#BSUB -e LOG/TSSM2D_%J_%I.out  # Appends std error to file %J.err.
#BSUB -q serial_30min            # queue
#-- LSF provided variable ---------------------------------
. FeoM_SBMT_ENV_VARS.sh
JOBDIR=${LS_SUBCWD}                   # directory of this script
JOBIDN=$( echo ${LSB_JOBID} | awk '{ printf("%08d\n", $1) }' ) # job-id
JOBNAM=${LSB_JOBNAME}                 # name of this script
ENSNO=$( echo ${LSB_JOBINDEX} | awk '{ printf("%02d\n", $1) }' )
ENSN4=$( echo ${LSB_JOBINDEX} | awk '{ printf("%04d\n", $1) }' )
ENSINFO=${ENSID}${ENSNO}
ENSDIR=${WRKDIR}/${ENSINFO}; cd ${ENSDIR}
#-- Set and Copy The Executables --------------------------
EXE=model_to_dart
M2DART=${DRTDIR}/${EXE}; ${COPY} ${M2DART} .
#-- DART I/O files ----------------------------------------
TMPNML=${DRTDIR}/input.nml.${TEMPLATE}
M2DOUT=filter_ics.${ENSN4}
D2MINP=filter_restart.${ENSN4}
#-- DART FEOM time handshaking ----------------------------
DART_INIT_TIME_DAYS=($(cat ${WRKDIR}/ENS01/FeoM_time|sed -n 1,1p|sed 's;_\|:\|=\|,; ;g'| sed 's/[A-Za-z]*//g'))
DART_INIT_TIME_SECS=($(cat ${WRKDIR}/ENS01/FeoM_time|sed -n 2,2p|sed 's;_\|:\|=\|,; ;g'| sed 's/[A-Za-z]*//g'))
DART_FIRST_OBS_DAYS=($(cat ${WRKDIR}/ENS01/FeoM_time|sed -n 3,3p|sed 's;_\|:\|=\|,; ;g'| sed 's/[A-Za-z]*//g'))
DART_FIRST_OBS_SECS=($(cat ${WRKDIR}/ENS01/FeoM_time|sed -n 4,4p|sed 's;_\|:\|=\|,; ;g'| sed 's/[A-Za-z]*//g'))
DART_LAST_OBS_DAYS=($(cat ${WRKDIR}/ENS01/FeoM_time|sed -n 5,5p|sed 's;_\|:\|=\|,; ;g'| sed 's/[A-Za-z]*//g'))
DART_LAST_OBS_SECS=($(cat ${WRKDIR}/ENS01/FeoM_time|sed -n 6,6p|sed 's;_\|:\|=\|,; ;g'| sed 's/[A-Za-z]*//g'))
ANALYSISFILE=${ENSDIR}/${ENSINFO}.${EXPYR}.oce.nc
echo "ANALYSIS FILE NAME: ${ANALYSISFILE}"
#-- Modify DART namelist template -------------------------
sed -e "s;^   model_analysis_filename.*$;   model_analysis_filename = '${ANALYSISFILE}';g" -e \
       "s/^   assimilation_period_days.*$/   assimilation_period_days = 0/g"  -e \
       "s/^   assimilation_period_seconds.*$/   assimilation_period_seconds = ${CYCLE}/g"  -e \
       "s/^   init_time_days.*$/   init_time_days = ${DART_INIT_TIME_DAYS}/g"              -e \
       "s/^   init_time_seconds.*$/   init_time_seconds = ${DART_INIT_TIME_SECS}/g"        -e \
       "s/^   first_obs_days.*$/   first_obs_days = ${DART_FIRST_OBS_DAYS}/g"              -e \
       "s/^   first_obs_seconds.*$/   first_obs_seconds = ${DART_FIRST_OBS_SECS}/g"        -e \
       "s/^   last_obs_days.*$/   last_obs_days = ${DART_LAST_OBS_DAYS}/g"                 -e \
       "s/^   last_obs_seconds.*$/   last_obs_seconds = ${DART_LAST_OBS_SECS}/g"           -e \
       "s/^   model_to_dart_output_file.*$/   model_to_dart_output_file = ${M2DOUT}/g"     -e \
       "s/^   dart_to_model_input_file.*$/   dart_to_model_input_file = ${D2MINP}/g"       -e \
       "s/^   num_output_state_members.*$/   num_output_state_members = ${MEMNO}/g"        -e \
       "s/^   num_output_obs_members.*$/   num_output_obs_members     = ${MEMNO}/g"        -e \
       "s/^   ens_size.*$/   ens_size      = ${MEMNO}/g"  \
        ${TMPNML} > input.nml
        ./${EXE} && ${LINK} ${ENSDIR}/${M2DOUT} ${FILDIR}/.
