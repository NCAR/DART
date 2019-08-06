#!/bin/bash
#BSUB -J ENSDART2M[1-ENSEMBLEMEMBERNO]         # Name of the job.
#BSUB -o LOG/TSSD2M_%J_%I.out  # Appends std output to file %J.out.
#BSUB -e LOG/TSSD2M_%J_%I.out  # Appends std error to file %J.err.
#BSUB -q serial_30min            # queue
JOBDIR=${LS_SUBCWD}                   # directory of this script
JOBIDN=$( echo ${LSB_JOBID} | awk '{ printf("%08d\n", $1) }' ) # job-id
JOBNAM=${LSB_JOBNAME}                 # name of this script
#-- Load Experiment Environment Variables -----------------     
. FeoM_SBMT_ENV_VARS.sh
#-- Ensemble Required Variables ---------------------------     
ENSNO=$( echo ${LSB_JOBINDEX} | awk '{ printf("%02d\n", $1) }' )
ENSN4=$( echo ${LSB_JOBINDEX} | awk '{ printf("%04d\n", $1) }' )
ENSINFO=${ENSID}${ENSNO}
ENSDIR=${WRKDIR}/${ENSINFO}; cd ${ENSDIR}
#-- Set and Copy The Executables --------------------------
EXE=dart_to_model
DART2M=${DRTDIR}/${EXE}; ${COPY} ${DART2M} .
#-- DART I/O files ----------------------------------------
M2DOUT=filter_ics.${ENSN4}
D2MINP=filter_restart.${ENSN4}
${LINK} ${FILDIR}/${D2MINP} .
        ./${EXE} 
