#!/bin/bash
#BSUB -a poe                     
#BSUB -J TSSFWDFEOM[1-ENSEMBLEMEMBERNO]%3     # Name of the job.
#BSUB -o LOG/TSSFWDFEOM_%I_%J.out  # Appends stdout to file %J.out.
#BSUB -e LOG/TSSFWDFEOM_%I_%J.out  # Appends stderr to file %J.err.
#BSUB -P fesom                     # Project ID.
#BSUB -q POEQUEUENAME              # queue
#BSUB -R "span[ptile=16]"          #
#BSUB -n NUMBEROFCORES             #
#BSUB -N                           #
####BSUB -x                        #
#--  BEGIN ATHENA CONFIG ----------------------------------
 MPIPROGINF=detail
 export MPIPROGINF
 export LSF_PJL_TYPE="poe"
 export MEMORY_AFFINITY=MCM
 export MP_WAIT_MODE=poll
 export MP_SINGLE_THREAD=yes
 export MP_TASK_AFFINITY=MCM
 export MP_PGMMODEL=mpmd
 export MP_WAIT_MODE=poll
 export MP_POLLING_INTERVAL=30000000
 export MP_SHARED_MEMORY=yes
 export MP_EUILIB=us
 export MP_EUIDEVICE=sn_all
 export LDR_CNTRL=TEXTPSIZE=64K@STACKPSIZE=64K@DATAPSIZE=64K
 export MP_TASK_AFFINITY=core
#--  END ATHENA CONFIG ------------------------------------
F_RSVTASK=1; export F_RSVTASK; #THIS CAN BE USEFUL FOR ENSEMBLE. CHECK!
F_ERRCNT=0; export F_ERRCNT
JOBDIR=${LS_SUBCWD}                   # directory of this script
JOBIDN=$( echo ${LSB_JOBID} | awk '{ printf("%08d\n", $1) }' ) # job-id
JOBNAM=${LSB_JOBNAME}                 # name of this script
#-- Load Experiment Environment Variables -----------------     
. FeoM_SBMT_ENV_VARS.sh
#-- Ensemble required variables ---------------------------
ENSNO=$( echo ${LSB_JOBINDEX} | awk '{ printf("%02d\n", $1) }' )
ENSN4=$( echo ${LSB_JOBINDEX} | awk '{ printf("%04d\n", $1) }' )
ENSINFO=${ENSID}${ENSNO}; 
ENSDIR=${WRKDIR}/${ENSINFO}; 
#-- RUN THE FORWARD MODEL ---------------------------------
cd ${ENSDIR}
${MPIEXEC} ./fesom.x
#-- Check if the forward model blowed up ------------------
CHECKRETURN=$(bpeek -J TSSFWDFEOM[${LSB_JOBINDEX}] | grep -ir "The model blows up")
CHECKRESULT=$(echo $?); echo ${CHECKRESULT}
if [ ${CHECKRESULT} -ne "0" ];  then
	echo "${ENSINFO} EXIT 0 :is ready to be resubmitted" >> ${CHECKFILE}
else
	echo "${ENSINFO} EXIT 1 :is blowed up   " >> ${CHECKFILE}
fi
