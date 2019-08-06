#!/bin/bash
#BSUB -J TSSENSFINAL         # Name of the job.
#BSUB -o LOG/TSSENSFINAL_%J.out        # Appends std output to file %J.out.
#BSUB -e LOG/TSSENSFINAL_%J.out        # Appends std error  to file %J.out.
#BSUB -q serial_30min            # queue
# LSF
JOBDIR=${LS_SUBCWD}                   # directory of this script
JOBIDN=${LSB_JOBID}                   # job-id
JOBNAM=${LSB_JOBNAME}                 # name of this script
#-- Load Experiment Environment Variables -----------------
. FeoM_SBMT_ENV_VARS.sh
#----------------------------------------------------------
cd ${WRKDIR}
#----------------------------------------------------------
SBMTFILE=${RUNDIR}/FeoM_SBMT_eNS_MeMBeRS_${EXPINFO}.lsf 
SECST=$(cat ENS01/ENS01.clock | sed -n 2,2p | awk '{printf("%d\n"), $1}')
DAYST=$(cat ENS01/ENS01.clock | sed -n 2,2p | awk '{print $2}')
RUNYR=$(cat ENS01/ENS01.clock | sed -n 2,2p | awk '{print $3}') 
if [ "${RUNYR}" -le "${ENDYR}" ] && [ "${DAYST}" -le "${ENDDY}" ]  && [ "${SECST}" -le 86400 ]; then
		bsub < ${SBMTFILE}; echo "EXPERIMENT CONTINUES..."
else
	echo "EXPERIMENT FINISHED..."; exit
fi
