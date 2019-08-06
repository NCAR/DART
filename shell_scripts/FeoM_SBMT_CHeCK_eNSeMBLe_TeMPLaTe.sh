#!/bin/bash
#BSUB -J TSSENSCHeCK             # Name of the job.
#BSUB -o LOG/TSSENSCHeCK_%J.out        # Appends std output to file %J.out.
#BSUB -e LOG/TSSENSCHeCK_%J.out        # Appends std error  to file %J.out.
#BSUB -q serial_30min            # queue
# LSF
JOBDIR=${LS_SUBCWD}                   # directory of this script
JOBIDN=${LSB_JOBID}                   # job-id
JOBNAM=${LSB_JOBNAME}                 # name of this script
function jobid {
  output=$($*)
  echo $output | head -n1 | cut -d'<' -f2 | cut -d'>' -f1
    }
#-- Load Experiment Environment Variables -----------------     
. FeoM_SBMT_ENV_VARS.sh
#-- Ensemble required variables ---------------------------
ENSINFO=${ENSID}${MEMNO}; cd ${WRKDIR}
ENSCHECK=$( cat ${CHECKFILE} | wc -l )
  if [ ${ENSCHECK} -eq ${MEMNO} ]; then 
	SUBMITCHECK=$( awk '{SUM += $3} END {print SUM}' ${CHECKFILE} ) 
	if [ ${SUBMITCHECK} -ne 0 ]; then 
		echo "One of the ensemble members stopped. Ensemble collapsed..."
		exit
	else 
	cd ${WRKDIR}
	rm ${CHECKFILE}.prev; mv ${CHECKFILE} ${CHECKFILE}.prev
               #-------------------------------------------------------------
               # RUN THE MODEL 2 DART HERE
               #-------------------------------------------------------------
               TMPLFILE=${RUNDIR}/FeoM_SBMT_MoDeL_To_DaRT_TeMPLaTe.sh
               SBMTFILE=${WRKDIR}/FeoM_SBMT_MoDeL_To_DaRT_${EXPINFO}.sh
               sed "s;ENSEMBLEMEMBERNO;${MEMNO};g" ${TMPLFILE} > ${SBMTFILE}; cd ${WRKDIR}; 
               ID=$( jobid bsub < ${SBMTFILE} ); echo ${ID}
               #-------------------------------------------------------------
               # RUN THE FILTER HERE
               #-------------------------------------------------------------
               TMPLFILE=${RUNDIR}/FeoM_SBMT_FiLTeR_TeMPLaTe.sh
               SBMTFILE=${WRKDIR}/FeoM_SBMT_FiLTeR_${EXPINFO}.sh
               sed "s;ENSEMBLEMEMBERNO;${MEMNO};g" ${TMPLFILE} > ${SBMTFILE}; 
	       cd ${WRKDIR}; ID=$( jobid bsub -w "done(${ID})" < ${SBMTFILE} ); echo ${ID}
               #-------------------------------------------------------------
               # RUN THE DART 2 MODEL HERE
               #-------------------------------------------------------------
               TMPLFILE=${RUNDIR}/FeoM_SBMT_DaRT_To_MoDeL_TeMPLaTe.sh
               SBMTFILE=${WRKDIR}/FeoM_SBMT_DaRT_To_MoDeL_${EXPINFO}.sh
               sed "s;ENSEMBLEMEMBERNO;${MEMNO};g" ${TMPLFILE} > ${SBMTFILE}; cd ${WRKDIR}; 
               ID=$( jobid bsub -w "done(${ID})" < ${SBMTFILE} ); echo ${ID}
               #-------------------------------------------------------------
               # FINALIZE OR RESUBMIT THE JOB HERE
               #-------------------------------------------------------------
               TMPLFILE=${RUNDIR}/FeoM_SBMT_FiNaLiZe_TeMPLaTe.sh
               SBMTFILE=${WRKDIR}/FeoM_SBMT_FiNaLiZe_${EXPINFO}.sh
               sed -e "s;^EXPID=.*$;EXPID=${EXPID};g" -e \
                      "s;^EXPNO=.*$;EXPNO=${EXPNO};g" \
               ${TMPLFILE} > ${SBMTFILE}; cd ${WRKDIR}; 
	       ID=$( jobid bsub -w "done(${ID})" < ${SBMTFILE} ); echo ${ID}
	fi
  fi
