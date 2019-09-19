#!/bin/bash
#
#--- check if any members of the ensemble crashed
#--- check if the experiment finished
#--- if not, call filter
#----------------------------------------------------------------------
# LSF options
#
#BSUB -J ens_check               # Name of the job.
#BSUB -o LOG/ens_check_%J.out    # Appends std output to file %J.out.
#BSUB -e LOG/ens_check_%J.out    # Appends std error  to file %J.out.
#BSUB -q serial_30min            # queue
#
#----------------------------------------------------------------------
# PBS options  (set for the NCAR machine "cheyenne")
#
#PBS -N ens_check
#PBS -l walltime=0:03:00
#PBS -q regular
#PBS -j oe
#PBS -A P86850054
#PBS -l select=1:ncpus=1:mpiprocs=1
#
#-- Load Experiment Environment Variables -----------------

. environment.load

#----------------------------------------------------------------------
# Translate the queueing-specific variables into a common tongue.

if [[ $SCHEDULER = "lsf" ]] ; then

   SUB_CMD='bsub < '
   DEP_CMD='bsub -w "done(XXXXXXXX)" < '

elif [[ ${SCHEDULER} = "pbs" ]] ; then

   TMPDIR=/glade/scratch/$USER/temp  # cheyenne-specific
   mkdir -p $TMPDIR                  # cheyenne-specific
   SUB_CMD=qsub
   DEP_CMD='qsub -W depend=afterok:XXXXXXXX '

fi


#-- Ensemble required variables ---------------------------

ENSINFO=${ENSID}${MEMNO}; cd ${WRKDIR}
ENSCHECK=$( cat ${CHECKFILE} | wc -l )

if [ ${ENSCHECK} -eq ${MEMNO} ]; then
   SUBMITCHECK=$( awk '{SUM += $3} END {print SUM}' ${CHECKFILE} )
   if [ ${SUBMITCHECK} -ne 0 ]; then
      echo "ERROR: One of the ensemble members failed to advance."
      exit 1
   else
      cd ${WRKDIR}
      rm ${CHECKFILE}.prev
      mv ${CHECKFILE} ${CHECKFILE}.prev
 
      #-------------------------------------------------------------
      # Submit a job to run the filter, capture the Job ID for next step
      #-------------------------------------------------------------
 
      TMPLFILE=${RUNDIR}/filter.template
      SBMTFILE=${WRKDIR}/filter.${EXPINFO}
      sed "s;ENSEMBLEMEMBERNO;${MEMNO};g" ${TMPLFILE} > ${SBMTFILE};
 
      cd ${WRKDIR}
      ID=$( jobid eval ${SUB_CMD} ${SBMTFILE} )
      echo "Filter job ID is ${ID}"
 
      #-------------------------------------------------------------
      # Submit a job to the queue that will not get executed until
      # the filter.${EXPINFO} job has finished successfully.
      # FINALIZE (which may) RESUBMIT THE JOB
      #-------------------------------------------------------------
 
      TMPLFILE=${RUNDIR}/finalize.sh
      SBMTFILE=${WRKDIR}/finalize.${EXPINFO}
      ${COPY} ${TMPLFILE} ${SBMTFILE}
      cd ${WRKDIR}

      DEP_STRING=`echo ${DEP_CMD} | sed "s/XXXXXXXX/${ID}/"`

      ID=$( jobid eval ${DEP_STRING} ${SBMTFILE} )
      echo "Finalize job ID is ${ID}"
   fi
fi
