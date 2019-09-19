#!/bin/bash
#
#--- Check if the experiment finished or not
#--- If finished, finalize.
#--- If not finished , resubmit.
#----------------------------------------------------------------------
# LSF options
#
#BSUB -J ens_finalize            # Name of the job.
#BSUB -o LOG/ens_finalize_%J.out # Appends std output to file %J.out.
#BSUB -e LOG/ens_finalize_%J.out # Appends std error  to file %J.out.
#BSUB -q serial_30min            # queue
#
#----------------------------------------------------------------------
# PBS options  (set for the NCAR machine "cheyenne")
#
#PBS -N ens_finalize
#PBS -l walltime=0:02:00
#PBS -q regular
#PBS -j oe
#PBS -A P86850054
#PBS -l select=1:ncpus=1:mpiprocs=1
#
#----------------------------------------------------------------------

#-- Load Experiment Environment Variables -----------------

. environment.load

# Translate the queueing-specific variables into a common tongue.

if [[ $SCHEDULER = "lsf" ]] ; then

   JOBDIR=${LS_SUBCWD}         # directory of this script
   JOBNAM=${LSB_JOBNAME}       # name of this script
   EXTENSION=lsf
   SUB_CMD='bsub < '

elif [[ ${SCHEDULER} = "pbs" ]] ; then

   JOBDIR=${PBS_O_WORKDIR}     # directory of this script
   JOBNAM=${PBS_JOBNAME}       # name of this script
   TMPDIR=/glade/scratch/$USER/temp  # cheyenne-specific
   mkdir -p $TMPDIR                  # cheyenne-specific
   EXTENSION=pbs
   SUB_CMD=qsub

fi

cd ${WRKDIR}

#----------------------------------------------------------
SBMTFILE=${RUNDIR}/ens_members.${EXPINFO}.sh

SECST=$(cat ${MEM01}/${MEM01}.clock | sed -n 2,2p | awk '{printf("%d\n"), $1}')
DAYST=$(cat ${MEM01}/${MEM01}.clock | sed -n 2,2p | awk '{print $2}')
RUNYR=$(cat ${MEM01}/${MEM01}.clock | sed -n 2,2p | awk '{print $3}')

if [ "${RUNYR}" -le "${ENDYR}" ] && [ "${DAYST}" -le "${ENDDY}" ]  && [ "${SECST}" -le 86400 ]; then
  eval ${SUB_CMD} ${SBMTFILE}
  echo "Experiment            at $RUNYR $DAYST $SECST"
  echo "Experiment continuing to $ENDYR $ENDDY 86400 ..."
else
  echo "Experiment Finished. Done."
  exit 0
fi
