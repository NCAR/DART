#!/bin/bash
#
#--- calls initialize the ensemble-------------------------
#--- calls forwarding the ensemble-------------------------
#--- calls check ensemble to decide if finished------------
#--- if not, submitted again until the experiment ends-----
#----------------------------------------------------------------------
# LSF options
#
#BSUB -J fesom_ens               # Name of the job.
#BSUB -o LOG/fesom_ens_%J.out    # Appends std output to file %J.out.
#BSUB -e LOG/fesom_ens_%J.out    # Appends std error  to file %J.out.
#BSUB -q serial_30min            # queue
#
#----------------------------------------------------------------------
# PBS options  (set for the NCAR machine "cheyenne")
#
#PBS -N fesom_ens
#PBS -l walltime=0:10:00
#PBS -q regular
#PBS -j oe
#PBS -A P86850054
#PBS -l select=1:ncpus=1:mpiprocs=1
#
#----------------------------------------------------------------------

#-- Load Experiment Environment Variables and Functions --

. environment.load


# Translate the queueing-specific variables into a common tongue.

if [[ $SCHEDULER = "lsf" ]] ; then

   JOBDIR=${LS_SUBCWD}         # directory of this script
   JOBNAM=${LSB_JOBNAME}       # name of this script
   SUB_CMD="bsub < "
   DEP_CMD='bsub -w "done(XXXXXXXX)" < '

   EXTENSION=lsf

elif [[ $SCHEDULER = "pbs" ]] ; then

   JOBDIR=${PBS_O_WORKDIR}     # directory of this script
   JOBNAM=${PBS_JOBNAME}       # name of this script
   TMPDIR=/glade/scratch/$USER/temp  # cheyenne-specific
   mkdir -p $TMPDIR                  # cheyenne-specific
   SUB_CMD=qsub
   DEP_CMD='qsub -W depend=afterok:XXXXXXXX '
   EXTENSION=pbs

fi

#-----------------------------------------------------------------------
# INITIALIZE ENSEMBLE IF NOT DONE YET
#-----------------------------------------------------------------------

if [ ! -d ${WRKDIR} ]; then

  # If the directory does not exist, there is some one-time setup
  # initialize.${EXPINFO} must be called. The initialize step is
  # embarassingly parallel, so a job array is used.

  DAYSTEP=1; RUNYEAR=2009

  mkdir ${WRKDIR}; mkdir ${LOGDIR}; mkdir ${FILDIR}

  ${COPY} ${RUNDIR}/environment.load ${WRKDIR}/.
  ${COPY} ${DRTDIR}/FESOM_time       ${FILDIR}

  TMPLFILE=${RUNDIR}/initialize.template
  SBMTFILE=${WRKDIR}/initialize.${EXPINFO}
  sed -e "s;ENSEMBLEMEMBERNO;${MEMNO};g" ${TMPLFILE} > ${SBMTFILE} || exit 1

  cd ${WRKDIR}
  ID=$( jobid eval ${SUB_CMD} ${SBMTFILE} )
  echo "Initialize job is ID ${ID}"

else
  cd ${WRKDIR}
  DAYSTEP=$(cat ${MEM01}/${MEM01}.clock | sed -n 2,2p | awk '{print $2}')
  RUNYEAR=$(cat ${MEM01}/${MEM01}.clock | sed -n 2,2p | awk '{print $3}')
  ${MOVE} ${CHECKFILE} ${CHECKFILE}.prev
  ID=0;
fi

#-----------------------------------------------------------------------
# ADVANCE THE MODEL
#-----------------------------------------------------------------------

TMPLFILE=${RUNDIR}/advance_model.template
SBMTFILE=${WRKDIR}/advance_model.${EXPINFO}
sed -e "s;ENSEMBLEMEMBERNO;${MEMNO};g" -e \
       "s;NUMBEROFCORES;${NPROC};"     -e \
       "s;POEQUEUENAME;${POENAME};" ${TMPLFILE} > ${SBMTFILE} || exit 1
if [ ${ID} = 0 ]; then
  cd ${WRKDIR}
  ID=$( jobid eval ${SUB_CMD} ${SBMTFILE} )
  echo "Model advance job is ID ${ID}"
else
  cd ${WRKDIR}
  DEP_STRING=`echo ${DEP_CMD} | sed "s/XXXXXXXX/${ID}/"`
  ID=$( jobid eval ${DEP_STRING} ${SBMTFILE} )
  echo "Model advance job is ID ${ID} ... queued and waiting."
fi

#-----------------------------------------------------------------------
# CHECK THE ENSEMBLE and RUN DART if all members advanced successfully.
# Resubmit the experiment if continues
#-----------------------------------------------------------------------

TMPLFILE=${RUNDIR}/check_ensemble.sh
SBMTFILE=${WRKDIR}/check_ensemble.${EXPINFO}
${COPY} ${TMPLFILE} ${SBMTFILE}

cd ${WRKDIR}
DEP_STRING=`echo ${DEP_CMD} | sed "s/XXXXXXXX/${ID}/"`
ID=$( jobid eval ${DEP_STRING} ${SBMTFILE} )
echo "Model advance check job is ${ID} ... queued and waiting.""

