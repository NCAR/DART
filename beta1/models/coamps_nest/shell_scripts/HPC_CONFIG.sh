#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################

function PBS_DIRECTIVE
{

usage="Usage: `basename $0` expname cdtg nprocs nnodes walltime logdir serial(y|n) job_array(y|n)" 

  EXPNAME_LOC=$1
  CDTG_LOC=$2
  NPROC_LOC=$3
  NNODES_LOC=$4
  WALLTIME_LOC=$5
  LOGDIR_LOC=$6
  serial=$7
  job_array=$8
  PBS_QUEUE=$9

: ${serial:='n'}
: ${job_array:='y'}

if [ $serial == 'y' ]; then 
  CPU_PER_NODE_LOC=${NNODES_LOC}
: ${PBS_QUEUE:=$PBS_QUEUE_TRANS}
else
  CPU_PER_NODE_LOC=${CPU_PER_NODE}
: ${PBS_QUEUE:=$PBS_QUEUE_FILTER}
fi
let "MPI_PROCS = (NPROC_LOC / NNODES_LOC) + (NPROC_LOC % CPU_PER_NODE_LOC)"

case ${PBS_QUEUE} in 
  debug)     ACCOUNT_LOC=${ACCOUNT_STANDARD}  ;;
  standard)  ACCOUNT_LOC=${ACCOUNT_STANDARD}  ;;
  challenge) ACCOUNT_LOC=${ACCOUNT_CHALLENGE} ;;
  transfer)  ACCOUNT_LOC=${ACCOUNT_STANDARD}  ;;
  *)         ACCOUNT_LOC=${ACCOUNT}           ;;
esac

  EXPNAME_CUT=`echo $EXPNAME_LOC | cut -c1-15`
# DIRECTIVES COMMON TO ALL SYSTEMS
  echo "#"
  echo "#### AUTOMATICALLY GENERATED PBS DIRECTIVES FOR ${HOST} ####"
  echo "#"
  echo "#PBS -S /bin/bash"
  echo "#PBS -V"
  echo "#PBS -N $EXPNAME_CUT"

  case $HOST in 
    (davinci)
    echo "#PBS -q ${PBS_QUEUE}"
    echo "#PBS -A ${ACCOUNT_LOC}"
    echo "#PBS -l walltime=${WALLTIME_LOC}"
    echo "#PBS -l select=${NNODES_LOC}:ncpus=${NPROC_LOC}:mpiprocs=${NPROC_LOC}"
    echo "#PBS -l place=scatter:excl"
#    if [ $job_array == 'y' ]; then echo "#PBS -J JOB_ARRAY"; fi
	;;
    (einstein)
    echo "#PBS -q ${PBS_QUEUE}"
    echo "#PBS -A ${ACCOUNT_LOC}"
    echo "#PBS -l walltime=$WALLTIME_LOC"
    echo "#PBS -l mppnppn=${CPU_PER_NODE_LOC}"
    if [ $serial == 'n' ]; then echo "#PBS -l mppwidth=${NPROC_LOC}"; fi
    if [ $job_array == 'y' ]; then echo "#PBS -J JOB_ARRAY"; fi
	;;
    (diamond)
    echo "#PBS -q ${PBS_QUEUE}"
    echo "#PBS -A ${ACCOUNT_LOC}"
    echo "#PBS -l walltime=$WALLTIME_LOC"
    if [ $serial == 'n' ]; then 
	  echo "#PBS -l select=${NNODES_LOC}:ncpus=${CPU_PER_NODE}:mpiprocs=${MPI_PROCS}"
      echo "#PBS -l place=scatter:excl"
    else
	  echo "#PBS -l select=${NNODES_LOC}:ncpus=${NNODES_LOC}"
    fi
    if [ $job_array == 'y' ]; then echo "#PBS -J JOB_ARRAY"; fi
	;;
    (hawk)
    echo "#PBS -q ${PBS_QUEUE}"
    echo "#PBS -A ${ACCOUNT_LOC}"
    echo "#PBS -l job_type=MPI"
    echo "#PBS -l walltime=$WALLTIME_LOC"
    echo "#PBS -l select=ncpus=${NPROC_LOC}"
    if [ $job_array == 'y' ]; then echo "#PBS -J JOB_ARRAY"; fi
	;;
    (maury)
    echo "#PBS -q ${PBS_QUEUE}"
    echo "#PBS -l nodes=${NNODES_LOC}:ppn=${CPU_PER_NODE_LOC}"
    if [ $job_array == 'y' ]; then echo "#PBS -t JOB_ARRAY"; fi
  esac

# MORE DIRECTIVES COMMON TO ALL SYSTEMS
  echo "#PBS -r y"
  echo "#PBS -m e"
  if [ $job_array == 'y' ]; then
    echo "#PBS -o ${LOGDIR_LOC}/${EXPNAME_LOC}.^array_index^.${CDTG_LOC}.out"
    echo "#PBS -e ${LOGDIR_LOC}/${EXPNAME_LOC}.^array_index^.${CDTG_LOC}.err"
  else
    echo "#PBS -o ${LOGDIR_LOC}/${EXPNAME_LOC}.${CDTG_LOC}.out"
    echo "#PBS -e ${LOGDIR_LOC}/${EXPNAME_LOC}.${CDTG_LOC}.err"
  fi
  echo "#"
}


function MACH_LAYOUT
{

machine=`uname -n`
kernal=`uname -s`
ACCOUNT_STANDARD='' # MUST HAVE AN ACCOUNT NUMBER
ACCOUNT_CHALLENGE=''# MUST HAVE AN ACCOUNT NUMBER
ACCOUNT=${ACCOUNT_STANDARD}

JOB_ARRAY='J'
JOB_ARRAY_AVAILABLE='y'

if [ $kernal = "AIX" ]; then
  JOB_ARRAY_AVAILABLE='n'
  PBS_QUEUE_ENS='challenge'
  PBS_QUEUE_FILTER='challenge'
  PBS_QUEUE_TRANS='transfer'
  ARCMACH=newton
  if [ `echo $machine | cut -c1` = "b" ]; then
    HOST=babbage
    CPU_PER_NODE=16
    FILTER_MPI_CMD='/usr/bin/poe'
    COAMM_MPI_CMD='/usr/bin/poe'
    COAMA_MPI_CMD='/usr/bin/poe'
    POST_MPI_CMD='/usr/bin/poe'
  elif [ `echo $machine | cut -c1` = "d" ]; then
    HOST=davinci
    LSMP='true'
    LSMP=`echo $LSMP | cut -c1`
    if [ $LSMP = 'T' -o $LSMP = 't' ]; then
      CPU_PER_NODE=64
      FILTER_MPI_CMD='/usr/bin/poe /site/bin/launch'
      COAMM_MPI_CMD='/usr/bin/poe /site/bin/launch'
      POST_MPI_CMD='/usr/bin/poe /site/bin/launch'
      COAMA_MPI_CMD='/usr/bin/poe -proc 1'
    else
      CPU_PER_NODE=32
      FILTER_MPI_CMD='/usr/bin/poe'
      COAMM_MPI_CMD='/usr/bin/poe'
      POST_MPI_CMD='/usr/bin/poe'
      COAMA_MPI_CMD='/usr/bin/poe -proc 1'
    fi
	MIN_PROCS=1
  fi
  MPI_CMD=`echo $MPI_CMD | sed 's/\//\\\\\//g'`
elif [ $kernal = "LINUX" -o $kernal = "Linux" ]; then
  cray=`echo $machine | cut -c1-3` 
  maury=`echo $machine | cut -c1-5` 
  hawk=`echo $machine | cut -c1-4` 
  diamond=`echo $machine | cut -c1-7 | sed 's/[^a-z]//g'` 
  if [ $diamond = "rin" ]; then diamond='diamond'; fi
  if [ $maury = "maury" ]; then
    JOB_ARRAY='t'
    HOST=maury
    PBS_QUEUE_ENS=debug
    PBS_QUEUE_FILTER=debug
    PBS_QUEUE_TRANS=debug
    ARCMACH=maury
    CPU_PER_NODE=4
	MIN_PROCS=1
	FILTER_MPI_CMD="mpirun -n ${FILTER_PROCS}"
	COAMM_MPI_CMD="mpirun -n ${COAMM_PROCS}"
	POST_MPI_CMD="mpirun -n ${POST_PROCS}"
	COAMA_MPI_CMD=""
  elif [ $cray = "nid" ]; then
    HOST=einstein
    PBS_QUEUE_ENS='standard'
	ACCOUNT_ENS=${ACCOUNT_STANDARD}
    PBS_QUEUE_FILTER='challenge'
    PBS_QUEUE_TRANS='transfer'
    ARCMACH=newton
    CPU_PER_NODE=8
	MIN_PROCS=1
	FILTER_MPI_CMD="aprun -n ${FILTER_PROCS} -N ${CPU_PER_NODE}"
	COAMM_MPI_CMD="aprun -n ${COAMM_PROCS} -N ${CPU_PER_NODE}"
	POST_MPI_CMD="aprun -n ${POST_PROCS} -N ${CPU_PER_NODE}"
	COAMA_MPI_CMD="aprun -n 1"
  elif [ $hawk = "hawk" ]; then
    HOST=hawk
    ARCMACH=${msas}
    PBS_QUEUE_ENS='challenge'
    PBS_QUEUE_FILTER='challenge'
    PBS_QUEUE_TRANS='challenge'
    CPU_PER_NODE=500
	MIN_PROCS=4
	FILTER_MPI_CMD="mpirun -n ${FILTER_PROCS}"
	COAMM_MPI_CMD="mpirun -n ${COAMM_PROCS}"
	POST_MPI_CMD="mpirun -n ${POST_PROCS}"
	COAMA_MPI_CMD="mpirun -n 1"
  elif [ $diamond = "diamond" ]; then
    JOB_ARRAY_AVAILABLE='n'
    HOST=diamond
    ARCMACH=${ARCHIVE_HOST}
    #PBS_QUEUE='R186875'
	PBS_QUEUE='standard'
    PBS_QUEUE_ENS=${PBS_QUEUE}
    PBS_QUEUE_FILTER=${PBS_QUEUE}
    PBS_QUEUE_TRANS='transfer'
    #PBS_QUEUE_TRANS=${PBS_QUEUE}
    CPU_PER_NODE=8
	MIN_PROCS=1

	FILTER_MPI_CMD="mpiexec_mpt -np ${FILTER_PROCS}"
	COAMM_MPI_CMD="mpiexec_mpt -np ${COAMM_PROCS}"
	POST_MPI_CMD="mpiexec_mpt -np ${POST_PROCS}"
	COAMA_MPI_CMD=""

	export MPI_GROUP_MAX=256
  fi
fi

export ARCMACH
export PBS_QUEUE_ENS
export PBS_QUEUE_FILTER
export HOST
export CPU_PER_NODE
export MIN_PROCS
export FILTER_MPI_CMD
export POST_MPI_CMD
export COAMM_MPI_CMD
export COAMA_MPI_CMD
export ACCOUNT
export ACCOUNT_STANDARD
export ACCOUNT_CHALLENGE
export JOB_ARRAY_AVAILABLE

}

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

