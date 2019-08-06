#!/bin/bash
#----------------------------------------------------------
#-- Environment variables common to all SBMT_FILES --------
#-- They can be common directories, common variables etc. -
#----------------------------------------------------------
#-- Some bash commands are assigned to variables ----------
export TEMPLATE=TeMPLaTe; export COPY='cp -f' 
export REMOVE='rm -f'; export LINK='ln -sf' 
export MOVE='mv -f'
export USRNAM=ans051
#-- LSF required variables --------------------------------
export NPROC=512           # How many cores to run FEOM
export POENAME=poe_short   # Name of the lsf queue to submit
export MPIEXEC=/users/home/opt/lsf/8.0/linux2.6-glibc2.3-x86_64/bin/mpirun.lsf
#-- Set experiment parameters here  -----------------------
export EXPID=FB1           # id of the experiment
export EXPNO=02            # nb of the experiment 
export EXPYR=2009          # experiment initial year
export EXPINFO=${EXPID}${EXPNO};
export MEMNO=30            # nb of ensemble members
export ENSID=ENS           # id of ensemble members
export ENDYR=2009          # year to end the experiment
export ENDDY=7             # day to end the experiment
# RNLEN is number of time step; use for FEOM
export RNLEN=1800           
# TSTEP is seconds in time step 
export TSTEP=12
# RNSEC is seconds to run; use for DART
export RNSEC=$( echo "${RNLEN} * ${TSTEP}" | bc )
export CYCLE=${RNSEC}
#-- Set the common directories here  ----------------------
export USRHOM=/users/home/${USRNAM}
export DRTDIR=${USRHOM}/DART/FEOM/models/FeoM/work
export WRKDIR=/work/${USRNAM}/TSS/${EXPINFO}
export RUNDIR=${USRHOM}/DART/FEOM/models/FeoM/shell_scripts
export DIADIR=${USRHOM}/FEOM_POSTPROC/MESH_READ
export LOGDIR=${WRKDIR}/LOG
export FILDIR=${WRKDIR}/FILTER
export CHECKFILE=${WRKDIR}/submitcheck.lst
export FSMHOM=${USRHOM}/FEOM
export MODELHOM=${FSMHOM}/FETSSOM.ENS01
export FSMPRE=${USRHOM}/FEOM_PREPROC
export FSMINI=${FSMPRE}/ENSEMBLE_IC 
export OBSSEQ=${DRTDIR}/obs_seq.ferrybox
#-- Mesh directories depending on the partitioning --------
if   [ ${NPROC} -eq 256 ]; then 
export MESHDIR=mesh-T2G1.5L110b
elif [ ${NPROC} -eq 1024 ]; then 
export MESHDIR=mesh-T2G1.5L110b.V2
elif [ ${NPROC} -eq 512 ]; then 
export MESHDIR=mesh-T2G1.5L110b.V3
elif [ ${NPROC} -eq 128 ]; then 
export MESHDIR=mesh-T2G1.5L110b.V4
fi
