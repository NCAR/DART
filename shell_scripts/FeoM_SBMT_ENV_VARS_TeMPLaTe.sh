#!/bin/bash
#----------------------------------------------------------
#-- Environment variables common to all SBMT_FILES --------
#-- They can be common directories, common variables etc. -
#----------------------------------------------------------
#-- Some bash commands are assigned to variables ----------
export TEMPLATE=TeMPLaTe; export COPY='cp -f' 
export REMOVE='rm -f'; export LINK='ln -sf' 
export MOVE='mv -f'
export USRNAM=             #username
#-- LSF required variables --------------------------------
export NPROC=512           # How many cores to run FEOM
export POENAME=poe_short   # Name of the lsf queue to submit
export MPIEXEC=/users/home/opt/lsf/8.0/linux2.6-glibc2.3-x86_64/bin/mpirun.lsf
#-- Set experiment parameters here  -----------------------
export EXPID= # id of the experiment
export EXPNO= # nb of the experiment 
export EXPYR= # experiment initial year
export EXPINFO=${EXPID}${EXPNO};
export MEMNO= # nb of ensemble members
export ENSID= # id of ensemble members
export ENDYR= # year to end the experiment
export ENDDY= # day to end the experiment
# RNLEN is number of time step; use for FEOM
export RNLEN= #             
# TSTEP is seconds in time step 
export TSTEP= #
# RNSEC is seconds to run; use for DART
export RNSEC=$( echo "${RNLEN} * ${TSTEP}" | bc )
export CYCLE=${RNSEC}
#-- Set the common directories here  ----------------------
export USRHOM= #user home directory
export DRTDIR= #dart work directory
export WRKDIR= #experiment directory 
export RUNDIR= #dart shell scripts directory
export DIADIR= #diagnostics directory
export LOGDIR= #lsf log file outputs
export FILDIR= #filter directory
export CHECKFILE= #checkfile for cycling FEOM+DART
export FSMHOM= #FEOM home directory
export MODELHOM= #FEOM code directory
export FSMPRE= #FEOM preproc directory
export FSMINI= #Initial conditions directory
export OBSSEQ= #obs_seq file
#-- Mesh directories depending on the partitioning --------
export MESHDIR= #FEOM mesh directory
