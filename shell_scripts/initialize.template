#!/bin/bash
#
#--- Initializes the ensemble if not yet done
#--- creates directories, copies initial ensemble, namelists ...
#--- Executed in a JOB ARRAY environment, only uses 1 task per job
#----------------------------------------------------------------------
# LSF options
#
#BSUB -J init_ens[1-ENSEMBLEMEMBERNO]   # Name of the job (array).
#BSUB -o LOG/init_ens_%J_%I.out         # Appends stdout to file %J.out.
#BSUB -e LOG/init_ens_%J_%I.out         # Appends stderr to file %J.err.
#BSUB -q serial_30min                   # queue
#
#----------------------------------------------------------------------
# PBS options  (set for the NCAR machine "cheyenne")
#
#PBS -N init_ens
#PBS -J 1-ENSEMBLEMEMBERNO
#PBS -l walltime=0:10:00
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
   JOBIDN=${LSB_JOBINDEX}      # job array index
   EXTENSION=lsf

elif [[ ${SCHEDULER} = "pbs" ]] ; then

   JOBDIR=${PBS_O_WORKDIR}     # directory of this script
   JOBNAM=${PBS_JOBNAME}       # name of this script
   JOBIDN=${PBS_ARRAY_INDEX}   # job array index
   TMPDIR=/glade/scratch/$USER/temp  # cheyenne-specific
   mkdir -p $TMPDIR                  # cheyenne-specific
   EXTENSION=pbs

fi

#-- Ensemble required variables ---------------------------

ENSNO=$( echo ${JOBIDN} | awk '{ printf("%02d\n", $1) }' )
ENSINFO=${ENSID}${ENSNO};
ENSDIR=${WRKDIR}/${ENSINFO}

#----------------------------------------------------------
#-------- CREATE ENSEMBLE MEMBER DIRECTORY ----------------
#----------------------------------------------------------

[ ! -d ${ENSDIR} ] && mkdir ${ENSDIR}
cd ${ENSDIR}

#----------------------------------------------------------
#-------- LINK INITIAL CONDITIONS -------------------------
#----------------------------------------------------------

if [ ! -f  ${ENSINFO}.clock ]; then
   INIYR=$[ ${EXPYR} - 1 ] # assuming experiment starts on Jan, 1
   cp ${FSMINI}/ENSHC.${INIYR}.clock ${ENSINFO}.clock
   ln -sf ${FSMINI}/ENSHC.${INIYR}.forcing.diag.nc .
   ln -sf ${FSMINI}/ENSHC.${INIYR}.oce.diag.nc .
   ln -sf ${FSMINI}/ENSHC.${INIYR}.ice.diag.nc .
   ln -sf ${FSMINI}/ENSHC.${INIYR}.oce.mean.nc .
   ln -sf ${FSMINI}/ENSHC.${INIYR}.ice.mean.nc .
   ln -sf ${FSMINI}/ENSHC.${INIYR}.ice.nc .
   ln -sf ${FSMINI}/${ENSINFO}.${INIYR}.oce.nc .

   #----------------------------------------------------------
   #-------- SET FESOM NAMELIST PARAMETERS -------------------
   #----------------------------------------------------------

   for File in ENSHC.${INIYR}.*.nc; do
     NFile=$( echo ${File} | sed 's/ENSHC/'${ENSINFO}'/' );
     mv ${File} ${NFile};
   done

   sed -e "s/EXPNUM/${EXPNO}/" -e "s/EXPDEF/${EXPID}/"    -e \
          "s/ENSNUM/${ENSNO}/" -e "s/ENSDEF/${ENSID}/"    -e \
          "s/TIMESTEP/7200/"   -e "s/MESHDIR/${MESHDIR}/" -e \
          "s/RUNLENGTH/${RNLEN}/" \
          ${MODELHOM}/namelist.config.template > ${ENSDIR}/namelist.config

   sed -e "s/WNDFORCING/MFS/"  -e "s/RADFORCING/MFS/"  -e \
          "s/PRCFORCING/NCEP/" -e "s/ROFFORCING/KARA/" -e \
          "s/SSSFORCING/MFS/" \
          ${MODELHOM}/namelist.forcing.template > ${ENSDIR}/namelist.forcing

   cp ${MODELHOM}/namelist.diag namelist.diag
   cp ${MODELHOM}/namelist.ice namelist.ice
   cp ${MODELHOM}/namelist.oce namelist.oce

   sed -e  "s/ENSNUM/${ENSNO}/" -e "s/ENSDEF/${ENSID}/"  \
           ${MODELHOM}/namelist.ensemble > ${ENSDIR}/namelist.ensemble

   cp ${MODELHOM}/fesom.x fesom.x
fi
