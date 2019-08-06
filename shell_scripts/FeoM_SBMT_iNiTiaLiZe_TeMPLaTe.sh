#!/bin/bash
#BSUB -J TSSENSINI[1-ENSEMBLEMEMBERNO]       # Name of the job.
#BSUB -o LOG/TSSENSINI_%J_%I.out  # Appends stdout to file %J.out.
#BSUB -e LOG/TSSENSINI_%J_%I.out  # Appends stderr to file %J.err.
#BSUB -q serial_30min            # queue
############################ LSF ###################################
JOBDIR=${LS_SUBCWD}                   # directory of this script
JOBIDN=$( echo ${LSB_JOBID} | awk '{ printf("%08d\n", $1) }' ) # job-id
JOBNAM=${LSB_JOBNAME}                 # name of this script
#-- Load Experiment Environment Variables -----------------     
. FeoM_SBMT_ENV_VARS.sh
#-- Ensemble required variables ---------------------------
ENSNO=$( echo ${LSB_JOBINDEX} | awk '{ printf("%02d\n", $1) }' )
ENSINFO=${ENSID}${ENSNO};
JOBID=$(bjobs | grep -ir TSS${EXPINFO} | awk '{print $1}')
ENSDIR=${WRKDIR}/${ENSINFO}
SUBMITFILE=${WRKDIR}/FeoM_SBMT_iNiTiaLiZe_${EXPINFO}
#----------------------------------------------------------
#-------- CREATE ENSEMBLE MEMBER DIRECTORY ----------------
#----------------------------------------------------------
if [ ! -d ${ENSDIR} ]; then
	mkdir ${ENSDIR}
fi
cd ${ENSDIR}
#----------------------------------------------------------
#-------- LINK INITIAL CONDITIONS -------------------------
#----------------------------------------------------------
if [ ! -f  ${ENSINFO}.clock ]; then 
	if [ ${EXPYR} -eq 2009 ]; then
		INIYR=2008
	cp ${FSMINI}/ENSHC.${INIYR}.clock ${ENSINFO}.clock
	ln -sf ${FSMINI}/ENSHC.${INIYR}.forcing.diag.nc .
	ln -sf ${FSMINI}/ENSHC.${INIYR}.oce.diag.nc .
	ln -sf ${FSMINI}/ENSHC.${INIYR}.ice.diag.nc .
	ln -sf ${FSMINI}/ENSHC.${INIYR}.oce.mean.nc .
	ln -sf ${FSMINI}/ENSHC.${INIYR}.ice.mean.nc .
	ln -sf ${FSMINI}/ENSHC.${INIYR}.ice.nc .

	ln -sf ${FSMINI}/${ENSINFO}.${INIYR}.oce.nc .
	fi
#----------------------------------------------------------
#-------- SET FEOM NAMELIST PARAMETERS --------------------
#----------------------------------------------------------
	for File in ENSHC.${INIYR}.*.nc; do
		NFile=$( echo ${File} | sed 's/ENSHC/'${ENSINFO}'/' );
		mv ${File} ${NFile};
	done
	sed -e  "s/EXPNUM/${EXPNO}/" -e "s/EXPDEF/${EXPID}/" -e \
	        "s/ENSNUM/${ENSNO}/" -e "s/ENSDEF/${ENSID}/" -e \
	        "s/TIMESTEP/7200/" -e  "s/MESHDIR/${MESHDIR}/" -e \
	        "s/RUNLENGTH/${RNLEN}/" \
		${MODELHOM}/namelist.config.template > ${ENSDIR}/namelist.config
	sed -e  "s/WNDFORCING/MFS/" -e "s/RADFORCING/MFS/" -e \
		"s/PRCFORCING/NCEP/" -e "s/ROFFORCING/KARA/" -e \
		"s/SSSFORCING/MFS/" \
		${MODELHOM}/namelist.forcing.template > ${ENSDIR}/namelist.forcing
	cp ${MODELHOM}/namelist.diag namelist.diag
	cp ${MODELHOM}/namelist.ice namelist.ice
	cp ${MODELHOM}/namelist.oce namelist.oce
	cp ${DRTDIR}/FeoM_time ${ENSDIR}/.
	sed -e  "s/ENSNUM/${ENSNO}/" -e "s/ENSDEF/${ENSID}/"  \
	        ${MODELHOM}/namelist.ensemble > ${ENSDIR}/namelist.ensemble
        cp ${MODELHOM}/fesom.x fesom.x
fi
