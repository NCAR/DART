#!/bin/bash
#BSUB -a poe
#BSUB -J ENSFILTER             # Name of the job.
#BSUB -o LOG/TSSFILTER_%J.out  # Appends std output to file %J.out.
#BSUB -e LOG/TSSFILTER_%J.out  # Appends std error to file %J.err.
#BSUB -q poe_short             # queue
#BSUB -R "span[ptile=16]"      #
#BSUB -n 16                    #
JOBDIR=${LS_SUBCWD}                   # directory of this script
JOBIDN=$( echo ${LSB_JOBID} | awk '{ printf("%08d\n", $1) }' ) # job-id
JOBNAM=${LSB_JOBNAME}                 # name of this script
#-- Load Experiment Environment Variables -----------------     
. FeoM_SBMT_ENV_VARS.sh
#-- Ensemble Required Variables ---------------------------     
ENSNO=$( echo 1 | awk '{ printf("%02d\n", $1) }' )
ENSN4=$( echo 1 | awk '{ printf("%04d\n", $1) }' )
ENSINFO=${ENSID}${ENSNO}
cd ${FILDIR}
#-- Set and Copy the Executables --------------------------
EXE=filter
FILTER=${DRTDIR}/${EXE}; ${COPY} ${FILTER} .
#-- DART I/O files ----------------------------------------
TMPNML=${DRTDIR}/input.nml.${TEMPLATE}
${COPY} ${OBSSEQ} obs_seq.out
M2DOUT=filter_ics.${ENSN4}
D2MINP=filter_restart.${ENSN4}
${COPY} ../ENS01/{input.nml,namelist.config} .
#-- DART FEOM time handshaking ----------------------------
DARTIME=($(grep 'FeoM    stop at :' ${WRKDIR}/ENS01/FeoM_time|sed 's;_\|:\|=\|,; ;g'| sed 's/[A-Za-z]*//g'))
DARTDAY=$( echo ${DARTIME[0]} | awk '{ printf("%06d\n", $1) }' )
DARTSEC=$( echo ${DARTIME[1]} | awk '{ printf("%05d\n", $1) }' )
${MPIEXEC} ./${EXE} 
	${COPY} obs_seq.final obs_seq.final.${DARTDAY}_${DARTSEC}
	${COPY} Posterior_Diag.nc Posterior_Diag_${DARTDAY}_${DARTSEC}.nc
	${COPY} Prior_Diag.nc Prior_Diag_${DARTDAY}_${DARTSEC}.nc
