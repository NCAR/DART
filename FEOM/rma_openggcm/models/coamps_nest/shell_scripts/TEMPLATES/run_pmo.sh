#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
###############################################################################
#
# AUTHOR:	P. A. Reinecke
#           Naval Research Laboratory
#
###############################################################################
#
usage="Usage: `basename $0` -d dtg -i icycle" 
while getopts ":d:i:W:" option
do
  case $option in
      W  ) depend=$OPTARG;;
      d  ) dtg=$OPTARG;;
      i  ) icycle_loc=$OPTARG;;
      \? ) echo $usage ; exit 65;;
      *  ) echo $usage ; exit 65;;
  esac
done
BATCH="qsub"
SED='sed -i -e'

if [[ $depend ]]; then BATCH_DEPEND="${BATCH} -W depend=afterok:${depend}"; fi
: ${BATCH_DEPEND:=$BATCH}

PATH_CONFIG=CONFIG/FILE/GOES/HERE
. ${PATH_CONFIG}
: ${icycle_loc:=$icycle}

scrDir=${COAMPS_RUN_BASE_DIR}
LOG_DIR=${scrDir}/log
scriptDir=${scrDir}/scripts
DART_BIN=${scrDir}/bin

let icycle_sec='icycle * 3600'
tau_offset=`expr ${fcstlen_init} - ${PMO_FCST_HR}`
: ${dtg:=`echo ${cdtg_beg} ${tau_offset} | ${DART_BIN}/advance_time`}
: ${dtg_fob:=`echo ${dtg} ${PMO_FCST_HR} | ${DART_BIN}/advance_time`}

ntimes=0
dtg_cnt=${dtg_fob}
while [ $dtg_cnt -le $cdtg_end ]; do 
  ntimes=`expr $ntimes + 1` 
  dtg_cnt=`echo ${dtg_cnt} ${icycle} | ${DART_BIN}/advance_time`
done

stageIN=${scriptDir}/stageIN_pmo.${dtg}.sh
run_pmo=${scriptDir}/run_pmo.${dtg}.sh
queue_log=${LOG_DIR}/queue_${EXPNAME}.log

. ${scriptDir}/HPC_CONFIG.sh
MACH_LAYOUT
let "PMO_NODES = (PMO_PROCS / CPU_PER_NODE) + (PMO_PROCS % CPU_PER_NODE)"

if [ -e ${run_pmo} ]; then rm ${run_pmo}; fi
cat > ${run_pmo} << EOF
#!/bin/bash
`PBS_DIRECTIVE RUN_PMO ${dtg} ${PMO_PROCS} ${PMO_NODES} ${WALLTIME_PMO} ${LOG_DIR} n n` 
EOF

cat >> ${run_pmo} << EOF
scrDir=${scrDir}
scriptDir=\${scrDir}/scripts
DART_BIN=\${scrDir}/bin

dtg_fob=${dtg_fob}
dtg=${dtg}
icycle_sec=${icycle_sec}
ntimes=${ntimes}

# Decide if we're running an interactive job or not
# Modify this to allow for the possibility that it just doesn't exist!
if [[ \${PBS_ENVIRONMENT} && \${PBS_ENVIRONMENT-_} ]]; then
  case \$PBS_ENVIRONMENT in
    PBS_BATCH)
      echo "Batch mode detected."
      WORKDIR=\$PBS_O_WORKDIR ;;
    PBS_INTERACTIVE)
      echo "Interactive mode detected."
      WORKDIR=\`pwd\` ;;
    *)
      echo "\$PBS_ENVIRONMENT found, but not set!"
      exit 1 ;;
    esac
else
    WORKDIR=\`pwd\`
    echo "PBS_ENVIRONMENT is not defined, setting workdir to here"
fi
cd \${WORKDIR}

\${scriptDir}/update_perfect.sh -n 't' -d \${dtg} 

\${DART_BIN}/create_obs_sequence < ./coseq.dat

ds_time=( \`echo \${dtg_fob} 0 -g | \${DART_BIN}/advance_time\` )
CFN_FILE='./cfn.dat'
cat > \${CFN_FILE} << EOF_CFN
set_def.out
1
\${ntimes}
\${ds_time[0]} \${ds_time[1]}
0 \${icycle_sec}
obs_seq.in
EOF_CFN

\${DART_BIN}/create_fixed_network_seq < \${CFN_FILE}
rm \${CFN_FILE}

\${scriptDir}/advance_perfect.sh -d \${dtg}

\${DART_BIN}/perfect_model_obs

# Clean up a bit

rm -f \${WORKDIR}/filter_control*
rm -f \${WORKDIR}/assim_model_state_ic.0001*
rm -f \${WORKDIR}/assim_model_state_ud.0001*
rm -f \${WORKDIR}/perfect_ics*

rm -f \${WORKDIR}/datahd*\${dtg_beg}*
rm -f \${WORKDIR}/terrht*\${dtg_beg}* 

exit 0
EOF

# Script to stage files onto workspace
if [ -e ${stageIN} ]; then rm -f ${stageIN}; fi
cat > ${stageIN} << EOF_STAGE
#!/bin/bash
`PBS_DIRECTIVE STAGE_IN_PMO ${dtg} ${TRANS_PROCS} 1 ${WALLTIME_TRANS} ${LOG_DIR} y n` 
EOF_STAGE

cat >> ${stageIN} << EOF_STAGE

scrDir=${scrDir}
scriptDir=\${scrDir}/scripts
DART_BIN=\${scrDir}/bin
OSSE_FREE_RUN=${OSSE_FREE_RUN}

PMO_ARCH=${PMO_ARCH}
ARCMACH=${ARCMACH}

#dtg_beg=${cdtg_beg}
dtg_beg=${dtg}
dtg_end=${cdtg_end}
icycle=${icycle}
fcst_lead=${PMO_FCST_HR}
dtg=\${dtg_beg}

MEMBER_DIR=${scrDir}/perfect/data
if [ ! -e \${MEMBER_DIR} ]; then
  echo "  Creating directory \${MEMBER_DIR}..."
  mkdir -p \${MEMBER_DIR}
fi

if [ \${OSSE_FREE_RUN} == 'T' ]; then
 cd \${MEMBER_DIR}
 ARCHFILES=(\${dtg}_sgl.tar \${dtg}_sfc.tar)
 for ARCHFILE in \${ARCHFILES[@]}; do
   rsh \${ARCMACH} "stage -x \${ARCHFILE}"
   rsh \${ARCMACH} "stage \${ARCHFILE}"
   rcp \${ARCMACH}:\${PMO_ARCH}/\${ARCHFILE} .
   tar xvf \${ARCHFILE}
   rm -f \${ARCHFILE}
 done
else
 while [ \${dtg} -le \${dtg_end} ]; do
   cd \${scrDir}
   dtg_init=\`echo \${dtg} -\${fcst_lead} | \${DART_BIN}/advance_time\`

   cd \${MEMBER_DIR}
   ARCHFILE=\${dtg_init}_sgl\${fcst_lead}.tar
   rsh \${ARCMACH} "stage -x \${ARCHFILE}"
   rsh \${ARCMACH} "stage \${ARCHFILE}"
   rcp \${ARCMACH}:\${PMO_ARCH}/\${ARCHFILE} .
   tar xvf \${ARCHFILE}
   rm -f \${ARCHFILE}

   cd \${scrDir}
   dtg=\`echo \${dtg} \${icycle} | \${DART_BIN}/advance_time\`
 done
fi

cp \${MEMBER_DIR}/datahd*\${dtg_beg}* \${scrDir}
cp \${MEMBER_DIR}/terrht*\${dtg_beg}* \${scrDir}

exit 0
EOF_STAGE

chmod 777 ${stageIN} ${run_pmo}

TRN_JOB_ID=`${BATCH_DEPEND} ${stageIN}`
PMO_JOB_ID=`${BATCH} -W depend=afterok:${TRN_JOB_ID} ${run_pmo}`

if [ ! -e ${queue_log} ]; then touch ${queue_log}; fi
cat  >> ${queue_log} << EOF_QUEUE
${EXPNAME} ${dtg} RUN_PMO
  Submitted at `date`
  TRN_JOB_ID : ${TRN_JOB_ID}
  PMO_JOB_ID : ${PMO_JOB_ID}
EOF_QUEUE

rm -f ${run_pmo} ${stageIN}

echo ${PMO_JOB_ID}

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

