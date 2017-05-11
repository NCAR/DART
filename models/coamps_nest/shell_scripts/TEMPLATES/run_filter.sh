#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
###############################################################################
#
# AUTHOR:   P. A. Reinecke
#           Naval Research Laboratory
#
###############################################################################
#
usage="Usage: `basename $0` -d dtg -i icycle -W depend" 
while getopts ":d:i:W:" option
do
  case $option in
      d  ) dtg=$OPTARG;;
      i  ) icycle_loc=$OPTARG;;
      W  ) depend=$OPTARG;;
      \? ) echo $usage ; exit 65;;
      *  ) echo $usage ; exit 65;;
  esac
done
BATCH="qsub"
SED='sed -i -e'

#if [[ $depend ]]; then BATCH_DEPEND="${BATCH} -W depend=afterok:${depend}"; fi
if [[ $depend ]]; then depend="-W depend=afterok:${depend}"; fi
: ${BATCH_DEPEND:=$BATCH}

PATH_CONFIG=CONFIG/FILE/GOES/HERE
. ${PATH_CONFIG}
: ${dtg:=$cdtg_beg}
: ${init_ens:='f'}
: ${init_assim:='f'}

scrDir=${COAMPS_RUN_BASE_DIR}
LOG_DIR=${scrDir}/log
scriptDir=${scrDir}/scripts
DART_BIN=${scrDir}/bin

# Test if this is the initial ensemble
if [ ${dtg} == ${cdtg_beg} ]; then init_ens='t'; fi
if [ ${init_ens} == 'f' ]; then
: ${icycle_loc:=$icycle}
else
: ${icycle_loc:=$fcstlen_init}
fi

# Test if this is the initial data assimilation
cd ${scrDir}
cdtg_beg_p1=`echo ${cdtg_beg} ${fcstlen_init} | ${DART_BIN}/advance_time`
if [ ${dtg} == ${cdtg_beg_p1} ]; then init_assim='t'; fi
if [ ${init_assim} == 'f' ]; then
 icycle_prev=${icycle}
else
 icycle_prev=${fcstlen_init}
fi

stageIN=${scriptDir}/stageIN_ens.${dtg}.sh
run_filter=${scriptDir}/run_filter.${dtg}.sh
queue_log=${LOG_DIR}/queue_${EXPNAME}.log

. ${scriptDir}/HPC_CONFIG.sh
MACH_LAYOUT
let "FILTER_NODES = (FILTER_PROCS / CPU_PER_NODE) + (FILTER_PROCS % CPU_PER_NODE)"

if [ -e ${run_filter} ]; then rm ${run_filter}; fi
if [ ${init_ens} == 't' ]; then
cat > ${run_filter} << EOF
#!/bin/bash
`PBS_DIRECTIVE INIT_FILTER ${dtg} ${TRANS_PROCS} ${MIN_PROCS} ${WALLTIME_FILTER} ${LOG_DIR} y n` 
EOF
else
cat > ${run_filter} << EOF
#!/bin/bash
`PBS_DIRECTIVE RUN_FILTER ${dtg} ${FILTER_PROCS} ${FILTER_NODES} ${WALLTIME_FILTER} ${LOG_DIR} n n` 
EOF
fi

cat >> ${run_filter} << EOF
#

# Set how many processors we will use for filter
MPI_CMD="${FILTER_MPI_CMD}"

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
#
dtg=${dtg}
init_ens=${init_ens}
is_real=${IS_REAL_DATA}
#
WORKDIR=\${WORKDIR}/\${dtg}
if [ ! -d \${WORKDIR} ]; then mkdir -p \${WORKDIR}; fi
#
scrDir=${scrDir}
DART_BIN=\${scrDir}/bin
LOG_DIR=\${scrDir}/log
scriptDir=\${scrDir}/scripts
#
cd \${WORKDIR}
#
echo "#################### BEGIN run_filter.sh ####################"
echo \`date\`

# Update the namelists to the new date-time group.
\${scriptDir}/update_namelists.sh -n \${init_ens} -d ${dtg} -i ${icycle_loc}

if [ \${init_ens} == 'f' ]; then
 
  # Real Data Simulations
  if [ \${is_real} == 'T' ]; then
    # Preprocess observations with navdas obs preprocessor 
    \${scriptDir}/navdas_preproc_obs.sh -d ${dtg} -m 1 -i ${icycle_prev}
    if [ \$? -ne 0 ]; then
      echo "Error in \${scriptDir}/navdas_preproc_obs.sh"
      exit 1
    fi

    # Translate the navdas innovations file to a dart obs_seq file
	\${DART_BIN}/innov_to_obs_seq
    if [ \$? -ne 0 ]; then
      echo "Error in \${DART_BIN}/innov_to_obs_seq"
      exit 1
    fi

  fi
  \${MPI_CMD} \${DART_BIN}/filter | tee filter.log
  if [ \$? -ne 0 ]; then
    echo "Error in \${DART_BIN}/filter"
    exit 1
  fi

# change the names of the prior and posterior netCDF files.
  mv ./preassim.nc ./preassim.\${dtg}.nc
  mv ./analysis.nc ./analysis.\${dtg}.nc

  if [ ! -e \${scrDir}/obs_seq.list ]; then touch \${scrDir}/obs_seq.list; fi
  echo "./\${dtg}/obs_seq.final" >> \${scrDir}/obs_seq.list

fi

\${scriptDir}/advance_wrapper.sh -d ${dtg} -i ${icycle_loc} -n \${init_ens}

echo \`date\`
echo "#################### END run_filter.sh ####################"

exit 0
EOF

if [ -e $stageIN ]; then rm ${stageIN}; fi
cat > ${stageIN} << EOF_STAGE
#!/bin/bash
`PBS_DIRECTIVE STAGE_IN_ENS ${dtg} ${TRANS_PROCS} ${MIN_PROCS} ${WALLTIME_TRANS} ${LOG_DIR} y n` 
EOF_STAGE

cat >> ${stageIN} << EOF_STAGE
NOGAPS_ARC=${NOGAPS_ARC}
ARCMACH=${ARCMACH}

NOGAPS_PUT=${NOGAPS_PATH}/e${dtg}
file1=${dtg:2:8}.outp.tar

if [ ! -e \${NOGAPS_PUT} ]; then mkdir -p \${NOGAPS_PUT}; fi
cd \${NOGAPS_PUT}

rsh ${ARCMACH} "stage -x \${NOGAPS_ARC}/\${file1}"
rsh ${ARCMACH} "stage \${NOGAPS_ARC}/\${file1}"
rcp ${ARCMACH}:\${NOGAPS_ARC}/\${file1} .
(tar xf \${file1} ; rm -f \${file1})

for mem in \`seq 1 ${ENSEMBLE_SIZE}\`; do
  file2=\`printf "outp%03d.tar" \${mem}\`
  (tar xf \${file2} ; rm -f \${file2})
done

exit 0
EOF_STAGE

chmod 777 ${stageIN} ${run_filter}

get_nogaps='f'
if [ $get_nogaps == 't' ]; then
  TRNS_JOB_ID=`${BATCH_DEPEND} ${stageIN}`
  FLTR_JOB_ID=`${BATCH} -W depend=afterok:${TRNS_JOB_ID} ${run_filter}`
else
  TRNS_JOB_ID=TRNS_JOB_ID
  if [ ${init_ens} == 'f' ]; then
#FLTR_JOB_ID=`${BATCH_DEPEND} ${run_filter}`
  FLTR_JOB_ID=`${BATCH} ${run_filter} ${depend}`
  fi
fi
#if [ ${init_ens} == 'f' ]; then
#   ${scriptDir}/advance_wrapper.sh -d ${dtg} -i ${icycle_loc} -n ${init_ens} -W ${FLTR_JOB_ID}
#else
#   ${run_filter}
#fi

if [ ${init_ens} != 'f' ]; then
   ${run_filter}
fi

if [ ! -e ${queue_log} ]; then touch ${queue_log}; fi
cat  >> ${queue_log} << EOF_QUEUE
${EXPNAME} ${dtg} RUN_FILTER
  Submitted at `date`
  TRNS_JOB_ID: ${TRNS_JOB_ID}
  FLTR_JOB_ID: ${FLTR_JOB_ID}
EOF_QUEUE

rm -f ${stageIN} 
rm -f ${run_filter}

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

