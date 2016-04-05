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
######
#
usage="Usage: `basename $0` [-d dtg -i icycle -n init_ens]" 
while getopts ":d:i:n:m:W:" option
do
  case $option in
      d  ) dtg=$OPTARG;;
      i  ) icycle_loc=$OPTARG;;
      n  ) init_ens=$OPTARG;;
      m  ) mems_to_run=$OPTARG;;
      W  ) depend=$OPTARG;;
      \? ) echo $usage ; exit 65;;
      *  ) echo $usage ; exit 65;;
  esac
done

echo "#"
echo "#################### BEGIN forecast_wrapper.csh ####################"
echo "#"

BATCH="qsub"
SED='sed -i -e'

if [[ $depend ]]; then depend="-W depend=afterok:${depend}"; fi
: ${BATCH_DEPEND:=$BATCH}

PATH_CONFIG=/usr/local/u/reinecke/RUN_DART/ALPS_NOV2007_32MEM_REAL/paths.config
. ${PATH_CONFIG}

: ${dtg:=$cdtg_beg}
: ${icycle_loc:=$icycle}
: ${mems_to_run:=-1}
: ${is_fcp_bndy:='f'}

WALLTIME_COAMM_LOC=${WALLTIME_FCST}

scrDir=${COAMPS_RUN_BASE_DIR}
LOG_DIR=${scrDir}/log
LOG_DTG=${LOG_DIR}/${dtg}
scriptDir=${scrDir}/scripts
: ${scriptDIAG:=$HOME/diagnostic_cdf}

forecast_group=${scriptDir}/forecast_group.${dtg}.sh

cd ${scrDir}
nnest=`grep nnest ${scrDir}/namelist | sed 's/,//g' | awk -F= '{print $2}'`

if [ ! -d  ${LOG_DIR} ]; then mkdir -p ${LOG_DIR}; fi
if [ ! -d  ${LOG_DTG} ]; then mkdir -p ${LOG_DTG}; fi

# Get the machine layout
. ${scriptDir}/HPC_CONFIG.sh
MACH_LAYOUT
let "COAMM_NODES = (COAMM_PROCS / CPU_PER_NODE) + (COAMM_PROCS % CPU_PER_NODE)"

mems_to_run_in=${mems_to_run}
mems_to_run=( `echo ${mems_to_run} | sed 's/,/ /g'` )
if [ $mems_to_run -eq -1 ]; then whole_ens='t'; else whole_ens='f'; fi 

if [ ${JOB_ARRAY_AVAILABLE} == 'n' -a ${mems_to_run} -eq -1 ]; then
  mems_to_run=( `seq 1 ${ENSEMBLE_SIZE}` )
fi

for member in ${mems_to_run[@]}; do

if [ ${JOB_ARRAY_AVAILABLE} == 'y' -a ${mems_to_run} -eq -1 ]; then

if [ -e ${forecast_group} ]; then rm ${forecast_group}; fi
forecast_group_array=${forecast_group}
cat > ${forecast_group} << EOF
#!/bin/bash
`PBS_DIRECTIVE ${EXPNAME} ${dtg} ${COAMM_PROCS} ${COAMM_NODES} ${WALLTIME_COAMM_LOC} ${LOG_DTG} n y ${PBS_QUEUE_ENS}` 
EOF
${SED} "s/JOB_ARRAY/1\-${ENSEMBLE_SIZE}/" ${forecast_group}

else

forecast_group_tmp=`echo ${forecast_group} | sed "s/forecast_group/forecast_group.${member}/"`
forecast_group_array=( ${forecast_group_array[@]} ${forecast_group_tmp} )

if [ -e ${forecast_group_tmp} ]; then rm ${forecast_group_tmp}; fi
cat > ${forecast_group_tmp} << EOF
#!/bin/bash
`PBS_DIRECTIVE M${member}.fcst.${EXPNAME} ${dtg} ${COAMM_PROCS} ${COAMM_NODES} ${WALLTIME_COAMM_LOC} ${LOG_DTG} n n` 
PBS_ARRAY_INDEX=${member}
EOF
fi

done

for forecast_group_script in ${forecast_group_array[@]}; do

cat >> ${forecast_group_script} << EOF
SED='${SED}'

scrDir=${scrDir}
WORKDIR=${scrDir}/${dtg}
LOG_DIR=\${scrDir}/log
scriptDir=\${scrDir}/scripts
COAMPS_BIN_DIR=${COAMPS_HOME}/bin
DART_BIN=\${scrDir}/bin
NL_STRIP=\${scriptDir}/strip_namelist.pl
scriptDIAG=${scriptDIAG}

echo Entering ${forecast_group_script}
echo   Running on host \`hostname\`
echo   Time is \`date\`
echo   Directory is \`pwd\`
echo   Working directory is: \$WORKDIR

MPI_FCST="${COAMM_MPI_CMD}"
MPI_BNDY="${COAMA_MPI_CMD}"

# The COAMPS ensemble configuration has directories associated with
# the ensemble member number, thus member numbers with digits less
# than the maximum must be zero-padded.  For ease of use, just force
# this to a fixed digit length for now.
element=\${PBS_ARRAY_INDEX}
digits=5
element_num=\`printf "%0\${digits}d" \$element\`
ensdir=\$element_num
is_fcp_bndy=${is_fcp_bndy}
COAMPS_ENS=\${scrDir}/\${ensdir}/fcst.${dtg}
icycle_loc=${icycle_loc}
fcstlen=${fcstlen}

logtime=$dtg

cd \${COAMPS_ENS}
# Get the boundary conditions
echo "Getting boundary conditions for \$element"
\${MPI_BNDY} \${COAMPS_BIN_DIR}/atmos_analysis.exe > log.a.\${logtime}
mpi_status=\$?
if  [ \$mpi_status -ne 0 ]; then
    echo "Unable to get boundaries \$element!"
    exit 1
fi

if [ \$is_fcp_bndy == 'true' ]; then
  cd \${scrDir}/\${ensdir}/data

  \${NL_STRIP} \${COAMPS_ENS}/namelist ./perturb.vars pert_bndy
  \${SED} 's|\(^&pert_bndy.*\)|\1\n  ens_size=$ENSEMBLE_SIZE,|' ./pert_bndy.nml
  \${SED} 's|\(^&pert_bndy.*\)|\1\n  alpha=$alpha_fcp,|' ./pert_bndy.nml
  \${SED} 's|\(^&pert_bndy.*\)|\1\n  npert_tot=$nbndy_perts,|' ./pert_bndy.nml
  \${SED} "s|\(^&pert_bndy.*\)|\1\n  dsnrff_pert=\'$BNDY_PERT_DIR\',|" ./pert_bndy.nml

  echo \$element | \${DART_BIN}/perturb_bndy
  status=\$?
  if  [ \$status -ne 0 ]; then
    echo "Unable to perturb boundaries for member \$element!"
    exit 1
  fi
  cd \${COAMPS_ENS}

fi

# Run the COAMPS forecast
echo "Running COAMPS ensemble member \$element"
\${MPI_FCST} \${COAMPS_BIN_DIR}/atmos_forecast.exe > log.m.\${logtime}
mpi_status=\$?
if [ \$mpi_status -ne 0 ]; then
    echo "Unsuccessful completion of member \$element!"
    exit 1
fi

for inest in \`seq 1 $nnest\`; do
  echo \${scriptDIAG}/diag_pv.sh ${dtg} \${inest} -name=\${logtime}.\${element_num} -indir=\${COAMPS_ENS}/data -odir=\${COAMPS_ENS}/data -begt=0 -endt=\${fcstlen} -intt=\${icycle_loc}
  \${scriptDIAG}/diag_pv.sh ${dtg} \${inest} -name=\${logtime}.\${element_num} -indir=\${COAMPS_ENS}/data -odir=\${COAMPS_ENS}/data -begt=0 -endt=\${fcstlen} -intt=\${icycle_loc}
done

echo Exiting ${forecast_group_script}
echo   Time is \`date\`

exit 0
EOF

done

echo "Submitting jobs to PBS queue..."
chmod 777 ${forecast_group_array[@]}

if [ ${whole_ens} == 't' ]; then
  ${scriptDir}/update_namelists.sh -d ${dtg} -f t -i ${icycle_loc}
else
  ${scriptDir}/update_namelists.sh -d ${dtg} -f t -i ${icycle_loc} -m ${mems_to_run_in}
fi

#for forecast_group_script in ${forecast_group_array[@]}; do
for forecast_group_script in ${forecast_group_array[@]}; do
  GROUP_JOB_ID=( ${GROUP_JOB_ID[@]} `${BATCH_DEPEND} ${forecast_group_script} ${depend} | sed 's/\..*$//'` )
done
GROUP_JOB_ID=`echo ${GROUP_JOB_ID[@]} | sed 's/ /:/g'`

cd ${scrDir}
MEAN_JOB_ID=`${scriptDir}/create_mean_std.sh -d ${dtg} -W ${GROUP_JOB_ID} -f t -n ${nnest}`


rm -f ${forecast_group_array[@]}

echo "#################### END forecast_wrapper.sh ####################"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

