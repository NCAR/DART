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
usage="Usage: `basename $0` [-d dtg -i icycle -m mems_to_run -W dependency]" 
while getopts ":d:i:m:W:" option
do
  case $option in
      d  ) dtg=$OPTARG;;
      i  ) icycle_loc=$OPTARG;;
      m  ) mems_to_run=$OPTARG;;
      W  ) depend=$OPTARG;;
      \? ) echo $usage ; exit 65;;
      *  ) echo $usage ; exit 65;;
  esac
done

echo "#"
echo "#################### BEGIN restart_wrapper.csh ####################"
echo "#"

BATCH="qsub"
SED='sed -i -e'

if [[ $depend ]]; then depend="-W depend=afterok:${depend}"; fi
: ${BATCH_DEPEND:=$BATCH}

PATH_CONFIG=CONFIG/FILE/GOES/HERE
. ${PATH_CONFIG}

: ${dtg:=$cdtg_beg}
: ${icycle_loc:=$icycle}
: ${is_restart:='f'}
: ${mems_to_run:=-1}
: ${CYCLE_DA:='t'}
CYCLE_DA=${CYCLE_DA:0:1}

: ${WALLTIME_RESTART:="00:02:00"}
WALLTIME_RESTART_LOC=${WALLTIME_RESTART}

scrDir=${COAMPS_RUN_BASE_DIR}
DART_BIN=${scrDir}/bin
scriptDir=${scrDir}/scripts

dtgp1=${dtg}
dtg=`echo ${dtgp1} -${icycle_loc} | ${DART_BIN}/advance_time`
dtgm1=`echo ${dtg} -${icycle_loc} | ${DART_BIN}/advance_time`

LOG_DIR=${scrDir}/log
LOG_DTG=${LOG_DIR}/${dtg}

restart_group=${scriptDir}/restart_group.${dtg}.sh
resubmitDA=${scriptDir}/resub_enkf.${dtg}.sh
queue_log=${LOG_DIR}/queue_${EXPNAME}.log

cd ${scrDir}
nnest=`grep nnest ${scrDir}/namelist | sed 's/,//g' | awk -F= '{print $2}'`

if [ ! -d  ${LOG_DIR} ]; then mkdir -p ${LOG_DIR}; fi
if [ ! -d  ${LOG_DTG} ]; then mkdir -p ${LOG_DTG}; fi

# Get the machine layout
. ${scriptDir}/HPC_CONFIG.sh
MACH_LAYOUT

if [ ${mems_to_run} -eq -1 ]; then
  ${scriptDir}/update_namelists.sh -d ${dtg}
else
  ${scriptDir}/update_namelists.sh -d ${dtg} -m ${mems_to_run}
fi

mems_to_run=( `echo ${mems_to_run} | sed 's/,/ /g'` )
if [ ${JOB_ARRAY_AVAILABLE} == 'n' -a ${mems_to_run} -eq -1 ]; then
  mems_to_run=( `seq 1 ${ENSEMBLE_SIZE}` )
fi

for member in ${mems_to_run[@]}; do

if [ ${JOB_ARRAY_AVAILABLE} == 'y' -a ${mems_to_run} -eq -1 ]; then

if [ -e ${restart_group} ]; then rm ${restart_group}; fi
restart_group_array=${restart_group}
cat > ${restart_group} << EOF
#!/bin/bash
`PBS_DIRECTIVE ${EXPNAME} ${dtg} ${CPU_PER_NODE} 1 ${WALLTIME_RESTART_LOC} ${LOG_DTG} n y ${PBS_QUEUE_ENS}` 
EOF
${SED} "s/JOB_ARRAY/1\-${ENSEMBLE_SIZE}/" ${restart_group}

else

restart_group_tmp=`echo ${restart_group} | sed "s/restart_group/restart_group.${member}/"`
restart_group_array=( ${restart_group_array[@]} ${restart_group_tmp} )

if [ -e ${restart_group_tmp} ]; then rm ${restart_group_tmp}; fi
cat > ${restart_group_tmp} << EOF
#!/bin/bash
`PBS_DIRECTIVE ${EXPNAME}.${member} ${dtg} ${CPU_PER_NODE} 1 ${WALLTIME_RESTART_LOC} ${LOG_DTG} n n` 
PBS_ARRAY_INDEX=${member}
EOF
fi

done

for restart_group_script in ${restart_group_array[@]}; do

cat >> ${restart_group_script} << EOF
SED='${SED}'

scrDir=${scrDir}
WORKDIR_NEW=${scrDir}/${dtgp1}
LOG_DIR=\${scrDir}/log
scriptDir=\${scrDir}/scripts
DART_BIN=\${scrDir}/bin

COAMPS_BIN_DIR=${COAMPS_HOME}/bin
COAMPS_DART_DIR=${DART_HOME}
NL_STRIP=\${scriptDir}/strip_namelist.pl

echo Entering ${restart_group_script}
echo   Running on host \`hostname\`
echo   Time is \`date\`
echo   Directory is \`pwd\`

# The COAMPS ensemble configuration has directories associated with
# the ensemble member number, thus member numbers with digits less
# than the maximum must be zero-padded.  For ease of use, just force
# this to a fixed digit length for now.
element=\${PBS_ARRAY_INDEX}
digits=5
element_num=\`printf "%0\${digits}d" \$element\`
ensdir=\$element_num
COAMPS_ENS=\${scrDir}/\${ensdir}
icycle_loc=${icycle_loc}

logtime=$dtg
cd \${COAMPS_ENS}/data
rm -f dart_vector

# Convert the resulting COAMPS file to a dart file
\${NL_STRIP} ../namelist convert.vars convert
\${SED} 's/\(^&convert.*\)/\1\n  is_pmo = \.false\.,/' ./convert.nml
\${SED} "s/\(^&convert.*\)/\1\n  is_first = \.true\.,/" ./convert.nml

\${DART_BIN}/trans_coamps_to_dart
if [ \$? -ne 0 ]; then
    echo "Unsuccessful translation to DART of member \$element!"
    exit 1
fi

if [ ! -e \${WORKDIR_NEW} ]; then mkdir -p \${WORKDIR_NEW}; fi
# Move the dart file back to the working directory -  only now can we
# safely delete the initial condition file.
ud_file=\`printf coamps_state_fcst.%04d \${element}\`
echo "mv dart_vector \${WORKDIR_NEW}/\${ud_file}"
mv dart_vector \${WORKDIR_NEW}/\${ud_file}
if [ \${element} -eq 1 ]; then
  cp datahd*_${dtg}_* \${WORKDIR_NEW}
  cp terrht*_${dtg}_* \${WORKDIR_NEW}
fi

cd \${COAMPS_ENS}

echo Exiting ${restart_group_script}
echo   Time is \`date\`

exit 0
EOF

done

if [ -e ${resubmitDA} ]; then rm ${resubmitDA}; fi
cat > ${resubmitDA} << EOF
#!/bin/bash
`PBS_DIRECTIVE RESUB_ENKF ${dtg} ${TRANS_PROCS} ${MIN_PROCS} 00:01:00 ${LOG_DIR} y n` 
${scriptDir}/run_filter.sh -d ${dtgp1}
exit 0
EOF

echo "Submitting jobs to PBS queue..."
chmod 777 ${resubmitDA} ${restart_group_array[@]}

for restart_group_script in ${restart_group_array[@]}; do
  GROUP_JOB_ID=( ${GROUP_JOB_ID[@]} `${BATCH_DEPEND} ${restart_group_script} ${depend} | sed 's/\..*$//'` )
done
GROUP_JOB_ID=`echo ${GROUP_JOB_ID[@]} | sed 's/ /:/g'`

cd ${scrDir}
if [ $dtgp1 -le  $cdtg_end ]; then
  RSUB_JOB_ID=`${BATCH} -W depend=afterok:${GROUP_JOB_ID} ${resubmitDA}`
fi

if [ ! -e ${queue_log} ]; then touch ${queue_log}; fi
cat  >> ${queue_log} << EOF_QUEUE
${EXPNAME} ${dtg} RESTART_WRAPPER
  Submitted at `date`
  CYCL_JOB_ID: ${GROUP_JOB_ID}
  RSUB_JOB_ID: ${RSUB_JOB_ID}
EOF_QUEUE

rm -f ${restart_group_array[@]}
rm -f ${resubmitDA}

echo "#################### END restart_wrapper.sh ####################"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

