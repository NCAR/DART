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
#
usage="Usage: `basename $0` [-d dtg -i icycle -n init_ens -c cycle_da]" 
while getopts ":d:i:n:m:W:c:" option
do
  case $option in
      d  ) dtg=$OPTARG;;
      i  ) icycle_loc=$OPTARG;;
      n  ) init_ens=$OPTARG;;
      m  ) mems_to_run=$OPTARG;;
      c  ) CYCLE_DA=$OPTARG;;
      W  ) depend=$OPTARG;;
      \? ) echo $usage ; exit 65;;
      *  ) echo $usage ; exit 65;;
  esac
done

echo "#"
echo "#################### BEGIN advance_wrapper.csh ####################"
echo "#"

BATCH="qsub"
SED='sed -i -e'

if [[ $depend ]]; then depend="-W depend=afterok:${depend}"; fi
: ${BATCH_DEPEND:=$BATCH}

PATH_CONFIG=CONFIG/FILE/GOES/HERE
. ${PATH_CONFIG}

: ${dtg:=$cdtg_beg}
: ${icycle_loc:=$icycle}
: ${init_ens:='f'}
: ${was_first_assim:='f'}
: ${mems_to_run:=-1}
: ${CYCLE_DA:='t'}
: ${is_fcp_bndy:='f'}
CYCLE_DA=${CYCLE_DA:0:1}

if [ ${init_ens} == 'f' ]; then 
  is_first_ens='false'
  WALLTIME_COAMM_LOC=${WALLTIME_COAMM}
else
  is_first_ens='true'
  WALLTIME_COAMM_LOC=${WALLTIME_COAMM_INIT}
fi

: ${scriptDIAG:=$HOME/diagnostic_cdf}

scrDir=${COAMPS_RUN_BASE_DIR}
LOG_DIR=${scrDir}/log
LOG_DTG=${LOG_DIR}/${dtg}
scriptDir=${scrDir}/scripts
DART_BIN=${scrDir}/bin

advance_group=${scriptDir}/advance_group.${dtg}.sh
resubmitDA=${scriptDir}/resub_enkf.${dtg}.sh
fcst_group=${scriptDir}/fcst_group.${dtg}.sh
queue_log=${LOG_DIR}/queue_${EXPNAME}.log

cd ${scrDir}
nnest=`grep nnest ${scrDir}/namelist | sed 's/,//g' | awk -F= '{print $2}'`
dtgp1=`echo ${dtg} ${icycle_loc} | ${DART_BIN}/advance_time`
dtgm1=`echo ${dtg} -${icycle_loc} | ${DART_BIN}/advance_time`
cdtg_beg_p1=`echo ${cdtg_beg}  ${fcstlen_init} | ${DART_BIN}/advance_time`
if [ ${dtg} -eq ${cdtg_beg_p1} ]; then was_first_assim='t'; fi

if [ ! -d  ${LOG_DIR} ]; then mkdir -p ${LOG_DIR}; fi
if [ ! -d  ${LOG_DTG} ]; then mkdir -p ${LOG_DTG}; fi

# Get the machine layout
. ${scriptDir}/HPC_CONFIG.sh
MACH_LAYOUT
let "COAMM_NODES = (COAMM_PROCS / CPU_PER_NODE) + (COAMM_PROCS % CPU_PER_NODE)"

mems_to_run=( `echo ${mems_to_run} | sed 's/,/ /g'` )
if [ ${JOB_ARRAY_AVAILABLE} == 'n' -a ${mems_to_run} -eq -1 ]; then
  mems_to_run=( `seq 1 ${ENSEMBLE_SIZE}` )
fi

for member in ${mems_to_run[@]}; do

if [ ${JOB_ARRAY_AVAILABLE} == 'y' -a ${mems_to_run} -eq -1 ]; then

if [ -e ${advance_group} ]; then rm ${advance_group}; fi
advance_group_array=${advance_group}
cat > ${advance_group} << EOF
#!/bin/bash
`PBS_DIRECTIVE ${EXPNAME} ${dtg} ${COAMM_PROCS} ${COAMM_NODES} ${WALLTIME_COAMM_LOC} ${LOG_DTG} n y ${PBS_QUEUE_ENS}` 
EOF
${SED} "s/JOB_ARRAY/1\-${ENSEMBLE_SIZE}/" ${advance_group}

else

advance_group_tmp=`echo ${advance_group} | sed "s/advance_group/advance_group.${member}/"`
advance_group_array=( ${advance_group_array[@]} ${advance_group_tmp} )

if [ -e ${advance_group_tmp} ]; then rm ${advance_group_tmp}; fi
cat > ${advance_group_tmp} << EOF
#!/bin/bash
`PBS_DIRECTIVE M${member}.${EXPNAME} ${dtg} ${COAMM_PROCS} ${COAMM_NODES} ${WALLTIME_COAMM_LOC} ${LOG_DTG} n n` 
PBS_ARRAY_INDEX=${member}
EOF
fi

done

for advance_group_script in ${advance_group_array[@]}; do

cat >> ${advance_group_script} << EOF
SED='${SED}'

scrDir=${scrDir}
WORKDIR=${scrDir}/${dtg}
WORKDIR_NEW=${scrDir}/${dtgp1}
LOG_DIR=\${scrDir}/log
scriptDir=\${scrDir}/scripts
DART_BIN=\${scrDir}/bin
scriptDIAG=${scriptDIAG}

COAMPS_BIN_DIR=${COAMPS_HOME}/bin
COAMPS_DART_DIR=${DART_HOME}
NL_STRIP=\${scriptDir}/strip_namelist.pl

echo Entering ${advance_group_script}
echo   Running on host \`hostname\`
echo   Time is \`date\`
echo   Directory is \`pwd\`
echo   Working directory is: \$WORKDIR

MPI_OCNA="${COAMM_MPI_CMD}"
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
COAMPS_ENS=\${scrDir}/\${ensdir}
is_first_ens=${is_first_ens}
is_fcp_bndy=${is_fcp_bndy}
was_first_assim=${was_first_assim}
icycle_loc=${icycle_loc}
fcstlen_init=${fcstlen_init}

if [ \$is_first_ens == false ]; then
  cd \${COAMPS_ENS}/data
  ic_file=\`printf coamps_state_anal.%04d \${element}\`
  cp \${WORKDIR}/\${ic_file} \${COAMPS_ENS}/data/dart_vector

  if [ \${was_first_assim} == t ]; then
    \${SED} "s/\(icycle\s*=\s*\).*/\1\${fcstlen_init},/" ../namelist
  fi
  \${NL_STRIP} ../namelist convert.vars convert
  \${SED} 's/\(^&convert.*\)/\1\n  is_pmo = \.false\.,/' ./convert.nml
  \${SED} "s/\(^&convert.*\)/\1\n  is_first = \.\${is_first_ens}\.,/" ./convert.nml

  # Convert the dart initial condition file to a COAMPS restart file
  rm -f dart.kstart
  \${DART_BIN}/trans_dart_to_coamps
  if [ \$? != 0 ]; then
    echo "Unsuccessful translation from DART of member \$element!"
    exit 1
  fi
fi

logtime=$dtg
cd \${COAMPS_ENS}

# Run ncoda ocean analysis
echo "Running ocean analysis for \$element"
\${COAMPS_BIN_DIR}/ncoda_prep 2D namelist ${dtg} > log.o.\${logtime}
mpi_status=\$?
if  [ \$mpi_status -ne 0 ]; then
    echo "Unable to run ncoda_prep for member \$element!"
    exit 1
fi

\${MPI_OCNA} \${COAMPS_BIN_DIR}/ncoda 2D namelist ${dtg} >> log.o.\${logtime}
if  [ \$mpi_status -ne 0 ]; then
    echo "Unable to run ncoda for member \$element!"
    exit 1
fi

\${COAMPS_BIN_DIR}/ncoda_post 2D namelist ${dtg} >> log.o.\${logtime}
if  [ \$mpi_status -ne 0 ]; then
    echo "Unable to run ncoda_post for member \$element!"
    exit 1
fi

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
  cd \${COAMPS_ENS}/data
  if [ \$is_first_ens == 'true' ]; then
    \${NL_STRIP} ../namelist ./perturb.vars pert_init
    \${SED} 's|\(^&pert_init.*\)|\1\n  ens_size=$ENSEMBLE_SIZE,|' ./pert_init.nml
    \${SED} "s|\(^&pert_init.*\)|\1\n  dsnrff_pert=\'$INIT_PERT_DIR\',|" ./pert_init.nml
    echo \$element | \${DART_BIN}/perturb_init
    status=\$?
    if  [ \$status -ne 0 ]; then
      echo "Unable to perturb initial conditions for member \$element!"
      exit 1
    fi
  fi

  \${NL_STRIP} ../namelist ./perturb.vars pert_bndy
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

# Copy the current data header file to the previous dtg
# tricking COAMPS to think this is not a cold start
if [ \${was_first_assim} == t ]; then
  \${SED} "s/\(icycle\s*=\s*\).*/\1\${icycle_loc},/" ./namelist
  datahd_cur=\`ls ./data/datahd*1a*${dtg}*\` 
  datahd_prev=\`echo \${datahd_cur} | sed "s/${dtg}/${dtgm1}/"\`
  cp \${datahd_cur} \${datahd_prev}   
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
  echo \${scriptDIAG}/diag_pv.sh ${dtg} \${inest} -name=\${logtime}.\${element_num} -indir=\${COAMPS_ENS}/data -odir=\${COAMPS_ENS}/data -begt=0 -endt=\${icycle_loc} -intt=\${icycle_loc}
  \${scriptDIAG}/diag_pv.sh ${dtg} \${inest} -name=\${logtime}.\${element_num} -indir=\${COAMPS_ENS}/data -odir=\${COAMPS_ENS}/data -begt=0 -endt=\${icycle_loc} -intt=\${icycle_loc}
done

# Now that the we've run the model, move into the COAMPS data directory 
# and move things back to the DART directory. 
cd \${COAMPS_ENS}/data
rm -f dart_vector

# Convert the resulting COAMPS file to a dart file
\${NL_STRIP} ../namelist convert.vars convert
\${SED} 's/\(^&convert.*\)/\1\n  is_pmo = \.false\.,/' ./convert.nml
#\${SED} "s/\(^&convert.*\)/\1\n  is_first = \.\${is_first_ens}\.,/" ./convert.nml
\${SED} "s/\(^&convert.*\)/\1\n  is_first = \.true\.,/" ./convert.nml
if [ \${is_first_ens} == true ]; then
  \${SED} "s/\(icycle\s*=\s*\).*/\1\${icycle_loc},/" ./convert.nml
fi
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

if [[ \${ic_file} ]]; then rm -f \${WORKDIR}/\${ic_file} ; fi

cd \${COAMPS_ENS}

echo Exiting ${advance_group_script}
echo   Time is \`date\`

exit 0
EOF

done

if [ -e ${resubmitDA} ]; then rm ${resubmitDA}; fi

if [ ${PBS_QUEUE:0:1} = 'R' ]; then

cat > ${resubmitDA} << EOF
`PBS_DIRECTIVE RESUB_ENKF ${dtg} ${CPU_PER_NODE} 1 00:02:00 ${LOG_DTG} n n`
${scriptDir}/run_filter.sh -d ${dtgp1}
exit 0
EOF

else

cat > ${resubmitDA} << EOF
#!/bin/bash
`PBS_DIRECTIVE RESUB_ENKF ${dtg} ${TRANS_PROCS} ${MIN_PROCS} 00:02:00 ${LOG_DIR} y n` 
${scriptDir}/run_filter.sh -d ${dtgp1}
exit 0
EOF

fi

echo "Submitting jobs to PBS queue..."
chmod 777 ${resubmitDA} ${advance_group_array[@]}

for advance_group_script in ${advance_group_array[@]}; do
  GROUP_JOB_ID=( ${GROUP_JOB_ID[@]} `${BATCH_DEPEND} ${depend} ${advance_group_script} | sed 's/\..*$//'` )
done
GROUP_JOB_ID=`echo ${GROUP_JOB_ID[@]} | sed 's/ /:/g'`

cd ${scrDir}
if [ $CYCLE_DA == "t" -o $CYCLE_DA == "T" ]; then
  # RESUBMIT FILTER JOB
  if [ $dtgp1 -le  $cdtg_end ]; then
     RSUB_JOB_ID=`${BATCH} -W depend=afterok:${GROUP_JOB_ID} ${resubmitDA}`
  fi

  # Submit archive job.
  TRNS_JOB_ID=`${scriptDir}/archive_coamps_ens.sh -d ${dtg} -W ${GROUP_JOB_ID}`

  # Submit forecasting job.
  #FCST_JOB_ID=`qsub -W depend=afterok:${GROUP_JOB_ID} ${fcst_group}`

  if [ ${init_ens} == 'f' ]; then 
    INCR_JOB_ID=`${scriptDir}/create_increment.sh -d ${dtg} -W ${GROUP_JOB_ID} -n ${nnest}`
  fi
  MEAN_JOB_ID=`${scriptDir}/create_mean_std.sh -d ${dtg} -W ${GROUP_JOB_ID} -n ${nnest}`
fi

if [ ! -e ${queue_log} ]; then touch ${queue_log}; fi
cat  >> ${queue_log} << EOF_QUEUE
${EXPNAME} ${dtg} ADVANCE_WRAPPER
  Submitted at `date`
  CYCL_JOB_ID: ${GROUP_JOB_ID}
  TRNS_JOB_ID: ${TRNS_JOB_ID}
  INCR_JOB_ID: ${INCR_JOB_ID}
  MEAN_JOB_ID: ${MEAN_JOB_ID}
  RSUB_JOB_ID: ${RSUB_JOB_ID}
  FCST_JOB_ID: ${FCST_JOB_ID}
EOF_QUEUE

rm -f ${advance_group_array[@]}
rm -f ${resubmitDA}
rm -f ${fcst_group}

echo "#################### END advance_wrapper.sh ####################"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

