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
# Runs the program that will calculate the incr and standard deviation from 
# all the flat files available in a given directctory structure.
#
###############################################################################
#
usage="Usage: `basename $0` -d dtg -W depend -q pbs_job [t|f]" 
while getopts ":d:W:q:n:" option
do
  case $option in
      d  ) dtg=$OPTARG;;
      q  ) pbs_job=$OPTARG;;
      W  ) depend=$OPTARG;;
      n  ) nnest=$OPTARG;;
      \? ) echo $usage ; exit 65;;
      *  ) echo $usage ; exit 65;;
  esac
done

PATH_CONFIG=CONFIG/FILE/GOES/HERE
. ${PATH_CONFIG}
: ${dtg:=$cdtg_beg}
: ${pbs_job:='t'}
: ${nnest:=1}

BATCH="qsub"
if [[ $depend ]]; then BATCH="${BATCH} -W depend=afterok:${depend}"; fi

scrDir=${COAMPS_RUN_BASE_DIR}
scriptDir=${scrDir}/scripts
LOG_DIR=${scrDir}/log
DART_BIN=${scrDir}/bin
create_incr=${scriptDir}/create_increment.${dtg}.sh

dtgm1=`echo $dtg -$icycle | ${DART_BIN}/advance_time`

. ${scriptDir}/HPC_CONFIG.sh
MACH_LAYOUT
let "NODES = (POST_PROCS / CPU_PER_NODE) + (POST_PROCS % CPU_PER_NODE)"
MPI_CMD=${POST_MPI_CMD}

if [ -e ${create_incr} ]; then rm ${create_incr}; fi
if [ ${pbs_job} == 't' ]; then
cat > ${create_incr} << EOF
#!/bin/bash
`PBS_DIRECTIVE COMP_INCR ${dtg} ${POST_PROCS} ${NODES} ${WALLTIME_POST} ${LOG_DIR} n n` 
EOF
else
BATCH=''
cat > ${create_incr} << EOF
#!/bin/bash
EOF
fi

cat >> ${create_incr} << EOF
#
MPI_CMD='${MPI_CMD}'
dtg=${dtg}
dtgm1=${dtgm1}
icycle=${icycle}
scrDir=${scrDir}
WORKDIR=\${scrDir}/${dtg}
COAMPS_DATA=${COAMPS_DATA}
#
scriptDir=\${scrDir}/scripts
DART_BIN=\${scrDir}/bin
NL_STRIP=\${scriptDir}/strip_namelist.pl

incrDir=\${scrDir}/incr
incrDat=\${incrDir}/data

if [ ! -e \${incrDir} ]; then mkdir -p \${incrDir}; fi
if [ ! -e \${incrDat} ]; then mkdir -p \${incrDat}; fi
if [ ! -e \${WORKDIR} ]; then mkdir -p \${WORKDIR}; fi
if [ ! -e \${WORKDIR}/input.nml ]; then cp \${scrDir}/input.nml \${WORKDIR}; fi

icycle_str=\`printf %04d 0\`

ensdir=\`printf \${COAMPS_DATA} 1\`
#ensdir_analyses=\${ensdir}/analyses
ensdir_analyses=\${ensdir}

cd \${ensdir_analyses}
#flist=( \`ls *\${dtg}*_analfld\` )
flist=( \`ls *sig*\${dtg}_\${icycle_str}*_fcstfld\` )
nfiles=\${#flist[@]}

cd \${ensdir}
datahd_dtgm1=\`ls datahd*_\${dtgm1}_*\`
terrht_dtgm1=\`ls terrht*_\${dtgm1}_*\`

meta_files=\`echo \${datahd_dtgm1[@]} \${terrht_dtgm1[@]}\`
for file in \${meta_files[@]}; do
  target_file=\`echo \$file | sed "s/\${dtgm1}/\${dtg}/"\`
  cp \${file} \${incrDat}/\${target_file}
done

cd \${WORKDIR}
var_file='./compute_incr.vars'
if [ -e \${var_file} ]; then rm \${var_file}; fi
cat > \${var_file} << EOF_NL
kka
n
m
EOF_NL
\${NL_STRIP} ${scrDir}/namelist \${var_file} compute_incr
rm -f \${var_file}

cat >> compute_incr.nml << EOF_NL
&field_proc
  cdtg   = '\${dtg}', 
  cdtgm1 = '\${dtgm1}', 
  icycle = \${icycle}, 
  ens_size = ${ENSEMBLE_SIZE}, 
  nfiles = \${nfiles}, 
  dsnrff = '\${scrDir}/', 
  flist  = 
EOF_NL
for coamps_file in \${flist[@]}; do
cat >> compute_incr.nml << EOF_NL
           '\${coamps_file}',
EOF_NL
done
cat >> compute_incr.nml << EOF_NL
/
EOF_NL

\${MPI_CMD} \${DART_BIN}/create_increment 
rm -f compute_incr.nml

for inest in \`seq 1 $nnest\`; do
  ${HOME}/diagnostic_cdf/diag_cdf.sh \${dtg} \${inest} -name=incr_\${dtg} -indir=\${incrDir}/data -odir=\${incrDir} -begt=0 -endt=0
done

exit 0
EOF

chmod 777 ${create_incr}

${BATCH} ${create_incr}
rm -f ${create_incr}

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

