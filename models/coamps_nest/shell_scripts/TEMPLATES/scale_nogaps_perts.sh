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
# Runs the program that will calculate the mean and standard deviation from 
# all the flat files available in a given directctory structure.
#
###############################################################################
#
usage="Usage: `basename $0` -c path -p pert_scale -d dtg_beg -e dtg_end -W depend" 
while getopts ":d:e:c:p:W:" option
do
  case $option in
      p  ) pert_scale=$OPTARG;;
      c  ) PATH_CONFIG=$OPTARG;;
      d  ) dtg_beg=$OPTARG;;
      e  ) dtg_end=$OPTARG;;
      W  ) depend=$OPTARG;;
      \? ) echo $usage ; exit 65;;
      *  ) echo $usage ; exit 65;;
  esac
done

PATH_CONFIG=CONFIG/FILE/GOES/HERE
. ${PATH_CONFIG}
: ${dtg:=$cdtg_beg}
: ${pert_scale:=0.5}

BATCH="qsub"
if [[ $depend ]]; then BATCH="${BATCH} -W depend=afterok:${depend}"; fi

scrDir=${COAMPS_RUN_BASE_DIR}
scriptDir=${scrDir}/scripts
LOG_DIR=${scrDir}/log
scale_perts=${scriptDir}/scale_nogaps_perts.${dtg}.sh

. ${scriptDir}/HPC_CONFIG.sh
MACH_LAYOUT
let "NODES = (POST_PROCS / CPU_PER_NODE) + (POST_PROCS % CPU_PER_NODE)"
MPI_CMD=${POST_MPI_CMD}

if [ -e ${scale_perts} ]; then rm ${scale_perts}; fi
cat > ${scale_perts} << EOF
#!/bin/bash
`PBS_DIRECTIVE SCALE_NGP_PERTS ${dtg} ${POST_PROCS} ${NODES} ${WALLTIME_POST} ${LOG_DIR} n n` 
EOF

cat >> ${scale_perts} << EOF
#
MPI_CMD='${MPI_CMD}'
dtg=${dtg}
scrDir=${scrDir}
WORKDIR=\${scrDir}/${dtg}
COAMPS_DATA=${COAMPS_DATA}
#
scriptDir=\${scrDir}/scripts
DART_BIN=\${scrDir}/bin
NL_STRIP=\${scriptDir}/strip_namelist.pl

meanDir=\${scrDir}/mean
meanDat=\${meanDir}/data

sdevDir=\${scrDir}/sdev
sdevDat=\${sdevDir}/data

if [ ! -e \${meanDir} ]; then mkdir -p \${meanDir}; fi
if [ ! -e \${meanDat} ]; then mkdir -p \${meanDat}; fi
if [ ! -e \${sdevDir} ]; then mkdir -p \${sdevDir}; fi
if [ ! -e \${sdevDat} ]; then mkdir -p \${sdevDat}; fi
if [ ! -e \${WORKDIR} ]; then mkdir -p \${WORKDIR}; fi
if [ ! -e \${WORKDIR}/input.nml ]; then cp \${scrDir}/input.nml \${WORKDIR}; fi

ensdir=\`printf \${COAMPS_DATA} 1\`
cd \${ensdir}
flist=( \`ls *\${dtg}*fcstfld\` )
#flist=( \`ls *sig*\${dtg}*fcstfld\` )
nfiles=\${#flist[@]}

cp datahd*_\${dtg}_* \${sdevDat}
cp datahd*_\${dtg}_* \${meanDat}
cp terrht*_\${dtg}_* \${sdevDat}
cp terrht*_\${dtg}_* \${meanDat}

cd \${WORKDIR}
var_file='./scale_nogaps_perts.vars'
if [ -e \${var_file} ]; then rm \${var_file}; fi
cat > \${var_file} << EOF_NL
kka
n
m
EOF_NL
\${NL_STRIP} ${scrDir}/namelist \${var_file} scale_nogaps_perts
rm -f \${var_file}

cat >> scale_nogaps_perts.nml << EOF_NL
&field_proc
  cdtg   = '\${dtg}', 
  ens_size = ${ENSEMBLE_SIZE}, 
  nfiles = \${nfiles}, 
  dsnrff = '\${scrDir}/', 
  dsnrff_out = '\${scrDir}/', 
  pert_scale = ${pert_scale}
  flist  = 
EOF_NL
for coamps_file in \${flist[@]}; do
cat >> scale_nogaps_perts.nml << EOF_NL
           '\${coamps_file}',
EOF_NL
done
cat >> scale_nogaps_perts.nml << EOF_NL
/
EOF_NL

\${MPI_CMD} \${DART_BIN}/scale_nogaps_perts 
rm -f scale_nogaps_perts.nml

exit 0
EOF

chmod 777 ${scale_perts}

${BATCH} ${scale_perts}
rm -f ${scale_perts}

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

