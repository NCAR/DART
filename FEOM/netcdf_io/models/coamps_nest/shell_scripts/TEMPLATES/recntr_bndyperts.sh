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
usage="Usage: `basename $0` -c path -d dtg -n dtg_new -b bndyDir" 
while getopts ":d:r:c:b:" option
do
  case $option in
      c  ) PATH_CONFIG=$OPTARG;;
      d  ) dtg=$OPTARG;;
      r  ) dtg_new=$OPTARG;;
      b  ) bndyDir=$OPTARG;;
      \? ) echo $usage ; exit 65;;
      *  ) echo $usage ; exit 65;;
  esac
done

PATH_CONFIG=CONFIG/FILE/GOES/HERE
. ${PATH_CONFIG}
: ${dtg:=$cdtg_beg}
: ${dtg_new:=$dtg}

BATCH="qsub"

scrDir=${COAMPS_RUN_BASE_DIR}
scriptDir=${scrDir}/scripts
LOG_DIR=${scrDir}/log
scale_perts=${scriptDir}/recntr_bndyperts.${dtg}.sh

: ${bndyDir=$scrDir}

. ${scriptDir}/HPC_CONFIG.sh
MACH_LAYOUT
let "NODES = (POST_PROCS / CPU_PER_NODE) + (POST_PROCS % CPU_PER_NODE)"
#MPI_CMD=${POST_MPI_CMD}
MPI_CMD=

if [ -e ${scale_perts} ]; then rm ${scale_perts}; fi
cat > ${scale_perts} << EOF
#!/bin/bash
`PBS_DIRECTIVE SCALE_NGP_PERTS ${dtg} ${POST_PROCS} ${NODES} ${WALLTIME_POST} ${LOG_DIR} n n` 
EOF

cat >> ${scale_perts} << EOF
#
MPI_CMD='${MPI_CMD}'
dtg=${dtg}
dtg_new=${dtg_new}
scrDir=${scrDir}
bndyDir=${bndyDir}
WORKDIR=\${scrDir}/${dtg}
COAMPS_DATA=${COAMPS_DATA}
#
scriptDir=\${scrDir}/scripts
DART_BIN=\${scrDir}/bin

if [ ! -e \${WORKDIR} ]; then mkdir -p \${WORKDIR}; fi
if [ ! -e \${WORKDIR}/input.nml ]; then cp \${scrDir}/input.nml \${WORKDIR}; fi

ensdir=\`printf \${bndyDir}/data 1\`
cd \${ensdir}
flist=( \`ls *\${dtg}*bndyfld\` )
nfiles=\${#flist[@]}

cd \${WORKDIR}
if [ -e \${recntr_bndyperts.nml} ]; then rm \${recntr_bndyperts.nml}; fi
cat > recntr_bndyperts.nml << EOF_NL
&field_proc
  ens_size = ${ENSEMBLE_SIZE}, 
  nfiles = \${nfiles}, 
  dsnrff = '\${bndyDir}/', 
  dsnrff_out = '\${bndyDir}/', 
  dsnrff_new = '\${bndyDir}/dtrmn/', 
  dtg_new = '\${dtg_new}/', 
  flist  = 
EOF_NL
for coamps_file in \${flist[@]}; do
cat >> recntr_bndyperts.nml << EOF_NL
           '\${coamps_file}',
EOF_NL
done
cat >> recntr_bndyperts.nml << EOF_NL
/
EOF_NL

\${MPI_CMD} \${DART_BIN}/recntr_bndyperts 
rm -f recntr_bndyperts.nml

exit 0
EOF

chmod 777 ${scale_perts}

#${BATCH} ${scale_perts}
${scale_perts}
rm -f ${scale_perts}

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

