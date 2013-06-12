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
usage="Usage: `basename $0` -d dtg -W depend -f [t|f]" 
while getopts ":d:W:f:n:" option
do
  case $option in
      d  ) dtg=$OPTARG;;
      W  ) depend=$OPTARG;;
      f  ) is_fcst=$OPTARG;;
      n  ) nnest=$OPTARG;;
      \? ) echo $usage ; exit 65;;
      *  ) echo $usage ; exit 65;;
  esac
done

PATH_CONFIG=CONFIG/FILE/GOES/HERE
. ${PATH_CONFIG}
: ${dtg:=$cdtg_beg}
: ${is_fcst:='f'}
: ${nnest:=1}

BATCH="qsub"
if [[ $depend ]]; then BATCH="${BATCH} -W depend=afterok:${depend}"; fi
#if [[ $depend ]]; then depend="-W depend=afterok:${depend}"; fi

scrDir=${COAMPS_RUN_BASE_DIR}
scriptDir=${scrDir}/scripts
LOG_DIR=${scrDir}/log
create_mean=${scriptDir}/create_mean_std.${dtg}.sh

. ${scriptDir}/HPC_CONFIG.sh
MACH_LAYOUT
let "NODES = (POST_PROCS / CPU_PER_NODE) + (POST_PROCS % CPU_PER_NODE)"
MPI_CMD=${POST_MPI_CMD}

if [ -e ${create_mean} ]; then rm ${create_mean}; fi
cat > ${create_mean} << EOF
#!/bin/bash
`PBS_DIRECTIVE COMP_MEAN ${dtg} ${POST_PROCS} ${NODES} ${WALLTIME_POST} ${LOG_DIR} n n` 
EOF

cat >> ${create_mean} << EOF
#
MPI_CMD='${MPI_CMD}'
dtg=${dtg}
icycle=${icycle}
scrDir=${scrDir}
WORKDIR=\${scrDir}/${dtg}
COAMPS_DATA=${COAMPS_DATA}
is_fcst=${is_fcst:0:1}
fcstlen=${fcstlen}
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
flist=( \`ls *1a*\${dtg}*fcstfld\` )
for inest in \`seq 2 ${nnest}\`; do
  flist=( \${flist[@]}  \`ls *\${inest}a*\${dtg}*fcstfld\` )
done

nfiles=\${#flist[@]}

cp datahd*_\${dtg}_* \${sdevDat}
cp datahd*_\${dtg}_* \${meanDat}
cp terrht*_\${dtg}_* \${sdevDat}
cp terrht*_\${dtg}_* \${meanDat}

cd \${WORKDIR}
var_file='./compute_mean.vars'
if [ -e \${var_file} ]; then rm \${var_file}; fi
cat > \${var_file} << EOF_NL
kka
n
m
EOF_NL
\${NL_STRIP} ${scrDir}/namelist \${var_file} compute_mean
rm -f \${var_file}

cat >> compute_mean.nml << EOF_NL
&field_proc
  cdtg   = '\${dtg}', 
  ens_size = ${ENSEMBLE_SIZE}, 
  nfiles = \${nfiles}, 
  dsnrff = '\${scrDir}/', 
  flist  = 
EOF_NL
for coamps_file in \${flist[@]}; do
cat >> compute_mean.nml << EOF_NL
           '\${coamps_file}',
EOF_NL
done
cat >> compute_mean.nml << EOF_NL
/
EOF_NL

\${MPI_CMD} \${DART_BIN}/create_mean_std 
rm -f compute_mean.nml

for inest in \`seq 1 $nnest\`; do
  if [ \${is_fcst} == 't' -o \${is_fcst} == 'T' ]; then 
    ${HOME}/diagnostic_cdf/diag_cdf.sh \${dtg} \${inest} -name=mean_\${dtg} -indir=\${meanDir}/data -odir=\${meanDir} -begt=0 -endt=\${fcstlen} -intt=\${icycle}
    ${HOME}/diagnostic_cdf/diag_cdf.sh \${dtg} \${inest} -name=sdev_\${dtg} -indir=\${sdevDir}/data -odir=\${sdevDir} -begt=0 -endt=\${fcstlen} -intt=\${icycle}
  else
    ${HOME}/diagnostic_cdf/diag_cdf.sh \${dtg} \${inest} -name=mean_\${dtg} -indir=\${meanDir}/data -odir=\${meanDir} -begt=0 -endt=\${icycle} -intt=\${icycle}
    ${HOME}/diagnostic_cdf/diag_cdf.sh \${dtg} \${inest} -name=sdev_\${dtg} -indir=\${sdevDir}/data -odir=\${sdevDir} -begt=0 -endt=\${icycle} -intt=\${icycle}
  fi
done

exit 0
EOF

chmod 777 ${create_mean}

echo ${BATCH} ${create_mean}
${BATCH} ${create_mean}
rm -f ${create_mean}

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

