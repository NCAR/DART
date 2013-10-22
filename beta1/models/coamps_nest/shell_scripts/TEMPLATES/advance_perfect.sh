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
usage="Usage: `basename $0` [-d dtg -i icycle]" 
while getopts ":d:i:" option
do
  case $option in
      d  ) dtg=$OPTARG;;
      i  ) icycle_loc=$OPTARG;;
      \? ) echo $usage ; exit 65;;
      *  ) echo $usage ; exit 65;;
  esac
done

echo "#"
echo "#################### BEGIN advance_perfect.sh ####################"
echo "#"

BATCH="qsub"
SED='sed -i -e'

PATH_CONFIG=CONFIG/FILE/GOES/HERE
. ${PATH_CONFIG}
: ${OSSE_FREE_RUN:='F'}

scrDir=${COAMPS_RUN_BASE_DIR}
ensdir=perfect
COAMPS_ENS=${scrDir}/${ensdir}
DART_BIN=${scrDir}/bin
scriptDir=${scrDir}/scripts

. ${scriptDir}/HPC_CONFIG.sh
MACH_LAYOUT
let "COAMM_NODES = (COAMM_PROCS / CPU_PER_NODE) + (COAMM_PROCS % CPU_PER_NODE)"

dtg_nml=`grep cdtg ${COAMPS_ENS}/namelist | awk -F"'" '{print $2}'`
: ${dtg:=$dtg_nml}
: ${icycle_loc:=$icycle}
: ${init_pmo:='f'}
cd ${scrDir}
dtgp1=`echo ${dtg} ${icycle_loc} | ${DART_BIN}/advance_time`
if [ ${OSSE_FREE_RUN} == 'T' ]; then PMO_FCST_HR=${icycle}; fi

tau_offset=`expr ${fcstlen_init} - ${PMO_FCST_HR}`
: ${dtg_beg:=`echo ${cdtg_beg} ${tau_offset} | ${DART_BIN}/advance_time`}
if [ $dtg -eq $dtg_beg ]; then init_pmo='t'; fi

# Compute number of hours offset from intial time
# store value in tau_offset
t1=( `echo ${cdtg_beg} 0 -g  | ${DART_BIN}/advance_time` )
t2=( `echo ${dtg}      0 -g  | ${DART_BIN}/advance_time` )
dt_day=`expr ${t2[0]} - ${t1[0]}`
dt_hr=`expr ${t2[1]} - ${t1[1]}`
let 'dt_day*=24'
let 'dt_hr/=3600'
tau_offset=`expr ${dt_hr} + ${dt_day}`

advance_perfect=${scriptDir}/advance_perfect.${dtg}.sh
if [ -e ${advance_perfect} ]; then rm ${advance_perfect}; fi
cat > ${advance_perfect} << EOF
#!/bin/bash
#
SED='${SED}'

scrDir=${scrDir}
scriptDir=\${scrDir}/scripts
DART_BIN=\${scrDir}/bin

COAMPS_BIN_DIR=${COAMPS_HOME}/bin
COAMPS_DART_DIR=${DART_HOME}
NL_STRIP=\${scriptDir}/strip_namelist.pl

echo Entering ${advance_perfect}
echo   Running on host \`hostname\`
echo   Time is \`date\`
echo   Directory is \`pwd\`
echo   Working directory is: \$scrDir

MPI_FCST='${COAMM_MPI_CMD}'
MPI_BNDY='${COAMA_MPI_CMD}'

element=perfect
digits=5
ensdir=\$element
COAMPS_ENS=\${scrDir}/\${ensdir}
logtime=$dtg

if [ $PMO_EXISTS == 'F' ]; then
# NOT WORKING RIGHT NOW 
  cd \${COAMPS_ENS}

  # Get the boundary conditions
  echo "Getting boundary conditions for \$element"
  \${MPI_BNDY} \${COAMPS_BIN_DIR}/atmos_analysis.exe > log.a.\${logtime}
  mpi_status=\$?
  if  [ \$mpi_status -ne 0 ]; then
    echo "Unable to get boundaries \$element!"
    exit 1
  fi

  # Run the COAMPS forecast
  echo "Running COAMPS ensemble member \$element"
  \${MPI_FCST} \${COAMPS_BIN_DIR}/atmos_forecast.exe > log.m.\${logtime}
  mpi_status=\$?
  if [ \$mpi_status -ne 0 ]; then
    echo "Unsuccessful completion of member \$element!"
    exit 1
  fi
else
 if [ ${OSSE_FREE_RUN} == 'T' ]; then
   KTAUF=${tau_offset}
 else
   KTAUF=${PMO_FCST_HR}
 fi
fi

# Now that the we've run the model, move into the COAMPS data directory 
# and move things back to the DART directory. 
cd \${COAMPS_ENS}/data
rm -f dart_vector

# Convert the resulting COAMPS file to a dart file
\${NL_STRIP} ../namelist convert.vars convert
\${SED} 's/\(^&convert.*\)/\1\n  is_pmo = \.true\.,/' ./convert.nml
\${SED} 's/\(^&convert.*\)/\1\n  is_first = \.false\.,/' ./convert.nml
\${SED} "s/^\(\s*ktauf\s*=\s*\).*,\s*/\1 \${KTAUF},  0,  0,/" ./convert.nml

if [ ${OSSE_FREE_RUN} == 'T' ]; then
  \${SED} "s/^\(\s*cdtg\s*=\s*\).*,\s*/\1\'${cdtg_beg}\',/" ./convert.nml
fi

\${DART_BIN}/trans_coamps_to_dart
if [ \$? -ne 0 ]; then
    echo "Unsuccessful translation to DART of member \$element!"
    exit 1
fi

\${scriptDir}/update_perfect.sh -d ${dtgp1}

if [ ${init_pmo} == 't' ]; then
  ud_file=perfect_ics
else
  ud_file=assim_model_state_ud.0001
fi
mv dart_vector \${scrDir}/\${ud_file}

rm -f \${scrDir}/filter_control*

exit 0
EOF

cd ${scrDir}
chmod 777 ${advance_perfect}
${advance_perfect}
rm -f ${advance_perfect}

echo "#################### END advance_perfect.sh ####################"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

