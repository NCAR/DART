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
# Updates the coamps and the dart namelist files for the new dtg.
#
###############################################################################
#
# Argument should supply where we can find the path configuration file
usage="Usage: `basename $0` [-n init_pmo (t|f) -d dtg -i icycle]" 

while getopts ":d:n:i:" option
do
  case $option in
      d  ) dtg=$OPTARG;;
      n  ) init_pmo=$OPTARG;;
      i  ) icycle_loc=$OPTARG;;
      \? ) echo $usage ; exit 65;;
      *  ) echo $usage ; exit 65;;
  esac
done

PERL='perl -i -p -e'
SED='sed -i -e'
CP='cp -f'
LN='ln -sf'

PATH_CONFIG=CONFIG/FILE/GOES/HERE
. ${PATH_CONFIG}
: ${init_pmo:='f'}
: ${dtg:=$cdtg_beg}
: ${icycle_loc:=$icycle}
: ${NUM_MEMS=$ENSEMBLE_SIZE}

init_ens=`echo ${init_ens:0:1}`

scrDir=${COAMPS_RUN_BASE_DIR}
DART_BIN=${scrDir}/bin
scriptDir=${scrDir}/scripts
MEMBER_DIR=${scrDir}/perfect
MEMBER_DAT=${MEMBER_DIR}/data

COAMPS_PATH=${MEMBER_DAT}/
GLOBAL_PATH=`printf ${NOGAPS_DETERMIN} ${dtg}`

cd ${scrDir}

dtgm1=`echo $dtg   -$icycle_loc | ${DART_BIN}/advance_time`
dtgm2=`echo $dtgm1 -$icycle_loc | ${DART_BIN}/advance_time`
ds_time=(`echo $dtg 0 -g | ${DART_BIN}/advance_time`)

# Make sure that everything that needs to exist, exists
if [ ! -e ${MEMBER_DIR} ]; then mkdir -p ${MEMBER_DIR}; fi
if [ ! -e ${MEMBER_DAT} ]; then mkdir -p ${MEMBER_DAT}; fi
if [ ! -e ${MEMBER_DAT}/backward ]; then mkdir -p ${MEMBER_DAT}/backward; fi
if [ ! -e ${MEMBER_DAT}/forward ];  then mkdir -p ${MEMBER_DAT}/forward;  fi
if [ ! -e ${MEMBER_DAT}/analyses ]; then mkdir -p ${MEMBER_DAT}/analyses; fi
if [ ! -e ${MEMBER_DAT}/state.vars ];   then ${LN} ${scrDir}/state.vars   ${MEMBER_DAT}; fi
if [ ! -e ${MEMBER_DAT}/convert.vars ]; then ${LN} ${scrDir}/convert.vars ${MEMBER_DAT}; fi
if [ ! -e ${MEMBER_DAT}/input.nml ];    then ${LN} ${scrDir}/input.nml    ${MEMBER_DAT}; fi
if [ ! -e ${MEMBER_DIR}/namelist ];     then ${CP} ${scrDir}/namelist     ${MEMBER_DIR}/namelist; fi

echo "Updating COAMPS namelist for ensemble member ${MEMBER_DIR}"
echo "  dtg   - update to $dtg"
${PERL} "s:^(\s*cdtg\s*=\s*).\d+.(,).*:\1'${dtg}'\2:" ${MEMBER_DIR}/namelist
echo "  dsngff - update to ${GLOBAL_PATH}"
${PERL} "s:^(\s*dsngff\s*=\s*\').*(\',):\1${GLOBAL_PATH}\2:" ${MEMBER_DIR}/namelist
echo "  dsnrff - update to ${COAMPS_PATH}"
${PERL} "s:^(\s*dsnrff\s*=\s*\').*(\',):\1${COAMPS_PATH}\2:" ${MEMBER_DIR}/namelist

if [ $init_pmo = 't' ]; then
  echo "  first_obs_days   - update to ${ds_time[0]}"
  ${SED} "/&perfect_model_obs_nml/,/first_obs_days/ s/^\(\s*first_obs_days\s*\)=.*/\1=${ds_time[0]},/" ./input.nml

  echo "  first_obs_second - update to ${ds_time[1]}"
  ${SED} "/&perfect_model_obs_nml/,/first_obs_seconds/ s/^\(\s*first_obs_seconds\s*\)=.*/\1=${ds_time[1]},/" ./input.nml

  echo "  cdtg             - update to ${dtg}"
  ${SED} "/&model_nml/,/cdtg/ s/^\(\s*cdtg\s*\=\s*\).*/\1\'${dtg}\',/" ./input.nml
fi

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

