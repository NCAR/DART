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

function hms_replace
{
  SED='sed -i -e'
  hms=`printf %3d,%3d,%3d, ${1} 0 0`
  nnest_loc=${2}
  param_loc=${3}
  nlist=${4}
  hms_loc=''
  ${SED} "s/\(^.*${param_loc}\s*=\)\(\s*\).*/\1\2/" ${nlist}
  ${SED} "/${param_loc}/,/=/ s/^\s*[0-9],.*/####/" -e '/####/d' ${nlist}
  for n in `seq 1 ${nnest_loc}`; do
    hms_loc="${hms_loc} ${hms}"
    #${SED} "/${param_loc}\s*=\(\s*\).*/a\ ${hms}" ${nlist}
  done
  ${SED} "s/\(${param_loc}\s*=\s*\).*/\1${hms_loc}/" ${nlist}
}

# Argument should supply where we can find the path configuration file
usage="Usage: `basename $0` -n init_ens (t|f) -d dtg -i icycle -f is_fcst (t|f)"

while getopts ":d:n:i:f:m:" option
do
  case $option in
      d  ) dtg=$OPTARG;;
      n  ) init_ens=$OPTARG;;
      m  ) mems_to_run=$OPTARG;;
      i  ) icycle_loc=$OPTARG;;
      f  ) is_fcst=$OPTARG;;
      \? ) echo $usage ; exit 65;;
      *  ) echo $usage ; exit 65;;
  esac
done

PERL='perl -i -p -e'
SED='sed -i -e'
CP='cp -f'
MK='mkdir -p'
LN='ln -sf'
digits=5

PATH_CONFIG=CONFIG/FILE/GOES/HERE
. ${PATH_CONFIG}
: ${init_ens:='f'}
: ${dtg:=$cdtg_beg}
: ${icycle_loc:=$icycle}
: ${NUM_MEMS=$ENSEMBLE_SIZE}
: ${is_first_assim:='f'}
: ${IS_REAL_DATA:='F'}
: ${is_fcst:='f'}
: ${fcstlen:=$icycle}
: ${mems_to_run:=`seq 1 $NUM_MEMS`}
: ${OCNOBS_PATH:='${BIGDISK}/COAMPS_DATA/ocnqc/'}
: ${is_fcp_bndy:='false'}

init_ens=${init_ens:0:1}

mems_to_run=`echo ${mems_to_run} | sed 's/\,/ /g'`

scrDir=${COAMPS_RUN_BASE_DIR}
WORKDIR=${scrDir}/${dtg}
DART_BIN=${scrDir}/bin
scriptDir=${scrDir}/scripts

cd ${scrDir}

cdtg_beg_p1=`echo $cdtg_beg  ${fcstlen_init} | ${DART_BIN}/advance_time`
dtgm1=`echo $dtg   -$icycle_loc | ${DART_BIN}/advance_time`
dtgm2=`echo $dtgm1 -$icycle_loc | ${DART_BIN}/advance_time`
ds_time=(`echo $dtg 0 -g | ${DART_BIN}/advance_time`)

# Is this the first assimilation? 
if [ ${dtg} -eq ${cdtg_beg_p1} ]; then is_first_assim='t'; fi

if [ ${dtg:8:2} = '06' -o ${dtg:8:2} = '18' ]; then
  dtg_ngt=`echo $dtg   -6 | ${DART_BIN}/advance_time`
else
  dtg_ngt=${dtg}
fi

NGT_FILE=ngt${dtg_ngt}
NGT_DIR=`printf ${OBS_PATH} ${dtg_ngt}`
if [ -e ${NGT_DIR}/${NGT_FILE} ]; then ${CP} ${NGT_DIR}/${NGT_FILE} ${WORKDIR} ; fi

for dir in ${mems_to_run}; do
  # Generate the name of the directory currently holding COAMPS data.  
  # Information will be copied *FROM* this directory into MEMBER_DIR.
  # 
  # The printf statement allows the substitution of the ensemble member
  # number into MEMBER_DIR if needed, otherwise it leaves it as is.
  # It goes without saying that you should not use any character
  # string that printf will recognize in your directory names.
  
  # Recast the ensemble member number to a directory name - just a 
  # fixed number of digits is the default
  MEMBER_DIR=`printf "%0${digits}d" ${dir}`
  if [ ${is_fcst:0:1} == 't' -o ${is_fcst:0:1} == 'T' ]; then
    MEMBER_RUN=${MEMBER_DIR}/fcst.${dtg}
  else
    MEMBER_RUN=${MEMBER_DIR}
  fi
  MEMBER_DAT=${MEMBER_DIR}/data
  if [ $is_fcp_bndy = true ]; then
    GLOBAL_PATH=`printf ${NOGAPS_DETERMIN} ${dtg}`
  else
    GLOBAL_PATH=`printf ${NOGAPS_ENSEMBLE} ${dtg} ${dir}`
  fi
  COAMPS_PATH=`printf ${COAMPS_DATA}/ ${dir}`

  if [ ! -e ${MEMBER_RUN} ]; then 
    echo "  Creating directory for ensemble member ${MEMBER_RUN}..."
    mkdir -p ${MEMBER_RUN}
  fi

  ${CP} ${scrDir}/namelist ${MEMBER_RUN}
  if [   -e ${WORKDIR}/${NGT_FILE} ];     then ${CP} ${WORKDIR}/${NGT_FILE} ${MEMBER_RUN}; fi
  if [ ! -e ${MEMBER_DAT} ];              then ${MK} ${MEMBER_DAT}; fi
  if [ ! -e ${MEMBER_DAT}/backward ];     then ${MK} ${MEMBER_DAT}/backward; fi
  if [ ! -e ${MEMBER_DAT}/forward ];      then ${MK} ${MEMBER_DAT}/forward; fi
  if [ ! -e ${MEMBER_DAT}/analyses ];     then ${MK} ${MEMBER_DAT}/analyses; fi
  if [ ! -e ${MEMBER_DAT}/state.vars ];   then ${LN} ${scrDir}/state.vars   ${MEMBER_DAT}; fi
  if [ ! -e ${MEMBER_DAT}/convert.vars ]; then ${LN} ${scrDir}/convert.vars ${MEMBER_DAT}; fi
  if [ ! -e ${MEMBER_DAT}/perturb.vars -a ${is_fcp_bndy} == 'true' ]; then ${LN} ${scrDir}/perturb.vars ${MEMBER_DAT}; fi
  if [ ! -e ${MEMBER_DAT}/input.nml ];    then ${LN} ${scrDir}/input.nml    ${MEMBER_DAT}; fi
 
  # Copy our template namelist into each ensemble member's directory
  echo "Updating COAMPS namelist for ensemble member ${MEMBER_RUN}"
  echo "  dtg   - update to $dtg"
  ${PERL} "s:^(\s*cdtg\s*=\s*).\d+.(,).*:\1'${dtg}'\2:" ${MEMBER_RUN}/namelist
  echo "  dsngff - update to ${GLOBAL_PATH}"
  ${PERL} "s:^(\s*dsngff\s*=\s*\').*(\',):\1${GLOBAL_PATH}\2:" ${MEMBER_RUN}/namelist
  echo "  dsnrff - update to ${COAMPS_PATH}"
  ${PERL} "s:^(\s*dsnrff\s*=\s*\').*(\',):\1${COAMPS_PATH}\2:" ${MEMBER_RUN}/namelist

  if [ ${is_fcst:0:1} == 't' -o ${is_fcst:0:1} == 'T' ]; then
    echo "  icycle - update to ${icycle}"
    ${SED} "s/^\(\s*icycle\s*=\s*\).*,\s*/\1${icycle_loc},/" ${MEMBER_RUN}/namelist
  else
    echo "  icycle - update to ${icycle_loc}"
    ${SED} "s/^\(\s*icycle\s*=\s*\).*,\s*/\1${icycle_loc},/" ${MEMBER_RUN}/namelist
  fi

  nnest=`grep nnest ${MEMBER_RUN}/namelist | sed 's/,//g' | awk -F= '{print $2}'`
  if [ ${is_fcst:0:1} == 't' -o ${is_fcst:0:1} == 'T' ]; then
    echo "  ktaust - update to ${icycle_loc},  0,  0"
    hms_replace ${icycle_loc} 1 ktaust ${MEMBER_RUN}/namelist 
    echo "  ktauf  - update to ${fcstlen},  0,  0"
    hms_replace ${fcstlen} ${nnest} ktauf ${MEMBER_RUN}/namelist 
  else
    echo "  ktaust - update to 0,  0,  0"
    hms_replace 0 1 ktaust ${MEMBER_RUN}/namelist 
    echo "  ktauf  - update to ${icycle_loc},  0,  0"
    hms_replace ${icycle_loc} ${nnest} ktauf ${MEMBER_RUN}/namelist 
  fi

if [ ${is_fcst:0:1} != 't' -a ${is_fcst:0:1} != 'T' ]; then
# Write new ncoda namelists
cat << EOF > ${MEMBER_RUN}/odsetnl 
 &dsetnl
 dsoudat = '${OCNOBS_PATH}',
 dsngff  = '${GLOBAL_PATH}',
 dsnrff  = '${COAMPS_PATH}',
 dsorff  = '${COAMPS_PATH}',
 dsowrk  = '${COAMPS_PATH}',
 dsoclim = '${COAMPS_DATABASE}/codaclim/',
 dsomdas = '${COAMPS_DATABASE}/codaclim/',
 dsogdem = '${COAMPS_DATABASE}/gdem/',
 /
EOF

cat << EOF > ${MEMBER_RUN}/oanl 
 &oanl
 locn3d  = .false.,
 rscl    = 3., 4., 4., 4.,
 /
EOF
fi

done

if [ ${is_fcst:0:1} == 't' -o ${is_fcst:0:1} == 'T' ]; then exit 0; fi

# Update the dart namelist
echo "Updating DART namelist"

echo "  first_obs_days   - update to ${ds_time[0]}"
${SED} "/&filter_nml/,/first_obs_days/ s/^\(\s*first_obs_days\s*\)=.*/\1=${ds_time[0]},/" ./input.nml

echo "  first_obs_second - update to ${ds_time[1]}"
${SED} "/&filter_nml/,/first_obs_seconds/ s/^\(\s*first_obs_seconds\s*\)=.*/\1=${ds_time[1]},/" ./input.nml

echo "  last_obs_days    - update to ${ds_time[0]}"
${SED} "/&filter_nml/,/last_obs_days/ s/^\(\s*last_obs_days\s*\)=.*/\1= ${ds_time[0]},/" ./input.nml

echo "  last_obs_second  - update to ${ds_time[1]}"
${SED} "/&filter_nml/,/last_obs_seconds/ s/^\(\s*last_obs_seconds\s*\)=.*/\1= ${ds_time[1]},/" ./input.nml

if [ ${is_first_assim} == 't' ]; then
  echo "  cdtg             - update to $cdtg_beg"
  ${PERL} "s:^(\s*cdtg\s*=\s*).\d+.(,).*:\1'${cdtg_beg}'\2:" ./input.nml
else
  echo "  cdtg             - update to $dtgm1"
  ${PERL} "s:^(\s*cdtg\s*=\s*).\d+.(,).*:\1'${dtgm1}'\2:" ./input.nml
fi

if [ ! -e ${WORKDIR} ]; then mkdir -p ${WORKDIR}; fi
${CP} ./input.nml   ${WORKDIR}
${CP} ./state.vars  ${WORKDIR}
if [ ${IS_REAL_DATA} == 'T' ]; then
  if [ $nnest -eq 1 ]; then
    innov_file_name="./innovations/innov_1_${dtg}"
  else
    innov_file_name="./innovations/innov_a_${dtg}"
  fi
  echo "  innov_file_name  - update to ${innov_file_name}"
  ${SED} ":&navdas_innov_nml:,:innov_file_name: s:\(^\s*innov_file_name\s*=\).*:\1 \'${innov_file_name}\',:" ${WORKDIR}/input.nml
  ${SED} ":&navdas_innov_nml:,:ngt_file_name: s:\(^\s*ngt_file_name\s*=\).*:\1 \'${NGT_FILE}\',:" ${WORKDIR}/input.nml
else
  ${CP} ./obs_seq.out ${WORKDIR}
fi

if [ $init_ens == 'f' ]; then
  if [ ${is_first_assim} == 't' ]; then
    echo "  inf_initial_from_restart -  .false."
    ${SED} "s/^\(\s*inf_initial_from_restart\s*=\s*\).*/\1.false., .false.,/" ${WORKDIR}/input.nml
    echo "  inf_sd_from_restart      - .false."
    ${SED} "s/^\(\s*inf_sd_initial_from_restart\s*=\s*\).*/\1.false., .false.,/" ${WORKDIR}/input.nml
  else
    cd ${WORKDIR}
    inf_initial_from_restart=( `grep inf_initial_from_restart ./input.nml | sed 's/,//g' | awk -F= '{print $2}'` )
    inf_in_file_name=( `grep inf_in_file_name ./input.nml | sed 's/,//g' | sed "s/'//g" | awk -F= '{print $2}'` )
    inf_out_file_name=( `grep inf_out_file_name ./input.nml | sed 's/,//g' | sed "s/'//g" | awk -F= '{print $2}'` )
    for i in 0 1; do
      if [ ${inf_initial_from_restart[$i]} == '.true.' -a -e ${scrDir}/${dtgm1}/${inf_out_file_name[$i]} ]; then
        ${CP} ${scrDir}/${dtgm1}/${inf_out_file_name[$i]} ${WORKDIR}/${inf_in_file_name[$i]} 
      fi
    done
  fi
fi

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

