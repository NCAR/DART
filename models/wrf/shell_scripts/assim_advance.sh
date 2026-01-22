#!/usr/bin/env bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# datea, emember, paramfile are command-line arguments - OR -
# are set by a string editor (sed) command.

set -uo pipefail

datea="${1:-}"
emember="${2:-}"
paramfile="${3:-}"

if [[ -z "${datea}" || -z "${emember}" || -z "${paramfile}" ]]; then
  echo "usage: $0 <datea> <emember> <paramfile>"
  exit 2
fi

source "${paramfile}"

domains="${NUM_DOMAINS}"

start_time="$(date +%s)"
echo "host is $(hostname)"
echo "assim_advance.sh is running in $(pwd)"

cd "${RUN_DIR}"

read -r -a gdate  < <(echo "${datea} 0 -g" | "${RUN_DIR}/advance_time")

if (( ASSIM_INT_MINUTES <= 0 )); then
  read -r -a gdatef < <(echo "${datea} ${ASSIM_INT_HOURS} -g" | "${RUN_DIR}/advance_time")
else
  read -r -a gdatef < <(echo "${datea} ${ASSIM_INT_MINUTES}m -g" | "${RUN_DIR}/advance_time")
fi

yyyy="${datea:0:4}"
mm="${datea:4:2}"
dd="${datea:6:2}"
hh="${datea:8:2}"
nn="00"
ss="00"

#  Copy files to appropriate location
echo "${start_time}" > "${RUN_DIR}/start_member_${emember}"

# Go into member directory and generate the needed wrf.info file
cd "${RUN_DIR}/advance_temp${emember}"

icnum="$(printf "%04d" "${emember}")"
if [[ -e "${RUN_DIR}/advance_temp${emember}/wrf.info" ]]; then
  ${REMOVE} "${RUN_DIR}/advance_temp${emember}/wrf.info"
fi
: > wrf.info

if [[ "${SUPER_PLATFORM}" == "LSF queuing system" ]]; then

  cat > "${RUN_DIR}/advance_temp${emember}/wrf.info" <<EOF
${gdatef[1]}  ${gdatef[0]}
${gdate[1]}   ${gdate[0]}
$yyyy $mm $dd $hh $nn $ss
          $domains
mpirun.lsf ./wrf.exe
EOF

elif [[ "${SUPER_PLATFORM}" == "derecho" ]]; then

  # module load openmpi
  cat > "${RUN_DIR}/advance_temp${emember}/wrf.info" <<EOF
${gdatef[1]}  ${gdatef[0]}
${gdate[1]}   ${gdate[0]}
$yyyy $mm $dd $hh $nn $ss
           $domains
 mpiexec -n 4 -ppn 4 ./wrf.exe
EOF

fi

cd "${RUN_DIR}"

# filter_control accounts for multiple domains
# Appends input (filter_restart) and output (prior) in consecutive pairs
# Should be consistent with filter_control setup for first_advance.sh
# during intial forecast step

echo "${emember}" > "${RUN_DIR}/filter_control${icnum}"

dn=1
while (( dn <= domains )); do
  dchar="$(printf "%02d" "${dn}")"
  echo "filter_restart_d${dchar}.${icnum}" >> "${RUN_DIR}/filter_control${icnum}"
  echo "prior_d${dchar}.${icnum}"          >> "${RUN_DIR}/filter_control${icnum}"
  (( dn++ ))
done # loop through domains

#  integrate the model forward in time
"${RUN_DIR}/new_advance_model.sh" "${emember}" "${domains}" "filter_control${icnum}" "${paramfile}"
${REMOVE} "${RUN_DIR}/filter_control${icnum}"

end_time="$(date +%s)"
length_time=$(( end_time - start_time ))
echo "duration = ${length_time}"

exit 0

