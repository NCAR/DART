#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

set -uo pipefail

datea="$1"
emember="$2"
paramfile="$3"

source "$paramfile"

start_time="$(date +%s)"
echo "host is $(hostname)"

domains="$NUM_DOMAINS"

cd "${RUN_DIR}"

read -r -a gdate  < <(echo "$datea 0 -g"            | "${RUN_DIR}/advance_time")
read -r -a gdatef < <(echo "$datea $ASSIM_INT_HOURS -g" | "${RUN_DIR}/advance_time")
yyyy="${datea:0:4}"
mm="${datea:4:2}"
dd="${datea:6:2}"
hh="${datea:8:2}"
nn="00"
ss="00"

echo "$start_time" > "${RUN_DIR}/start_member_${emember}"

# Enter the member directory and generate the wrf.info file to prepare WRF integration
cd "${RUN_DIR}/advance_temp${emember}"

icnum="$(echo "$emember + 10000" | bc | cut -b2-5)"
if [[ -e "${RUN_DIR}/advance_temp${emember}/wrf.info" ]]; then
   ${REMOVE} "${RUN_DIR}/advance_temp${emember}/wrf.info"
fi

touch wrf.info

if [[ "${SUPER_PLATFORM}" == "LSF queuing system" ]]; then

   cat > "${RUN_DIR}/advance_temp${emember}/wrf.info" << EOF
 ${gdatef[1]}  ${gdatef[0]}
 ${gdate[1]}   ${gdate[0]}
 $yyyy $mm $dd $hh $nn $ss
           1
 mpirun.lsf ./wrf.exe
EOF

elif [[ "${SUPER_PLATFORM}" == "derecho" ]]; then

   export MPI_SHEPHERD=false

   cat > "${RUN_DIR}/advance_temp${emember}/wrf.info" << EOF
 ${gdatef[1]}  ${gdatef[0]}
 ${gdate[1]}   ${gdate[0]}
 $yyyy $mm $dd $hh $nn $ss
           $domains
 mpiexec -n 4 -ppn 4 ./wrf.exe
EOF

fi

cd "${RUN_DIR}"

# The filter_control file accounts for multiple domains and it appends
# input (filter_restart) and output (prior) in consecutive pairs for each domain..
# This is consistent with the filter_control setup within assim_advance.sh
# for assimilation steps after first_advance.

echo "$emember" > "${RUN_DIR}/filter_control${icnum}"

dn=1
while (( dn <= domains )); do
   dchar="$(echo "$dn + 100" | bc | cut -b2-3)"
   echo "filter_restart_d${dchar}.${icnum}" >> "${RUN_DIR}/filter_control${icnum}"
   echo "prior_d${dchar}.${icnum}"          >> "${RUN_DIR}/filter_control${icnum}"
   (( dn++ ))
done # loop through domains

# Call new_advance_model.sh to integrate WRF forward in time
"${RUN_DIR}/new_advance_model.sh" "${emember}" "${domains}" "filter_control${icnum}" "$paramfile"
${REMOVE} "${RUN_DIR}/filter_control${icnum}"

# Move the WRF forecast (prior) to the appropriate directory
mkdir -p "${OUTPUT_DIR}/${datea}/PRIORS"

dn=1
while (( dn <= domains )); do
   dchar="$(echo "$dn + 100" | bc | cut -b2-3)"	
   mv "${RUN_DIR}/prior_d${dchar}.${icnum}" "${OUTPUT_DIR}/${datea}/PRIORS/prior_d${dchar}.${icnum}"
   (( dn++ ))
done # loop through domains

end_time="$(date +%s)"
length_time=$(( end_time - start_time ))
echo "duration = $length_time"

exit 0

