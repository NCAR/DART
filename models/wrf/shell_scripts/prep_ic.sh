#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

set -uo pipefail

if [[ $# -gt 0 ]]; then
   n="$1"         # pass in the ensemble member number
   datep="$2"     # needed for correct path to file
   domains="$3"
   paramfile="$4"
else # values come from environment variables 
   n="${mem_num}"
   datep="${date}"
   domains="${domain}"
   paramfile="${paramf}"
fi

source "$paramfile"

echo "prep_ic.sh using n=$n datep=$datep domains=$domains paramfile=$paramfile"
echo "domain 1 using cycle_vars_a, any other nested domains using cycle_vars_b"

dn=1
while (( dn <= domains )); do
  dchar="$(echo "$dn + 100" | bc | cut -b2-3)"

   if (( dn == 1 )); then

      # cycle_vars_a defined in paramfile (bash array)
      num_vars="${#cycle_vars_a[@]}"
      cycle_str=''   # these are variables we want to cycle
      i=0
      while (( i < num_vars-1 )); do
         cycle_str+="${cycle_vars_a[$i]},"
         (( i++ ))
      done
      cycle_str+="${cycle_vars_a[$((num_vars-1))]}"
      echo "${cycle_str}"

   else   # larger (nested) domains can use a different list of cycled variables (e.g. radar)

      # cycle_vars_b defined in paramfile (bash array)
      num_vars="${#cycle_vars_b[@]}"
      cycle_str=''   # these are variables we want to cycle
      i=0
      while (( i < num_vars-1 )); do
         cycle_str+="${cycle_vars_b[$i]},"
         (( i++ ))
      done
      cycle_str+="${cycle_vars_b[$((num_vars-1))]}"
      echo "${cycle_str}"

   fi

   ensstring="$(printf "%04d" "$n")"
   dchar="$(printf "%02d" "$dn")"

   ncks -A -v "${cycle_str}" \
             "${OUTPUT_DIR}/${datep}/PRIORS/prior_d${dchar}.${ensstring}" \
             "${RUN_DIR}/advance_temp${n}/wrfinput_d${dchar}"

   touch "${RUN_DIR}/ic_d${dchar}_${n}_ready"

   (( dn++ ))
done  # loop through domains

exit 0

