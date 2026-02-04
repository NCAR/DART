#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# init_ensemble_var.sh - script that creates perturbed initial
#                         conditions from the WRF-VAR system.
#                         (perts are drawn from the perturbation bank)

set -euo pipefail

initial_date="$1"
paramfile="$(readlink -f "$2")" # Get absolute path for param.sh from command line arg
source "$paramfile"

cd "${RUN_DIR}"

# Generate the i/o lists in rundir automatically when initializing the ensemble
# Required to run filter during assimilation step
num_ens="${NUM_ENS}"      # set from param file
domains="${NUM_DOMAINS}"
dn=1

while (( dn <= domains )); do
     dchar="$(echo "$dn + 100" | bc | cut -b2-3)"
     input_file_name="input_list_d${dchar}.txt"
     input_file_path="./advance_temp"
     output_file_name="output_list_d${dchar}.txt"
     n=1

     [[ -e "$input_file_name"  ]] && rm -f "$input_file_name"
     [[ -e "$output_file_name" ]] && rm -f "$output_file_name"

     while (( n <= num_ens )); do

          ensstring="$(printf %04d "$n")"
          in_file_name="${input_file_path}${n}/wrfinput_d${dchar}"
          out_file_name="filter_restart_d${dchar}.${ensstring}"

          echo "$in_file_name"  >> "$input_file_name"
          echo "$out_file_name" >> "$output_file_name"

          (( n++ ))
     done   # loop through ensemble members
     (( dn++ ))
done   # loop through domains


read -r -a gdate  < <(echo "$initial_date 0h -g" | "${DART_DIR}/models/wrf/work/advance_time")
read -r -a gdatef < <(echo "$initial_date ${ASSIM_INT_HOURS}h -g" | "${DART_DIR}/models/wrf/work/advance_time")
wdate="$(echo "$initial_date 0h -w" | "${DART_DIR}/models/wrf/work/advance_time")"
yyyy="${initial_date:0:4}"
mm="${initial_date:4:2}"
dd="${initial_date:6:2}"
hh="${initial_date:8:2}"
nn="00"
ss="00"

${COPY} "${TEMPLATE_DIR}/namelist.input.meso" namelist.input
${REMOVE} "${RUN_DIR}/WRF"
${LINK} "${OUTPUT_DIR}/${initial_date}" WRF

n=1
while (( n <= NUM_ENS )); do

   echo "  QUEUEING ENSEMBLE MEMBER $n at $(date)"

   mkdir -p "${RUN_DIR}/advance_temp${n}"

   ${LINK} "${RUN_DIR}/WRF_RUN/"* "${RUN_DIR}/advance_temp${n}/."
   ${LINK} "${RUN_DIR}/input.nml" "${RUN_DIR}/advance_temp${n}/input.nml"

   ${REMOVE} script.sed
   cat > script.sed << EOF
/start_year/c\
start_year = ${yyyy},
/start_month/c\
start_month = ${mm},
/start_day/c\
start_day = ${dd},
/start_hour/c\
start_hour = ${hh},
/start_minute/c\
start_minute = ${nn},
/start_second/c\
start_second = ${ss},
/end_year/c\
end_year = ${yyyy},
/end_month/c\
end_month = ${mm},
/end_day/c\
end_day = ${dd},
/end_hour/c\
end_hour = ${hh},
/end_minute/c\
end_minute = ${nn},
/end_second/c\
end_second = ${ss},
/max_dom/c\
max_dom = ${NUM_DOMAINS},
EOF

   sed -f script.sed "${RUN_DIR}/namelist.input" > "${RUN_DIR}/advance_temp${n}/namelist.input"

   ${COPY} "${OUTPUT_DIR}/${initial_date}/wrfinput_d01_${gdate[0]}_${gdate[1]}_mean" \
           "${RUN_DIR}/advance_temp${n}/wrfvar_output.nc"
   sleep 3
   ${COPY} "${RUN_DIR}/add_bank_perts.ncl" "${RUN_DIR}/advance_temp${n}/."

   cmd3="ncl 'MEM_NUM=${n}' 'PERTS_DIR=\"${PERTS_DIR}/work/boundary_perts\"' ${RUN_DIR}/advance_temp${n}/add_bank_perts.ncl"
   ${REMOVE} "${RUN_DIR}/advance_temp${n}/nclrun3.out"
   cat > "${RUN_DIR}/advance_temp${n}/nclrun3.out" << EOF
${cmd3}
EOF

   cat > "${RUN_DIR}/rt_assim_init_${n}.sh" << EOF
#!/bin/bash
#=================================================================
#PBS -N first_advance_${n}
#PBS -j oe
#PBS -A ${COMPUTER_CHARGE_ACCOUNT}
#PBS -l walltime=${ADVANCE_TIME}
#PBS -q ${ADVANCE_QUEUE}
#PBS -l job_priority=${ADVANCE_PRIORITY}
#PBS -m ae
#PBS -M ${EMAIL}
#PBS -k eod
#PBS -l select=${ADVANCE_NODES}:ncpus=${ADVANCE_PROCS}:mpiprocs=${ADVANCE_MPI}
#=================================================================

set -uo pipefail

echo "rt_assim_init_${n}.sh is running in \$(pwd)"

cd "${RUN_DIR}/advance_temp${n}"

if [[ -e wrfvar_output.nc ]]; then
  echo "Running nclrun3.out to create wrfinput_d01 for member ${n} at \$(date)"

  chmod +x nclrun3.out
  ./nclrun3.out > add_perts.out 2>&1

  if [[ ! -s add_perts.err ]]; then
    echo "Perts added to member ${n}"
  else
    echo "ERROR! Non-zero status returned from add_bank_perts.ncl. Check ${RUN_DIR}/advance_temp${n}/add_perts.err."
    cat add_perts.err
    exit 1
  fi

  ${MOVE} wrfvar_output.nc wrfinput_d01

  # For nested domain setups only, downscale perturbations to inner domains
  # For single domains, next section is skipped, no downscaling applied

    dn="${domains}"
    while (( dn > 1 )); do

    # Prep domain files for ndown.exe
    # input files of wrfout_d01_[time] (parent domain), wrfndi_d02 (nested domain)
    # output files of wrfinput_d02 and wrfbdy_d02

    dchar="$(echo "$dn + 100" | bc | cut -b2-3)"

    ${LINK} wrfinput_d01 wrfout_d01_${yyyy}-${mm}-${dd}_${hh}:00:00

    ${COPY} "${OUTPUT_DIR}/${initial_date}/wrfinput_d${dchar}_${gdate[0]}_${gdate[1]}_mean" \
            "${RUN_DIR}/advance_temp${n}/wrfndi_d02"

    echo "Running ndown.exe to downscale perturbed wrfinput_d01 onto wrfinput_d${dchar} for member ${n}"

    # Downscale parent domain to nested domain (wrfinput_d{??})
    mpiexec -n 4 -ppn 4 ./ndown.exe  > ndown_d${dchar}.out

    ${MOVE} wrfinput_d02 wrfinput_d${dchar}
    ${MOVE} wrfbdy_d02   wrfbdy_d${dchar}    
    ${REMOVE} wrfndi_d02
         (( dn-- ))
    done   # loop through domains
  

fi

cd "${RUN_DIR}"

echo "Running first_advance.sh for member ${n} at \$(date)"
"${SHELL_SCRIPTS_DIR}/first_advance.sh" "${initial_date}" "${n}" "${paramfile}"
EOF

   chmod +x "${RUN_DIR}/rt_assim_init_${n}.sh"
   qsub "${RUN_DIR}/rt_assim_init_${n}.sh"

   (( n++ ))

done

exit 0

