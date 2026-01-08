#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# Script to save a 'perturbation bank' generated from the WRFDA CV3 option
#
# provide the following:
# 	namelist.input
# 	wrfinput_d01
# 	ensemble size
# 	list of perturbed variables
# 	wrfda executable and be.dat

set -uo pipefail

datea=2024051900   # need to start from a known valid date matching the wrfinput_d01 date
paramfile="/glade/derecho/scratch/bmraczka/WRFv4.5_nested_bash/scripts/param.sh"

echo "Sourcing parameter file"
source "$paramfile"

mkdir -p "${PERTS_DIR}/work/boundary_perts"
# this has all wrf and wrfda executables and support files
# /glade/work/bmraczka/WRF/WRFDAv4.5_git
wrfda_dir="${RUN_DIR}/WRF_RUN"    # set this appropriately
work_dir="${PERTS_DIR}/work"      # set this appropriately
# Put the final eperturbation files here for later use
save_dir="${PERTS_DIR}/work/boundary_perts"  # set this appropriately

# These scale variables are not used directly  -- default is to use hard coded values within
# wrfvar section of namelist.input.3dvar file
IC_PERT_SCALE=0.1
IC_HORIZ_SCALE=0.8
IC_VERT_SCALE=0.8


# Number of perturbations to generate, suggest 3-4X the ensemble size. 
# Recommended to test single member first to confirm functionality and desired performance.

num_ens=150


# Get wrfinput_d01 file directly from the 'mean' file generated from real.exe during gen_retro_icbc.sh
# set wrfin_dir = ${work_dir}/wrfin
ASSIM_INT_HOURS=6

module load nco

cd "$work_dir" || exit 1
cp "${TEMPLATE_DIR}/input.nml.template" input.nml

# get a wrfdate and parse
read -r -a gdate < <(echo "$datea 0h -g" | "${DART_DIR}/models/wrf/work/advance_time")
read -r -a gdatef < <(echo "$datea ${ASSIM_INT_HOURS}h -g" | "${DART_DIR}/models/wrf/work/advance_time")
wdate=$(echo "$datea 0h -w" | "${DART_DIR}/models/wrf/work/advance_time")

yyyy="${datea:0:4}"
mm="${datea:4:2}"
dd="${datea:6:2}"
hh="${datea:8:2}"

for ((n=1; n<=num_ens; n++)); do

    mkdir -p "${work_dir}/mem_${n}"
    cd "${work_dir}/mem_${n}" || exit 1
    cp "${wrfda_dir}"/* "${work_dir}/mem_${n}/"

    ln -sf "${OUTPUT_DIR}/${datea}/wrfinput_d01_${gdate[0]}_${gdate[1]}_mean" "${work_dir}/mem_${n}/fg"

    seed_array2=$((n*10))

    cat > script.sed << EOF
/run_hours/c\
run_hours = 0,
   /run_minutes/c\
   run_minutes                = 0,
   /run_seconds/c\
   run_seconds                = 0,
   /start_year/c\
   start_year                 = 1*${yyyy},
   /start_month/c\
   start_month                = 1*${mm},
   /start_day/c\
   start_day                  = 1*${dd},
   /start_hour/c\
   start_hour                 = 1*${hh},
   /start_minute/c\
   start_minute               = 1*00,
   /start_second/c\
   start_second               = 1*00,
   /end_year/c\
   end_year                   = 1*${yyyy},
   /end_month/c\
   end_month                  = 1*${mm},
   /end_day/c\
   end_day                    = 1*${dd},
   /end_hour/c\
   end_hour                   = 1*${hh},
   /end_minute/c\
   end_minute                 = 1*00,
   /end_second/c\
   end_second                 = 1*00,
   /analysis_date/c\
   analysis_date = \'${wdate}.0000\',
   s/PERT_SCALING/${IC_PERT_SCALE}/
   s/HORIZ_SCALE/${IC_HORIZ_SCALE}/
   s/VERT_SCALE/${IC_VERT_SCALE}/
   /seed_array1/c\
   seed_array1 = ${datea},
/seed_array2/c\
seed_array2 = ${seed_array2} /
EOF

   # namelist.input.3dvar must be set for single parent domain to work with gen_pert_bank.sh, 
   # contains all namelist options wrfvar1-14

   sed -f script.sed "${TEMPLATE_DIR}/namelist.input.3dvar" > "${work_dir}/mem_${n}/namelist.input"

    # Create PBS script
    cat > "${work_dir}/mem_${n}/gen_pert_${n}.sh" << EOF
#!/bin/sh
#=================================================================
#PBS -N gen_pert_bank_mem${n}
#PBS -j oe
#PBS -A ${COMPUTER_CHARGE_ACCOUNT}
#PBS -l walltime=0:05:00
#PBS -q ${ADVANCE_QUEUE}
#PBS -l job_priority=${ADVANCE_PRIORITY}
#PBS -o gen_pert_bank_mem${n}.out
#PBS -l select=1:ncpus=4:mpiprocs=4
#PBS -k eod
#PBS -V
#=================================================================

cd "${work_dir}/mem_${n}" || exit 3

mpiexec -n 1 -ppn 1 ./da_wrfvar.exe >& output.wrfvar
mv wrfvar_output wrfinput_d01

# Extract only the fields that are updated by wrfvar, then diff to generate the pert file for this member
ncks -h -F -A -a -v U,V,THM,QVAPOR,MU fg orig_data.nc
ncks -h -F -A -a -v U,V,THM,QVAPOR,MU wrfinput_d01 pert_data.nc
ncdiff pert_data.nc orig_data.nc "pert_bank_mem_${n}.nc"
mv "pert_bank_mem_${n}.nc" "${save_dir}/pert_bank_mem_${n}.nc"
EOF

    qsub "${work_dir}/mem_${n}/gen_pert_${n}.sh"

done

exit 0

