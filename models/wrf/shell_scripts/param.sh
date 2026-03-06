#!/bin/bash

# DART software - Copyright UCAR.
# This open source software is provided "as is" without charge.

# ADAPTIVE_INFLATION is disconnected from input.nml
# ASSIM_INT_HOURS is implicit in (ALL) the scripts except assim_advance.sh
# ASSIM_INT_MINUTES support needs to be added to param.sh

# -----------------------------------------------------------
# Environment setup (example for NCAR Derecho)
# -----------------------------------------------------------
module load nco          
module load ncl/6.6.2    
set -uo pipefail
# -----------------------------------------------------------
# Assimilation parameters
# -----------------------------------------------------------
NUM_ENS=20
ASSIM_INT_MINUTES=0    # 0 means use ASSIM_INT_HOURS
ASSIM_INT_HOURS=6      # ignored if ASSIM_INT_MINUTES > 0
IC_PERT_SCALE=0.25
# Set to 1 to enable adaptive inflation
# For pure forecast mode turn off adaptive inflation (set = 0)
ADAPTIVE_INFLATION=1   
NUM_DOMAINS=2          
# -----------------------------------------------------------
# Directory structure
# IMPORTANT: scripts rely on these relative names
# -----------------------------------------------------------
BASE_DIR=/glade/derecho/scratch/bmraczka/WRFv4.5_kansas   
RUN_DIR="${BASE_DIR}/rundir"
TEMPLATE_DIR="${BASE_DIR}/template"
OBSPROC_DIR="${BASE_DIR}/obsproc"
OUTPUT_DIR="${BASE_DIR}/output"
ICBC_DIR="${BASE_DIR}/icbc"
POST_STAGE_DIR="${BASE_DIR}/post"
OBS_DIAG_DIR="${BASE_DIR}/obs_diag"
PERTS_DIR="${BASE_DIR}/perts"
# -----------------------------------------------------------
# Component paths
# -----------------------------------------------------------
SHELL_SCRIPTS_DIR="${BASE_DIR}/scripts"
DART_DIR=/glade/work/bmraczka/DART                                
WRF_DM_SRC_DIR=/glade/work/bmraczka/WRF/WRFv4.5_git               
WPS_SRC_DIR=/glade/work/bmraczka/WRF/WPSv4.5_git                  
VAR_SRC_DIR=/glade/work/bmraczka/WRF/WRFDAv4.5_git               
if [[ ${NUM_DOMAINS} -gt 1 ]]; then
   echo 
   DART_DOM_DIR="${DART_DIR}/models/wrf/tutorial/template_nest"
   echo "NUM_DOMAINS = ${NUM_DOMAINS}" 
   echo "Assigning input.nml.template for multiple WRF domains"
   echo
else
   echo
   DART_DOM_DIR="${DART_DIR}/models/wrf/tutorial/template"
   echo "NUM_DOMAINS = ${NUM_DOMAINS}" 
   echo "Assigning input.nml.template for single WRF domain"
fi

# -----------------------------------------------------------
# Template / IC file sources
# -----------------------------------------------------------
GEO_FILES_DIR=/glade/u/home/wrfhelp/WPS_GEOG                      
GRIB_DATA_DIR="${ICBC_DIR}/grib_data"                             
GRIB_SRC='GFS'

# -----------------------------------------------------------
# Variable lists for extraction/cycling
# -----------------------------------------------------------

extract_vars_a=( U V PH THM MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN \
                 U10 V10 T2 Q2 PSFC TSLB SMOIS TSK RAINC RAINNC GRAUPELNC )

extract_vars_b=( U V W PH THM MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN \
                 U10 V10 T2 Q2 PSFC TSLB SMOIS TSK RAINC RAINNC GRAUPELNC )

cycle_vars_a=( U V PH THM MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN \
               U10 V10 T2 Q2 PSFC TSLB SMOIS TSK )

cycle_vars_b=( U V PH THM MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN \
               U10 V10 T2 Q2 PSFC TSLB SMOIS TSK )

increment_vars_a=( U V PH THM MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN U10 V10 T2 Q2 PSFC )
increment_vars_b=( U V PH THM MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN U10 V10 T2 Q2 PSFC )

# -----------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------
# OBS_VERIF_HRS=$(( ASSIM_INT_HOURS / 2 ))

# -----------------------------------------------------------
# Queueing / HPC system settings
# -----------------------------------------------------------
SUPER_PLATFORM="derecho"
COMPUTER_CHARGE_ACCOUNT=P86850054     
EMAIL="bmraczka@ucar.edu"            

if [[ "$SUPER_PLATFORM" == "derecho" ]]; then

    # PBS queue example
    FILTER_QUEUE="main"
    FILTER_PRIORITY="premium"
    FILTER_TIME="0:35:00"
    FILTER_NODES=2
    FILTER_PROCS=128
    FILTER_MPI=128

    ADVANCE_QUEUE="main"
    ADVANCE_PRIORITY="premium"
    ADVANCE_TIME="0:20:00"
    ADVANCE_NODES=1
    ADVANCE_PROCS=128
    ADVANCE_MPI=128

else

    # LSF/SLURM example values
    FILTER_QUEUE="regular"
    FILTER_TIME="0:25"
    FILTER_CORES=512

    ADVANCE_QUEUE="regular"
    ADVANCE_TIME="0:18"
    ADVANCE_CORES=64

    FILTER_PTILE=16
    ADVANCE_PTILE=16
fi

# -----------------------------------------------------------
# System commands
# -----------------------------------------------------------
export REMOVE='rm -rf'
export COPY='cp -pfr'
export MOVE='mv -f'
export LINK='ln -fs'
export WGET='/usr/bin/wget'
export LIST='ls'

return

