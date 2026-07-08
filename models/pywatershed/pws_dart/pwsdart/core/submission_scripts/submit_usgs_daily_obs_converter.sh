#!/bin/bash

#PBS -N JOB_NAME_TEMPLATE
#PBS -e JOB_NAME_TEMPLATE.stderr
#PBS -o JOB_NAME_TEMPLATE.stdout
#PBS -l PBS_SELECT_TEMPLATE
#PBS -l WALLTIME_TEMPLATE
#PBS -A ACCOUNT_TEMPLATE
#PBS -q QUEUE_TEMPLATE
#PBS -m EMAIL_WHEN_TEMPLATE
#PBS -M EMAIL_WHO_TEMPLATE

module purge
module load ncarenv/23.09
module load craype/2.7.23
module load intel/2023.2.1
module load ncarcompilers/1.0.0
module load cray-mpich/8.1.27
module load hdf5/1.12.2
module load netcdf/4.9.2
module load nco/5.1.9
export PALS_CPU_BIND=none

module load conda
conda activate /glade/work/arezoo/conda-envs/pws

python3 create_usgs_daily_obs_seq.py >& JOB_NAME_TEMPLATE.py.stdeo
rm -rf WAIT_FILE_TEMPLATE

exit 0
