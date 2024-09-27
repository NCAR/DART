#!/bin/bash -l
### Job Name
#PBS -N run_SIassimilations
### Charging account
#PBS -A UWAS0083
### Request one chunk of resources with 1 CPU and 10 GB of memory
#PBS -l select=1:ncpus=1:mem=4GB
### Allow job to run up to 30 minutes
#PBS -l walltime=00:50:00 
#PBS -l job_priority=economy
### Route the job to the casper queue
#PBS -q main
### Join output and error streams into single file
#PBS -j oe
#PBS -m ae
#PBS -M mmw906@uw.edu

export TMPDIR=/glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

### Load Python module and activate NPL environment
conda activate cice-scm-da

### Run analysis script
python ../python/04a_setup_1day_da_case.py FB_Barents_May mollyw spinup_Barents uniform_ice 2011 5 15 atm
python ../python/04b_cycle.py FB_Barents_May mollyw free_Barents 13 bounded 2011 5 15 2011 5 15 null SAT_SEAICE_AGREG_FREEBOARD
# python ../python/05_postprocess.py SIT_Barents_Jan mollyw 3 all
