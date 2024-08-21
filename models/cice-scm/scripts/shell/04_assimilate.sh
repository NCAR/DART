#!/bin/bash -l
### Job Name
#PBS -N run_SIassimilations
### Charging account
#PBS -A UWAS0083
### Request one chunk of resources with 1 CPU and 10 GB of memory
#PBS -l select=1:ncpus=1:mem=4GB
### Allow job to run up to 30 minutes
#PBS -l walltime=11:00:00 
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
python ../python/04a_setup_da_case.py SIT_test_bounded mollyw spinup_test open_water atm
python ../python/04b_cycle.py SIT_test_bounded mollyw free_test 3 bounded 2011 1 2 2011 12 31 null SAT_SEAICE_AGREG_THICKNESS
python ../python/05_postprocess.py SIT_test_bounded 3 all
