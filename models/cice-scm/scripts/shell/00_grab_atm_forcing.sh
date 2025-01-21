#!/bin/bash -l
### Job Name
#PBS -N extract_atmospheric_forcings
### Charging account
#PBS -A P93300065
### Request one chunk of resources with 1 CPU and 10 GB of memory
#PBS -l select=1:ncpus=1:mem=200GB
### Set job walltime
#PBS -l walltime=1:00:00
#PBS -l job_priority=economy
### Route the job to the economy queue
#PBS -q main
### Join output and error streams into single file
#PBS -j oe
### send emails on abort and exit
#PBS -m ae
#PBS -M mmw906@uw.edu

export TMPDIR=/glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

### Load Python module and activate NPL environment
# module load ncarenv python
conda activate cice-scm-da

## If you want spinup forcings (from 2000-2010) uncomment the next line...
python ../python/00_grab_atm_forcing.py 2000 2010 spinup 30
## If you want free/assimilation forcings (from 2011-2015) uncomment the next line...
# python ../python/00_grab_atm_forcing.py 2011 2015 free 30
## If you want forecast forcings (from 2012-2016) uncomment the next line...
# python ../python/00_grab_atm_forcing.py 2012 2016 branch_2012 30