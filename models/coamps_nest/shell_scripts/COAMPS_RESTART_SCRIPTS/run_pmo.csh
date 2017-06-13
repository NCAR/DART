#!/bin/tcsh 
#
# SCRIPT:	run_pmo.csh
# AUTHOR:	T. R. Whitcomb
#           Naval Research Laboratory
#
# Runs perfect_model_obs as a PBS job
######
#PBS -N perf_model_obs
#PBS -r n
#PBS -e perf_model_obs.err
#PBS -o perf_model_obs.out
#PBS -q long
#PBS -l nodes=4:ppn=2

# Import the job-specific resource commands
source ${PBS_O_WORKDIR}/job_setup.csh

cd $PBS_O_WORKDIR

set LOGFILE = pmo.dump

./perfect_model_obs | tee ${LOGFILE}

# Once the job completes, we can eliminate this file as we'll get a copy
# of standard output in perf_model_obs.out
rm ${LOGFILE}

