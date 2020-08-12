#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DESCRIPTION:
#
# Generate daily streamflow observation sequence files.
# Submit this script
# to your batch system and it will invoke the 'makedaily.csh'
# script once for each conversion day.
# this one does N conversions in parallel from a command script.
#
#==========================================================================
# PBS directives                qsub paralell_daily.batch
#
# qstat    information on running jobs
# qsub     submit a job
# qdel     kill a job
# qpeek    see output from a running job
#

#PBS -N JOB_NAME_TEMPLATE
#PBS -e JOB_NAME_TEMPLATE.stderr
#PBS -o JOB_NAME_TEMPLATE.stdout
#PBS -A ACCOUNT_TEMPLATE
#PBS -m EMAIL_WHEN_TEMPLATE
#PBS -M EMAIL_WHO_TEMPLATE
#PBS -q QUEUE_TEMPLATE
#PBS -l WALLTIME_TEMPLATE

# Do the batch processing on just one cheyenne node.
#PBS -l select=1:ncpus=36:mpiprocs=36

#
#==========================================================================
# Everything establish in experiment YAML namelist file under
# observation_preparation: USGS_daily

# set things that vary between batch systems here.

if ($?SLURM_JOB_ID) then

  echo 'running SLURM'
  setenv SLURM true
  source /glade/u/apps/opt/slurm_init/init.csh
  setenv RUNCMD "srun --multi-prog"
  setenv JOBNAME $SLURM_JOB_NAME

  # should match the mpiprocs=X setting above
  # todo ... there is an enviroment variable for this ...
  set njobs = 36

else if ($?PBS_NODEFILE) then

  echo 'running PBS'
  setenv PBS true
  setenv MPI_SHEPHERD true
  setenv RUNCMD "mpiexec_mpt launch_cf.sh"
  setenv JOBNAME $PBS_JOBNAME

  set njobs = `cat $PBS_NODEFILE | wc -l`
  echo "njobs from PBS_NODEFILE is $njobs"

else if ($?LSB_HOSTS) then

  echo 'running LSF'
  setenv LSF true
  setenv MP_PGMMODEL mpmd
  setenv RUNCMD "mpirun.lsf -cmdfile"
  setenv JOBNAME $LSB_OUTPUTFILE:ar
  # should match the mpiprocs=X setting above
  # todo ... there is an enviroment variable for this ...
  set njobs = 36

else

  echo 'running without a batch system'
  setenv NOBATCH true
  setenv RUNCMD "csh"
  setenv JOBNAME "streamflow"
  set njobs = 36

endif

setenv TMPDIR /glade/scratch/$USER/temp
mkdir -p $TMPDIR

echo job started at `date`

# The day_cmdfiles are generated in python in setup_usgs_daily.py

foreach cmdfile ( `ls day_cmdfile.*` )
 
   $RUNCMD ./${cmdfile} >& output_${cmdfile}
   mv $cmdfile done_${cmdfile}
   
end

# This file maybe created by the submitting process.... 
rm -f .parallel_daily_batch_not_complete

echo job ended at `date`

exit 0


