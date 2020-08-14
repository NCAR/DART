#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
#==========================================================================
#
# This utility launches a series of dependent jobs for the SLURM scheduler to
# accomodate a cycling experiment. Multiple jobs get queued, but only run if
# the previous job completes successfully.
#
# The first dependent job is a job array where each model advance is run in
# a separate job. When ALL of those jobs (i.e. the entire job array) is
# successfully complete, the assimilation job starts. When that job finishes,
# the next job array of advances starts ... and so on.
#
# The resources for the model advance jobs are specified in advance_ensemble.csh
# The resources for the assimilation job is specified in run_filter.csh
# Each resource can be tailored separately for maximum efficiency.
#
# This utility is designed to be run interactively and requires knowledge of
# the EXPERIMENT DIRECTORY, a job name, and how many jobs to submit.

set NCYCLES = 11
set rundir = EXPERIMENT_DIRECTORY

cd $rundir || exit 1

#--------------------------------------------------------------------------
# Overall strategy is to fire off a series of dependent jobs.
# Successful completion of the first filter job will free the queued model
# advances. That successful completion of that job array will free the
# next filter job ... and so on.
#--------------------------------------------------------------------------

set depstr = " "
set i = 1

while ( $i <= $NCYCLES )

   echo "queueing cycle $i"

   #-----------------------------------------------------------------------
   # launch job array of ensemble advances and capture job ID for dependency
   #-----------------------------------------------------------------------

   set submissionstring = `sbatch $depstr ./advance_ensemble.csh`
   set ensjob = `echo $submissionstring | awk '{print($4)}'`
   set depstr = "--dependency=afterok:$ensjob"

   #-----------------------------------------------------------------------
   # run Filter to generate the analysis
   #-----------------------------------------------------------------------

   set submissionstring = `sbatch $depstr ./run_filter.csh`
   set dajob = `echo $submissionstring | awk '{print($4)}'`
   set depstr = "--dependency=afterok:$dajob"

   @ i++

end

exit 0

