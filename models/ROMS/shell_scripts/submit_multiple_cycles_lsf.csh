#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
#==========================================================================
#
# This utility will launch a series of dependent jobs for the LSF scheduler to
# accomodate a cycling experiment. Multiple jobs get queued, but only run if
# the previous job completes successfully.
#
# This utility is designed to be run interactively and requires knowledge of
# the EXPERIMENT DIRECTORY, a job name, and how many jobs to submit.

set jobname = roms_dart
set njobs = 4
set rundir = EXPERIMENT_DIRECTORY

cd $rundir

@ n = 1

while ( $n <= $njobs )

  @ nm1 = $n - 1

  set scriptname = `printf multi_cycle_job_%04d.lsf $n`

  set thisjobname = `printf ${jobname}_%04d $n`
  set lastjobname = `printf ${jobname}_%04d $nm1`

  echo '#\!/bin/csh'                   >! $scriptname
  echo "#BSUB -J $thisjobname"        >> $scriptname
  echo "#BSUB -o $thisjobname.%J.log" >> $scriptname
  echo "#BSUB -P P8685nnnn"           >> $scriptname
  echo "#BSUB -q small"               >> $scriptname
  echo "#BSUB -n 16"                  >> $scriptname
  echo "#BSUB -R 'span[ptile=16]'"    >> $scriptname
  echo "#BSUB -W 1:00"                >> $scriptname
  echo "#BSUB -N -u ${USER}@ucar.edu" >> $scriptname

  if ($n > 1) then
     echo '#BSUB -w done("'$lastjobname'")' >> $scriptname
  endif

  echo " "                            >> $scriptname
  echo "cd $rundir"                   >> $scriptname
  echo "./cycle.csh"                  >> $scriptname
  echo " "                            >> $scriptname

  set submissionstring = `bsub < $scriptname`

  # submissionstring is of the form:
  # Job <584064> is submitted to queue <geyser>.

  set job_ID = `echo $submissionstring | grep -oE '[[:digit:]]+'`

  echo "submitted job $n of $njobs. $scriptname $job_ID"

  @ n++
end

exit 0


