#!/bin/tcsh
# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$

# Modified to use /scratch/local/raeder instead of /scratch/cluster for CAM runs
# See #sl
# Modified to use actual times/dates from DART; see #td

#   Multi-processor jobs must be submitted as batch, under PBS
### Job name
#PBS -N dart_cam
### Declare job non-rerunable
#PBS -r n
### Output files
#PBS -e dart_cam.err
#PBS -o dart_cam.log
### Queue name (small, medium, long, verylong)
#PBS -q small
# Mark Moore recommended 4 as the max # of nodes to use, but might allow more
# for short durations.
# filter; as many as available
#PBS -l nodes=2
# perfect_model_obs:
# #PBS -l nodes=1
### This job's working directory; must cd to it, or it will run in /home...
cd $PBS_O_WORKDIR
### Output to confirm job characteristics
echo Running $PBS_JOBNAME on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo This job runs on the following processors:
### Define number of processors; # of lines in PBS_NODEFILE
setenv NPROCS `wc -l < $PBS_NODEFILE`
cat "$PBS_NODEFILE"
echo This job has allocated $NPROCS nodes

#td
# create 'times' file from DART dates in assim_model_state_ic1
$PBS_O_WORKDIR/trans_time

# First line of filter_control should have number of model states to be integrated
set nensmbl = `head -1 filter_control`

# figure # batches of CAM runs to do, from # ensemble members and # processors
@ nbatch = $nensmbl / $NPROCS
if ($nensmbl % $NPROCS != 0 ) @ nbatch++
echo $nbatch batches will be executed
#
# Create a directory for each member to run in for namelists
# || set element = 1
# || while($element <= $nensmbl)
set element = 0
set batch = 1
while($batch <= $nbatch)
   foreach node ( `cat $PBS_NODEFILE` )
      @ element++
      if ($element > $nensmbl) goto all_elements_done
# Make a temporary directory for this elements run
#sl      mkdir tempdir$element
# Copy the appropriate ics file to the temp directory
#sl      mv assim_model_state_ic$element tempdir$element/temp_ic
# Change to temp directory and run integrate
# || ;must be done as part of rsh below
#      cd tempdir$element
# || packaging;
#    element is not available to advance_model; do multiple commands in 1 rsh 
#       (if there are too many; extract them to a separate script)
#    "" allows expansion of $element and $PBS... here, rather than on node
#sl         "cd $PBS_O_WORKDIR/tempdir$element; \
#sl          csh ../advance_model.csh >&! ../$element.out" &
      rsh $node "csh $PBS_O_WORKDIR/advance_model.csh $PBS_O_WORKDIR $element" &
#    echo ... could be removed, unless we want to keep track of progress
      echo rsh $node \
         "csh $PBS_O_WORKDIR/advance_model.csh $PBS_O_WORKDIR $element" &
#sl         "cd $PBS_O_WORKDIR/tempdir$element; \
#sl         "cd $PBS_O_WORKDIR/tempdir$element; \
#sl          csh ../advance_model.csh >&! ../$element.out" &
# Pop back up to main directory
#      cd $PBS_O_WORKDIR
# ||      cd ..
   end
# Another way to monitor progress.  batchflag has other info to start,
# so this echo can be removed and scripts will still work.
      echo waiting to finish batch $batch  >> $PBS_O_WORKDIR/batchflag
      wait
   @ batch++
end
all_elements_done:

# Wait for all *background* processes to finish up
wait

# ?   Should we use the disc associated with each node (tempdir# not needed?)
# ?   instead of all tempdirs on same big disc with diff names
#        /scratch/local has 14Gb available, none of it's being used
#        /scratch/cluster has 110Gb, 81% used
# ?   Should this loop be part of batch loop above; reduce number of tempdirs
#     active at any time.
#sl set element = 1
#sl while($element <= $nensmbl)
#sl   mv tempdir$element/temp_ud assim_model_state_ud$element
#sl    rm -rf tempdir$element
#sl    @ element++
#sl end
# signal to sync_submit.csh to continue
rm -f $PBS_O_WORKDIR/batchflag
