#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# set echo verbose

# Shell script to work with *syncronous* filter integration       ?
# This script needs to be piped to the filter program with the
# filter namelist async variable set to .true.

echo in submit_cam.csh  >> $1/submit_cam.err
pwd  >> $1/submit_cam.err

set PBS_O_WORKDIR = $1
set PBS_NODEFILE = $2

# If this is first of recursive calls need to get rid of old _ud and _ic files
if ($?First) then
   setenv First no
#  move _ud# files from PBS_O_WORKDIR, where advance_model put them
# see note at _ud below   mv $PBS_O_WORKDIR/assim_model_state_ud* .
   echo 'not first submit_cam' >> $PBS_O_WORKDIR/submit_cam.err
   ls -l >> $PBS_O_WORKDIR/submit_cam.err
else
   setenv First yes
   echo 'first submit_cam' >> $PBS_O_WORKDIR/submit_cam.err
# Clean up any assim_model_ic and ud files 
   rm -f $PBS_O_WORKDIR/assim_model_state_ic*
   rm -f $PBS_O_WORKDIR/assim_model_state_ud*
# Call the model's initialization script to allow it to set up if needed
   csh $PBS_O_WORKDIR/init_advance_model.csh
endif

while(1 == 1)
   rm -f .async_garb
   ls async_may_go > .async_garb
   if ($status == 0) break
   echo waiting_for_async_may_go_file
   sleep 15
end

# if perfect/filter says it's ok to go, then _ic files must have been created;
if (-e assim_model_state_ic1) then
# create 'times' file from DART dates in assim_model_state_ic1
   $PBS_O_WORKDIR/trans_time
# move files to place where all nodes can find them easily
   cp times $PBS_O_WORKDIR
   mv assim_model_state_ic* $PBS_O_WORKDIR
   echo "ls -l PBS_O_WORKDIR/assim " >> $PBS_O_WORKDIR/submit_cam.err
   ls -l $PBS_O_WORKDIR/assim* >> $PBS_O_WORKDIR/submit_cam.err
else
   stop 'no ic file available for trans_time'
endif


# create file to signal status of batch execution of ensemble
echo 'batch not done' > $PBS_O_WORKDIR/batchflag

# First line of filter_control should have number of model states to be integrated
echo filter_control existence  >> $PBS_O_WORKDIR/submit_cam.err
ls -l filter*  >> $PBS_O_WORKDIR/submit_cam.err

set nensmbl = `head -1 filter_control`
set nnodes = `wc -l < $PBS_NODEFILE`
# don't run cam on the filter node (first in list, highest # node)
@ nnodes = $nnodes - 1

echo "nensmbl, nnodes =  $nensmbl , $nnodes"   >> $PBS_O_WORKDIR/submit_cam.err

echo waiting_for_batch_advance_ens 

# batch execution of ensemble
# -----
# figure # batches of CAM runs to do, from # ensemble members and # processors
@ nbatch = $nensmbl / $nnodes
if ($nensmbl % $nnodes != 0 ) @ nbatch++
echo "$nbatch batches will be executed" >> $PBS_O_WORKDIR/submit_cam.err
echo $nbatch batches will be executed
#
set element = 0
set batch = 1
while($batch <= $nbatch)
#   foreach node ( `cat $PBS_NODEFILE` )
   foreach node ( `tail -$nnodes $PBS_NODEFILE` )
      @ element++
      if ($element > $nensmbl) goto all_elements_done
      rsh $node "csh $PBS_O_WORKDIR/advance_model.csh $PBS_O_WORKDIR $element" &
#      echo rsh $node "csh $PBS_O_WORKDIR/advance_model.csh $PBS_O_WORKDIR $element" >> \
#           $PBS_O_WORKDIR/submit_cam.err
      echo rsh $node "csh $PBS_O_WORKDIR/advance_model.csh $PBS_O_WORKDIR $element" >> \
           $PBS_O_WORKDIR/filter_cam.log
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
echo batch_is_done

# These must be moved just before filter takes over *every* iteration
# not just the first.  Doing it at the top of this script for later iterations
# is too late.
# if ($First == yes) mv $PBS_O_WORKDIR/assim_model_state_ud* .
mv $PBS_O_WORKDIR/assim_model_state_ud* .

# Remove the semaphore files
rm -f async_may_go
rm -f $PBS_O_WORKDIR/batchflag

# Cleaned up; let the filter know it can proceed
echo All_done:Please_proceed

# Doing recursive call
# || csh ./async_filter.csh
csh $PBS_O_WORKDIR/submit_cam.csh $PBS_O_WORKDIR $PBS_NODEFILE
