#!/bin/tcsh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# $Id$
#

### This job's working directory; must cd to it, or it will run in /home...
setenv PBS_O_WORKDIR `pwd`

### Output to confirm job characteristics
echo "Running on host "`hostname`
echo "Time is "`date`
echo "Directory is "`pwd`
echo "This job runs on the following processors:"

### Define number of processors; # of lines in PBS_NODEFILE
   setenv PBS_NODEFILE nodefile
   rm -f $PBS_NODEFILE
   set NPROCS = 4
   set iproc = 1
   while($iproc <= $NPROCS)
      echo proc$iproc >> $PBS_NODEFILE
      @ iproc ++
   end
cat "$PBS_NODEFILE"
echo This_job_has_allocated $NPROCS nodes

# First line of filter_control should have number of model states to be 
# integrated
set nensmbl = `head -1 filter_control`

# figure # batches of runs to do, from # ensemble members and # processors
@ nbatch = $nensmbl / $NPROCS
if ($nensmbl % $NPROCS != 0 ) @ nbatch++
echo $nbatch batches will be executed
#
# Create a directory for each member to run in for namelists
set element = 0
set batch = 1
while($batch <= $nbatch)
   foreach node ( `cat $PBS_NODEFILE` )
      @ element++
      if ($element > $nensmbl) goto all_elements_done

      $PBS_O_WORKDIR/advance_model.csh $PBS_O_WORKDIR $element ${USER}_tempdir${element} &

      sleep 0.1

   end
# Another way to monitor progress.  batchflag has other info to start,
# so this echo can be removed and scripts will still work.
      echo waiting to finish batch $batch  >> $PBS_O_WORKDIR/batchflag
      wait

      sleep 1

   @ batch++
end
all_elements_done:

# Wait for all *background* processes to finish up
wait

# 
# signal to async_filter.csh to continue
rm -f $PBS_O_WORKDIR/batchflag
