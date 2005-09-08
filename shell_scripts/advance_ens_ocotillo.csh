#!/bin/csh -f
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# $Id$
#
#   Multi-processor jobs must be submitted as batch, under PBS
### Job name
#PBS -N Exp32
### Declare job non-rerunable
#PBS -r n
### Output files
#PBS -e Exp32.err
#PBS -o Exp32.log
### Queue name (small(20min), medium(2hr), long(12hr), verylong(72hr))
#PBS -q medium
# #PBS -l nodes=5
# perfect_model_obs:
#PBS -l nodes=1

### This job's working directory; must cd to it, or it will run in /home...
if ($?PBS_O_WORKDIR) then
   cd $PBS_O_WORKDIR
else
   setenv PBS_O_WORKDIR `pwd`
endif
### Output to confirm job characteristics
if ($?PBS_JOBNAME) then
   echo "Running $PBS_JOBNAME on host "`hostname`
else
   echo "Running on host "`hostname`
endif
echo "Time is "`date`
echo "Directory is "`pwd`
echo "This job runs on the following processors:"

# First line of filter_control should have number of model states to be 
# integrated
set nensmbl = `head -1 filter_control`

### Define number of processors; # of lines in PBS_NODEFILE
if ($?PBS_NODEFILE) then
   setenv NPROCS `wc -l < $PBS_NODEFILE`
else
   setenv PBS_NODEFILE nodefile
   rm -f $PBS_NODEFILE
   set startnode = 1
   set endnode = 14
   set inode = $startnode
   set NPROCS = 28
   set iproc = 1
   while($iproc <= $NPROCS & $iproc <= $nensmbl)
      echo node$inode >> $PBS_NODEFILE
      if ($inode == $endnode) then
         set inode = $startnode
      else
         @ inode ++
      endif
      @ iproc ++
   end
endif
cat "$PBS_NODEFILE"

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

      rsh $node "csh $PBS_O_WORKDIR/advance_model.csh $PBS_O_WORKDIR $element /var/tmp/${USER}_tempdir${element}" &

      sleep 0.1

   end
# Another way to monitor progress.  batchflag has other info to start,
# so this echo can be removed and scripts will still work.
      echo waiting to finish batch $batch of $nbatch >> $PBS_O_WORKDIR/batchflag
      wait
   @ batch++
end
all_elements_done:

# Wait for all *background* processes to finish up
wait

# signal to async_filter.csh (if async=1) or to Aadvance_state (if async=2) to continue
rm -f $PBS_O_WORKDIR/batchflag
