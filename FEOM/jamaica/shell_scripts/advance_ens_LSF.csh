#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$

# The number of processors (NPROCS) is determined by one of two things.
# If this script is run interactively, NPROCS is unity.
# If this script is run in batch mode, NPROCS is user-defined by  
# the argument to 'bsub -n xxxx'
#
# Initial version to run on lightning IBM Linux cluster

### Job name
#BSUB -J advance_ens
### Declare job non-rerunable (default behavior with BSUB?)

### Output files
#BSUB -o adv_ens.%J.log
#BSUB -P 86850054
### Queue charging    cheapest ..... most expensive   (at least on lightning)
### Queue name (standby, economy, [regular,debug], premium)
#BSUB -q regular
#BSUB -n 28

# First line of filter_control should have number of model states to be integrated
set nensmbl = `head -1 filter_control`

# Determine number of processors
#
# list of hosts/machines is in $PROCNAMES
# the quoting is VERY IMPORTANT for PROCNAMES

if ($?LSB_HOSTS) then                       ;# batch
   set NPROCS = `echo $LSB_HOSTS | wc -w`
#   set PROCNAMES = "$LSB_HOSTS"
   set PROCNAMES = ($LSB_HOSTS)
else                                        ;# interactive
   set NPROCS = 1
   set PROCNAMES = $host
endif

### This job's working directory; must cd to it, or it will run in /home...
if ($?LS_SUBCWD) then
   cd $LS_SUBCWD
else
   setenv LS_SUBCWD `pwd`
endif

### Output to confirm job characteristics
if ($?LSB_JOBNAME) then
   echo Running $LSB_JOBNAME on host `hostname`
else
   echo "Running on host "`hostname`
endif
echo Time is `date`
echo Directory is `pwd`
echo This job runs on the following nodes:
echo $PROCNAMES

echo This job has allocated $NPROCS processors

# figure # batches of runs to do, from # ensemble members and # processors
@ nbatch = $nensmbl / $NPROCS
if ($nensmbl % $NPROCS != 0 ) @ nbatch++
echo $nbatch batches will be executed

# Send jobs to nodes
set element = 0
set batch = 1
while($batch <= $nbatch)
   foreach proc ( $PROCNAMES )
      @ element++
      if ($element > $nensmbl) goto all_elements_done

      ssh $proc "csh $LS_SUBCWD/advance_model.csh $LS_SUBCWD $element /ptmp/${user}/tmp$user$element " &

   end
# Another way to monitor progress.  batchflag has other info to start,
# so this echo can be removed and scripts will still work.
   echo waiting to finish batch $batch of $nbatch >> $LS_SUBCWD/batchflag
   wait
   @ batch++
end
all_elements_done:

# Wait for all *background* processes to finish up
wait

# Attempt to rerun members that did not advance successfully.
set rerun = ' '
set nrerun = 0
set goodprocs = ' '
set badprocs = ' '
set element = 0
set batch = 1
set NPROCS = 0
while($batch <= $nbatch)
   set iproc = 1
   foreach proc ( $PROCNAMES )
      @ element++
      if ($element > $nensmbl) goto all_elements_checked
      set iblo = `expr $element \+ 10000`
      set iblo = `echo $iblo | cut -c2-5`
      set blown = `grep $iblo blown_*.out | cat | wc -l`
      if ($blown == 0 && -e $LS_SUBCWD/assim_model_state_ud$element) then
#      if ($blown == 0 && -e /ptmp/${user}/tmp${user}${element}/dart_wrf_vector) then
         set goodprocs = ($goodprocs $PROCNAMES[$iproc])
         @ NPROCS++
      else
         set badprocs = ($badprocs $PROCNAMES[$iproc])
         set rerun = ($rerun $element)
         @ nrerun++
      endif
      @ iproc++
   end
   @ batch++
end
all_elements_checked:

rm -f blown_*.out

if ($nrerun > 0) then
   echo $nrerun members will be rerun
   @ nbatch = $nrerun / $NPROCS
   if ($nrerun % $NPROCS != 0 ) @ nbatch++
   echo $nbatch batches will be executed
   set element = 0
   set batch = 1
   while($batch <= $nbatch)
      foreach proc ( $goodprocs )
         @ element++
         if ($element > $nrerun) goto all_elements_rerun

         ssh $proc "csh $LS_SUBCWD/advance_model.csh $LS_SUBCWD $rerun[$element] /ptmp/${user}/tmp${user}$rerun[$element] " &

      end
# Another way to monitor progress.  batchflag has other info to start,
# so this echo can be removed and scripts will still work.
      echo waiting to finish rerun batch $batch of $nbatch >> $LS_SUBCWD/batchflag
      wait
      @ batch++
   end
all_elements_rerun:
endif

# Wait for all *background* processes to finish up
wait

mkdir -p FAILURES
mv blown_*.out FAILURES/

# signal to async_filter.csh (if async=1) or to Aadvance_state (if async=2) to continue
rm -f $LS_SUBCWD/batchflag
