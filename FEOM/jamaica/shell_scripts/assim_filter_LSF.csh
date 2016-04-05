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

# Initial version to run on lightning IBM Linux cluster

### Job name
#BSUB -J assim_filter
### Declare job non-rerunable (default behavior with BSUB?)

### Output files
#BSUB -o assim_filter.%J.o
### Queue name (economy, regular, premium)
#BSUB -q regular
set NPROCS = 9
#BSUB -n 9

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
echo $LSB_HOSTS
echo This job has allocated $NPROCS processors

# First line of assim_region_control should have number of regions to be assimilated
set nregions = `head -1 assim_region_control`

# figure # batches of runs to do, from # regions and # processors
@ nbatch = $nregions / $NPROCS
if ($nregions % $NPROCS != 0 ) @ nbatch++
echo $nbatch batches will be executed

# Send jobs to nodes
set element = 0
set batch = 1
while($batch <= $nbatch)
   foreach node ( $LSB_HOSTS )
      @ element++
      if ($element > $nregions) goto all_elements_done

      ssh $node "csh $LS_SUBCWD/assim_region.csh $LS_SUBCWD $element /ptmp/${user}/tmp$user$element " &

   end
# Another way to monitor progress.  batchflag has other info to start,
# so this echo can be removed and scripts will still work.
   echo waiting to finish batch $batch  >> $LS_SUBCWD/batchflag
   wait
   @ batch++
end
all_elements_done:

# Wait for all *background* processes to finish up
wait

# signal to filter_assim to continue
rm -f $LS_SUBCWD/batchflag
