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

# Get rid of any pre-existing go_advance or go_filter files
# Need to guarantee this gets done before filter gets to advance
rm -f go_advance_model
rm -f go_end_filter
rm -f go_assim_regions

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


# Hang around forever for now and wait for go_advance file to appear
while(1 == 1)
   # If go_end_filter exists then stop this process
   ls go_end_filter > .option3_garb
   if($status == 0) exit
 

#--------------------------------------------------------------------------------------- 
   # Check to see if the go_advance_model file exists
   ls go_advance_model > .option3_garb
   if($status == 0) then


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
      rm -f $PBS_O_WORKDIR/batchflag

      # finished with advance_model so remove the go_advance_model file
      rm -f go_advance_model
   endif

#--------------------------------------------------------------------------------------- 
   # Check to see if the go_advance_model file exists
   ls go_assim_regions > .option3_garb
   if($status == 0) then

# First line of assim_region_control should have number of regions to be assimilated
set nregions = `head -1 assim_region_control`

# figure # batches of runs to do, from # ensemble members and # processors
@ nbatch = $nregions / $NPROCS
if ($nregions % $NPROCS != 0 ) @ nbatch++
echo $nbatch batches will be executed
#
# Create a directory for each member to run in for namelists
set element = 0
set batch = 1
while($batch <= $nbatch)
   foreach node ( `cat $PBS_NODEFILE` )
      @ element++
      if ($element > $nregions) goto all_regions_done

      $PBS_O_WORKDIR/assim_region.csh $PBS_O_WORKDIR $element ${USER}_tempdir${element} &

      sleep 0.1

   end
# Another way to monitor progress.  batchflag has other info to start,
# so this echo can be removed and scripts will still work.
      echo waiting to finish batch $batch  >> $PBS_O_WORKDIR/batchflag
      wait

      sleep 1

   @ batch++
end
all_regions_done:

# Wait for all *background* processes to finish up
wait

#
# signal to async_filter.csh to continue
rm -f $PBS_O_WORKDIR/batchflag





      rm -f go_assim_regions

   endif


#--------------------------------------------------------------------------------------- 

# No files found, wait and check again
   sleep 1
end

