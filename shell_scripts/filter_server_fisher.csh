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

rm -f filter_server.log

### Output to confirm job characteristics
echo "filter_server Running on host "`hostname` >> filter_server.log
echo "Initialized at " `date` >> filter_server.log
echo "This job runs on the following processors:" >> filter_server.log

### Define number of processors; # of lines in PBS_NODEFILE
setenv PBS_NODEFILE nodefile
rm -f $PBS_NODEFILE
set NPROCS = 3
set iproc = 1
while($iproc <= $NPROCS)
   echo proc$iproc >> $PBS_NODEFILE
   echo proc$iproc >> filter_server.log
   @ iproc ++
end
echo This_job_has_allocated $NPROCS nodes >> filter_server.log

# Hang around forever for now and wait for go_advance file to appear
while(1 == 1)
   # If go_end_filter exists then stop this process
   ls go_end_filter > .option3_garb
   if($status == 0) then
      echo "terminating normally at " `date` >> filter_server.log
      exit
   endif
 

#--------------------------------------------------------------------------------------- 
   # Check to see if the go_advance_model file exists
   ls go_advance_model > .option3_garb
   if($status == 0) then


      # First line of filter_control should have number of model states to be integrated
      set nensmbl = `head -1 filter_control`
      echo "advancing $nensmbl members at " `date`>> filter_server.log
     
      # figure # batches of runs to do, from # ensemble members and # processors
      @ nbatch = $nensmbl / $NPROCS
      if ($nensmbl % $NPROCS != 0 ) @ nbatch++
      echo $nbatch batches will be executed >> filter_server.log

      # Create a directory for each member to run in for namelists
      set element = 0
      set batch = 1
      while($batch <= $nbatch)
         foreach node ( `cat $PBS_NODEFILE` )
            @ element++
            if ($element > $nensmbl) goto all_elements_done
            $PBS_O_WORKDIR/advance_model.csh $PBS_O_WORKDIR $element ${USER}_tempdir${element} &

         end
         echo "waiting to finish batch $batch at" `date`  >> filter_server.log
         wait

         @ batch++
      end
      all_elements_done:

      # Need to have all backgrounds completed before continuing
      wait

      # finished with advance_model so remove the go_advance_model file
      echo "Completed this advance at " `date` >> filter_server.log
      echo --------- >> filter_server.log
      rm -f go_advance_model
   endif

#--------------------------------------------------------------------------------------- 
   # Check to see if the go_advance_model file exists
   ls go_assim_regions > .option3_garb
   if($status == 0) then

      # First line of assim_region_control should have number of regions to be assimilated
      set nregions = `head -1 assim_region_control`
      echo "assimilating $nregions regions at " `date`>> filter_server.log

      # figure # batches of runs to do, from # ensemble members and # processors
      @ nbatch = $nregions / $NPROCS
      if ($nregions % $NPROCS != 0 ) @ nbatch++
      echo $nbatch batches will be executed >> filter_server.log
      #
      # Create a directory for each member to run in for namelists
      set element = 0
      set batch = 1
      while($batch <= $nbatch)
         foreach node ( `cat $PBS_NODEFILE` )
            @ element++
            if ($element > $nregions) goto all_regions_done

            $PBS_O_WORKDIR/assim_region.csh $PBS_O_WORKDIR $element ${USER}_tempdir${element} &

         end
         echo "waiting to finish batch $batch" `date`  >> filter_server.log
         wait

         @ batch++
      end
      all_regions_done:

      # Need to have all backgrounds completed before continuing
      wait

      # signal to async_filter.csh to continue
      echo "Completed this assimilation at " `date` >> filter_server.log
      echo --------- >> filter_server.log
      rm -f go_assim_regions

   endif


#--------------------------------------------------------------------------------------- 

# No files found, wait and check again
   sleep 1
end

