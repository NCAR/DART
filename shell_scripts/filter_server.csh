#!/bin/tcsh

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
#PBS -q medium
#PBS -l nodes=4


### This job's working directory; must cd to it, or it will run in /home...
if ($?PBS_O_WORKDIR) then
   cd $PBS_O_WORKDIR
else
   setenv PBS_O_WORKDIR `pwd`
endif

# Clear out the log file
rm -f filter_server.log

### Output to confirm job characteristics
if ($?PBS_JOBNAME) then
   echo Running $PBS_JOBNAME on host `hostname` >> filter_server.log
else
   echo "Running on host "`hostname` >> filter_server.log
endif
echo Initialized at `date` >> filter_server.log
echo Directory is `pwd` >> filter_server.log
echo This job runs on the following processors: >> filter_server.log

### Define number of processors; # of lines in PBS_NODEFILE
if ($?PBS_NODEFILE) then
   setenv NPROCS `wc -l < $PBS_NODEFILE`
   cat $PBS_NODEFILE >> filter_server.log
else
   setenv PBS_NODEFILE nodefile
   rm -f $PBS_NODEFILE
   set NPROCS = 10
   set inode = 1
   while($inode <= $NPROCS)
      echo node$inode >> $PBS_NODEFILE
      echo node$inode >> filter_server.log
      @ inode ++
   end
endif
echo This job has allocated $NPROCS nodes >> filter_server.log

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
      echo "advancing $nensmbl members at " `date` >> filter_server.log

      # figure # batches of CAM runs to do, from # ensemble members and # processors
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

            rsh $node "csh $PBS_O_WORKDIR/advance_model.csh $PBS_O_WORKDIR $element /scratch/local/tmp$user$$$element " &

         end
         echo waiting to finish batch $batch at `date` >> filter_server.log
         wait

         @ batch++
      end
      all_elements_done:

      # finished with advance_model so remove the go_advance_model file
      echo "Completed this advance at " `date` >> filter_server.log
      echo --------- >> filter_server.log
      rm -f go_advance_model

   endif

#--------------------------------------------------------------------------------------- 
   # Check to see if the go_advance_regions file exists
   ls go_assim_regions > .option3_garb
   if($status == 0) then

      # First line of filter_control should have number of regions to be assimilated
      set nregions = `head -1 assim_region_control`
      echo "assimilating $nregions regions at " `date`>> filter_server.log

      # figure # batches of CAM runs to do, from # regions and # processors
      @ nbatch = $nregions / $NPROCS
      if ($nregions % $NPROCS != 0 ) @ nbatch++
      echo $nbatch batches will be executed >> filter_server.log

      # Create a directory for each member to run in for namelists
      set element = 0
      set batch = 1
      while($batch <= $nbatch)
         foreach node ( `cat $PBS_NODEFILE` )
            @ element++
            if ($element > $nregions) goto all_regions_done

            rsh $node "csh $PBS_O_WORKDIR/assim_region.csh $PBS_O_WORKDIR $element /scratch/local/tmp$user$$$element " &

         end
         echo "waiting to finish batch $batch" `date`  >> filter_server.log
         wait

         @ batch++
      end
      all_regions_done:

      # signal to async_filter.csh to continue
      echo "Completed this assimilation at " `date` >> filter_server.log
      echo --------- >> filter_server.log
      rm -f go_assim_regions

   endif


#--------------------------------------------------------------------------------------- 

# No files found, wait and check again
   #sleep 1
end
