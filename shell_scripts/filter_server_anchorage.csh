#!/bin/tcsh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$
#
#---------------------------
# Script to manage the execution of CAM and assim_region.
# Submitted to batch queue by job.csh (the main control script).
# Calls advance_model.csh and assim_region.csh when signals are received
#     from filter, via job.csh.
# Executes on a compute node, which also hosts one of the ensemble member advances
#     and one (or more) of the regions being assimilated.
# last modified 11/11/04 Kevin Raeder
#---------------------------
#
# It is wise to coordinate the number of nodes requested here with the number
# of regions defined in filter (num_domains), keeping in mind the number of
# processors/ node.  Each region can runs on a single processor, while each
# CAM is set to run on a single node, using all processors on that node.
#
# Note that the PBS -l nodes only defines the number of nodes when the
# script is submitted as a batch job by qsub.
# Otherwise, set NPROCS (not NNODES) further down in the script.

#   Multi-processor jobs must be submitted as batch, under PBS
### Job name
#PBS -N filter_server
### Declare job non-rerunable
#PBS -r n
### Output files
#PBS -e filter_server.err
#PBS -o filter_server.log
### Queue name (small, medium, long, verylong)
#PBS -q medium
#PBS -l nodes=3:ppn=1


### This job's working directory; must cd to it, or it will run in /home...
if ($?PBS_O_WORKDIR) then
   # it's a batch job and the central experiment directory is automatically defined
   cd $PBS_O_WORKDIR
else
   setenv PBS_O_WORKDIR `pwd`
endif

### Output to confirm job characteristics
echo ' ' >> $PBS_O_WORKDIR/run_job.log
if ($?PBS_JOBNAME) then
   echo server- Running $PBS_JOBNAME on host `hostname`  >> $PBS_O_WORKDIR/run_job.log
else
   echo server- "Running on host "`hostname`  >> $PBS_O_WORKDIR/run_job.log
endif
echo Initialized at `date`  >> $PBS_O_WORKDIR/run_job.log
echo Directory is `pwd`  >> $PBS_O_WORKDIR/run_job.log
echo This job runs on the following processors:  >> $PBS_O_WORKDIR/run_job.log

### Define number of processors; # of lines in PBS_NODEFILE
# Machines with > 1 processor/node will have corresponding repetitions of
# the node name in PBS_NODEFILE
if ($?PBS_NODEFILE) then
   # batch job submitte by qsub
   setenv NPROCS `wc -l < $PBS_NODEFILE`

   set prev_proc = none
   set NPROCS_PER_NODE = 1
   foreach proc ( `cat $PBS_NODEFILE` )
      if ($proc == $prev_proc) then
         # This fails on the first one of a multi-proc machine, which is fine 
         # because the first is alread counted.  The next test will make sure it
         # proceeds to the second one.
         @ NPROCS_PER_NODE++
      else if ($prev_proc != none ) then
         goto end_count_procs
      endif
      set prev_proc = $proc
   end
   end_count_procs:

   @ NNODES = $NPROCS / $NPROCS_PER_NODE

   cat $PBS_NODEFILE  >> $PBS_O_WORKDIR/run_job.log
   echo PBS_NODEFILE is $PBS_NODEFILE >> $PBS_O_WORKDIR/run_job.log
else
   # not a batch job
   setenv PBS_NODEFILE nodefile
   rm -f $PBS_NODEFILE
   set NPROCS = 10
   set inode = 1
   while($inode <= $NPROCS)
      echo node$inode >> $PBS_NODEFILE
      echo node$inode 
      @ inode ++
   end
endif
echo server- This job has allocated $NPROCS procs  >> $PBS_O_WORKDIR/run_job.log

# Hang around forever for now and wait for go_{advance,assim_regions} file to appear
while(1 == 1)

#-----------------------------------------------------------------------------------
   # Check to see if the go_advance_model file exists
   ls go_advance_model >& .option3_garb
   if($status == 0) then

      # First line of filter_control should have number of model states to be integrated
      set nensemble = `head -1 filter_control`
      echo "server- advancing $nensemble members at " `date`  >> $PBS_O_WORKDIR/run_job.log

      # figure # batches of CAM runs to do, from # ensemble members and # processors
      @ nbatch = $nensemble / $NNODES
      if ($nensemble % $NNODES != 0 ) @ nbatch++

      # Create a directory for each member to run in for namelists
      set element = 0
      set batch = 1
      # mechanism for submitting CAM only once to a node, not once to each processor
      set prev_proc = none
      while($batch <= $nbatch)
#        don't wait for all the initial files, just the last one for this batch
         @ last_num = $batch * $NNODES + 1
         set nsec = 1
         while (! -e assim_model_state_ic${last_num} && $batch != $nbatch)
            sleep $nsec
            if ($nsec < 8) @ nsec = 2 * $nsec
         end

         foreach proc ( `cat $PBS_NODEFILE` )
            echo ' ' >> $PBS_O_WORKDIR/run_job.log
            if ($proc != $prev_proc || $NPROCS == 1) then
               @ element++
               if ($element > $nensemble) goto all_elements_done

               if ($NPROCS_PER_NODE == 1) then
                  rsh $proc "csh $PBS_O_WORKDIR/advance_model.csh $PBS_O_WORKDIR $element /scratch/local/tmp$$$user$element " &
               else
                  echo $proc >! $PBS_O_WORKDIR/machine$element
                  set n = 2
                  while ($n <= $NPROCS_PER_NODE)
                     echo $proc >> $PBS_O_WORKDIR/machine$element
                     @ n++
                  end
                  rsh $proc "csh $PBS_O_WORKDIR/advance_model.csh $PBS_O_WORKDIR $element /scratch/local/tmp$$$user$element $PBS_O_WORKDIR/machine$element " &

               endif 
               echo rsh $proc "csh $PBS_O_WORKDIR/advance_model.csh $PBS_O_WORKDIR $element /scratch/local/tmp$$$user$element machine$element "   >> $PBS_O_WORKDIR/run_job.log

            endif
            set prev_proc = $proc
         end
         echo server- waiting to finish batch $batch of $nbatch at `date`  >> $PBS_O_WORKDIR/run_job.log
         wait

         @ batch++
      end
      all_elements_done:

      # Make sure all processes complete before continuing
      wait

      rm $PBS_O_WORKDIR/machine*

      # finished with advance_model so remove the go_advance_model file
      echo "server- Completed this advance at " `date`  >> $PBS_O_WORKDIR/run_job.log
      echo ---------  >> $PBS_O_WORKDIR/run_job.log
      ls -lt assim_model_state_ud*  >> $PBS_O_WORKDIR/run_job.log
      set go = true
      set n = 0
      while ($n < $nensemble && $go == true)
         @ n++
         if (-z assim_model_state_ud$n || ! -e assim_model_state_ud$n) set go = false
      end
      if ($go == true) then
         rm -f go_advance_model
      else
         echo MISSING assim_model_state_ud$n and aborting >> $PBS_O_WORKDIR/run_job.log
         exit
      endif

   endif

#--------------------------------------------------------------------------------------- 
   # Check to see if the go_advance_regions file exists
   ls go_assim_regions >& .option3_garb
   if($status == 0) then

      # First line of filter_control should have number of regions to be assimilated
      set nregions = `head -1 assim_region_control`
      echo "server- assimilating $nregions regions at " `date` >> $PBS_O_WORKDIR/run_job.log

      # figure # batches of CAM runs to do, from # regions and # processors
      @ nbatch = $nregions / $NPROCS
      if ($nregions % $NPROCS != 0 ) @ nbatch++

      # Create a directory for each member to run in for namelists
      set element = 0
      set batch = 1
      while($batch <= $nbatch)
         foreach proc ( `cat $PBS_NODEFILE` )
            @ element++
            if ($element > $nregions) goto all_regions_done

            rsh $proc "csh $PBS_O_WORKDIR/assim_region.csh $PBS_O_WORKDIR $element /scratch/local/tmp$$$user$element " &

         end
         echo "server- waiting to finish batch $batch of $nbatch regions" `date`   >> $PBS_O_WORKDIR/run_job.log
         wait

         @ batch++
      end
      all_regions_done:

      # Make sure all (background) processes complete before continuing
      wait

      # signal to async_filter.csh to continue
      rm -f go_assim_regions
      echo "server- Completed this assimilation at " `date`  >> $PBS_O_WORKDIR/run_job.log
      echo ---------  >> $PBS_O_WORKDIR/run_job.log

   endif

   # If go_end_filter_server exists then stop this process
   # No wait is necessary because the script is inactive while advance_model or assim_region
   # is active, and a wait is undesirable because when they're done we want to exit immediately 
   # when go_end_filter appears
   if(-e go_end_filter) then
      echo "server- terminating normally at " `date`  >> $PBS_O_WORKDIR/run_job.log
      exit
   endif

#-----------------------------------------------------------------------------------
end
