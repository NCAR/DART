#!/bin/csh

# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$

#----------------------------------------------------------------------
# Script to manage the execution of advance_model and assim_region.
# Submitted to batch queue by job.csh (the main control script).
#
# Calls advance_model.csh and assim_region.csh when signals are received
#     from filter, via job.csh.
#
# Executes on a set of compute processors; 
#    each processor advances an ensemble member and (possibly) a region.
#
# 11/11/04 Kevin Raeder
# modified by committee 3/31/05 Kevin Raeder, Alain Caya, typing by Tim Hoar
#----------------------------------------------------------------------
#
# It is wise to coordinate the number of nodes requested here with the number
# of regions defined in filter (num_domains).
# Each region and model advance runs on a single processor.
#
# The number of processors (NPROCS) is determined by one of two things.
# If this script is run interactively, NPROCS is unity.
# If this script is run in batch mode, NPROCS is user-defined by  
# the argument to 'bsub -n xxxx'
#
#### LSF options for BSUB
### -J      job name
### -o      output listing filename 
### -q      queue
### Queue charging    cheapest ..... most expensive   (at least on lightning)
### Queue name (standby, economy, [regular,debug], premium)
#
#BSUB -J filter_server
#BSUB -o filter_server.%J.log
#BSUB -q debug
#BSUB -n 1

if ($?LSB_MCPU_HOSTS) then 
   echo "LSB_MCPU_HOSTS is $LSB_MCPU_HOSTS"
endif
if ($?LS_SUBCWD ) then 
   echo "LS_SUBCWD     is $LS_SUBCWD"
endif
if ($?LSB_HOSTS) then 
   echo "LSB_HOSTS      is $LSB_HOSTS"
endif
if ($?LSB_EXECHOSTS) then 
   echo "LSB_EXECHOSTS  is $LSB_EXECHOSTS"
endif
if ($?LSB_JOBNAME) then 
   echo "LSB_JOBNAME    is $LSB_JOBNAME"
endif
if ($?LSB_JOBFILENAME) then 
   echo "LSB_JOBFILENAME is $LSB_JOBFILENAME"
endif

# Try to get a batch-system independent working directory name variable 
set SCRATCHDIR = /ptmp/$user/filter_server

### This job's working directory; must cd to it, or it will run in /home...
# this is the directory FROM WHICH you submitted the batch job.
if ( $?LS_SUBCWD ) then           ;# batch job
   # the central experiment directory is automatically defined
   set DARTWORKDIR = $LS_SUBCWD
   cd $DARTWORKDIR
else
   setenv DARTWORKDIR `pwd`
endif

# Determine number of processors
# each ensemble member or region is done entirely by one processor 
if ($?LSB_HOSTS) then                       ;# batch
   set NPROCS = `echo $LSB_HOSTS | wc -w`
   set PROCNAMES = $LSB_HOSTS
else                                        ;# interactive
   set NPROCS = 1
   set PROCNAMES = $host
endif

### Output to confirm job characteristics
echo " " >> $DARTWORKDIR/run_job.log
echo "server- Running filter_server on host "`hostname`  >> $DARTWORKDIR/run_job.log
echo "Initialized at "`date`                             >> $DARTWORKDIR/run_job.log
echo "DARTWORKDIR is "`pwd`                              >> $DARTWORKDIR/run_job.log
echo "This job has allocated $NPROCS procs"              >> $DARTWORKDIR/run_job.log
echo "they are: "                                        >> $DARTWORKDIR/run_job.log
echo $PROCNAMES                                          >> $DARTWORKDIR/run_job.log

# The script has two fundamental sections.
# The first part advances the ensemble members; 
# The second part advances the regions. 

set go_secs = 1
# Hang around forever for now and wait for go_{advance,assim_regions} file to appear
while (1 == 1)

   #--------------------------------------------------------------------------------
   # When filter wants the model to advance, it creates the file 'go_advance_model'
   # If go_advance_model exists, advance all the ensemble members.
   #--------------------------------------------------------------------------------
   if( -e go_advance_model ) then   ;# advance all ensemble members

      # First line of filter_control should have number of model states to be integrated
      if (-e filter_control) then
         set nensemble = `head -1 filter_control`
      else
         set nensemble = 1
      endif
      echo "server- advancing $nensemble members at " `date`  >> $DARTWORKDIR/run_job.log

      # figure # batches of model runs to do, from # ensemble members and # processors
      @ nbatch = $nensemble / $NPROCS
      if ($nensemble % $NPROCS != 0 ) @ nbatch++

      # Create a directory for each member to run in for namelists
      set element = 0
      set batch = 1
      while ($batch <= $nbatch)

         # don't wait for all the initial files, just the last one for this batch
         @ last_num = $batch * $NPROCS 
         set exist_secs = 1
         while (! -e assim_model_state_ic${last_num})
            echo "waiting for assim_model_state_ic${last_num} to exist" >> run_job.log
            sleep $exist_secs
            if ($exist_secs < 8) @ exist_secs = 2 * $exist_secs
         end

         # Ensure that the last batch does not start until its files are the right size.
         # read size of assim_model_state_ic from file written from perfect.csh or filter.csh
         # running on a compute node.  Compare size currently in DARTWORKDIR with the whole
         # file size.

         set size = `cat assim_size`
         set list = `ls -l assim_model_state_ic${last_num}`
         set size_secs = 1
         while ($list[5] != $size)
            echo "waiting for assim_model_state_ic${last_num} size $size list $list[5]"  \
                 >> $DARTWORKDIR/run_job.log 
            sleep $size_secs
            if ($size_secs < 8) @ size_secs = 2 * $size_secs
            set list = `ls -l assim_model_state_ic${last_num}`
         end

         # Advance the model for each ensemble member
         # advance_model has an additional optional arg for mpi "machines"
         foreach proc ( $PROCNAMES )
            echo " " >> $DARTWORKDIR/run_job.log
               @ element++
               if ($element > $nensemble) goto all_elements_done

               set ELEMENTDIR = $SCRATCHDIR/tmp$$_$element 

         echo "ssh $proc  csh $DARTWORKDIR/advance_model.csh $DARTWORKDIR $element $ELEMENTDIR" >> $DARTWORKDIR/run_job.log
               ssh $proc "csh $DARTWORKDIR/advance_model.csh $DARTWORKDIR $element $ELEMENTDIR" &

         end

         echo "server- waiting to finish batch $batch of $nbatch at "`date`  >> $DARTWORKDIR/run_job.log
         wait

         @ batch++
      end
      all_elements_done:

      # Make sure all processes complete before continuing
      wait
      \rm -rf ${HOME}/lnd*.rpointer   ;# leftover cam stuff

      echo "server- Completed this advance at " `date`  >> $DARTWORKDIR/run_job.log
      echo "---------"                                  >> $DARTWORKDIR/run_job.log
      pwd                                               >> $DARTWORKDIR/run_job.log
      ls -lt assim_model_state_ud*                      >> $DARTWORKDIR/run_job.log

      # wait till all ensemble members finish 
      # count all the full-size assim_model_state* files

      set go = true
      set n = 0
      while ($n < $nensemble && $go == true)
         @ n++
         if (-z assim_model_state_ud$n || ! -e assim_model_state_ud$n) set go = false
      end

      # finished with advance_model so remove the go_advance_model file
      # This signals 'filter' to proceed with the assimilation.

      if ($go == true) then
         rm -f go_advance_model
      else
         echo "MISSING assim_model_state_ud$n and aborting" >> $DARTWORKDIR/run_job.log
         echo "MISSING assim_model_state_ud$n and aborting" >> $DARTWORKDIR/run_job.log
         echo "MISSING assim_model_state_ud$n and aborting" >> $DARTWORKDIR/run_job.log
         exit
      endif

   endif

   #------------------------------------------------------------------------------------ 
   # Check to see if the go_advance_regions file exists
   #------------------------------------------------------------------------------------ 
   if( -e go_assim_regions ) then

      # First line of filter_control should have number of regions to be assimilated
      set nregions = `head -1 assim_region_control`
      echo "server- assimilating $nregions regions at " `date` >> $DARTWORKDIR/run_job.log

      # figure # batches of model runs to do, from # regions and # processors
      @ nbatch = $nregions / $NPROCS
      if ($nregions % $NPROCS != 0 ) @ nbatch++

      # Create a directory for each member to run in for namelists
      set element = 0
      set batch = 1
      while ($batch <= $nbatch)
         foreach proc ( $PROCNAMES )
            @ element++
            if ($element > $nregions) goto all_regions_done

            set ELEMENTDIR = $SCRATCHDIR/tmp$$_$element

      echo "ssh $proc  csh $DARTWORKDIR/assim_region.csh $DARTWORKDIR $element $ELEMENTDIR" >> $DARTWORKDIR/run_job.log
            ssh $proc "csh $DARTWORKDIR/assim_region.csh $DARTWORKDIR $element $ELEMENTDIR" &

         end
         echo "server- waiting to finish batch $batch of $nbatch regions" `date`   >> $DARTWORKDIR/run_job.log
         wait

         @ batch++
      end
      all_regions_done:

      echo "TJHdebug - before wait " `date`  >> $DARTWORKDIR/run_job.log
      # Make sure all (background) processes complete before continuing
      wait

      # Removing 'go_assim_regions' is the signal to filter.csh to continue
      echo "server- Completed this assimilation at " `date`  >> $DARTWORKDIR/run_job.log
      echo "---------"                                       >> $DARTWORKDIR/run_job.log
      set go = true
      set n = 0
      while ($n < $nregions && $go == true)
         @ n++
         if (-z filter_assim_region_out$n || ! -e filter_assim_region_out$n) set go = false
      end
      if ($go == true) then
         rm -f go_assim_regions
      else
         echo "MISSING filter_assim_region_out$n and aborting" >> $DARTWORKDIR/run_job.log
         echo "MISSING filter_assim_region_out$n and aborting" >> $DARTWORKDIR/run_job.log
         echo "MISSING filter_assim_region_out$n and aborting" >> $DARTWORKDIR/run_job.log
         exit
      endif

   endif

   # When filter is done with the obs_sequence file, 'go_end_filter' is created.
   # When filter.csh detects 'go_end_filter', 'go_end_filter_server' is created.
   # When 'go_end_filter_server' exists, we can stop.

   if(-e go_end_filter_server) then
      echo "server- terminating normally at " `date`  >> $DARTWORKDIR/run_job.log
      \rm -f go_end_filter*
      exit
   else
      # No files found, wait and check again
      sleep $go_secs
      if ($go_secs < 8) @ go_secs = 2 * $go_secs
      # date
   endif

#--------------------------------------------------------------------------------------- 
end
