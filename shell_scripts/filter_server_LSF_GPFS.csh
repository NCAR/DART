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
### -J      job name (master script job.csh presumes filter_server.xxxx.log)
### -o      output listing filename 
### -q      queue
### -P      account
### Queue charging    cheapest ..... most expensive   (at least on lightning)
### Queue name (standby, economy, [regular,debug], premium)
#
#BSUB -J filter_server
#BSUB -o filter_server.%J.log
#BSUB -P 86850054
#BSUB -q economy
#BSUB -n 16

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

# Determine number of processors -- one of three ways.
# 1) Batch jobs set a variable LSB_HOSTS
# 2) Interactive jobs can have a NPROCS environment variable defined.
# 3) Interactive jobs default to 1 (one).
#
# list of hosts/machines is in $PROCNAMES
# the quoting is VERY IMPORTANT for PROCNAMES
#
# each ensemble member or region is done entirely by one processor 
if ($?LSB_HOSTS) then                       ;# batch
   set NPROCS = `echo $LSB_HOSTS | wc -w`
   set PROCNAMES = "$LSB_HOSTS"
else if ($?NPROCS) then
   set PROCNAMES = $host
   set iproc = 2
   while($iproc <= $NPROCS)
      set PROCNAMES = "$PROCNAMES $host"
      @ iproc ++
   end
else                                        ;# interactive
   set NPROCS = 1
   set PROCNAMES = $host
endif

### This job's working directory; must cd to it, or it will run in /home...
# this is the directory FROM WHICH you submitted the batch job.
if ( $?LS_SUBCWD ) then           ;# batch job
   # the central experiment directory is automatically defined
   set CENTRALDIR = $LS_SUBCWD
   cd $CENTRALDIR
else
   setenv CENTRALDIR `pwd`
endif

# Set Variable for a 'master' logfile
set MASTERLOG = ${CENTRALDIR}/run_job.log

# Try to get a batch-system independent working directory name variable 
set SCRATCHDIR = ${CENTRALDIR}/filter_server$$

### Output to confirm job characteristics
echo " "                                                 >> $MASTERLOG
echo "server- Running filter_server on host "`hostname`  >> $MASTERLOG
echo "Initialized at "`date`                             >> $MASTERLOG
echo "CENTRALDIR is "`pwd`                               >> $MASTERLOG
echo "This job has allocated $NPROCS procs"              >> $MASTERLOG
echo "they are: "                                        >> $MASTERLOG
echo $PROCNAMES                                          >> $MASTERLOG

echo " "
echo "server- Running filter_server on host "`hostname`
echo "Initialized at "`date`
echo "CENTRALDIR is "`pwd`
echo "This job has allocated $NPROCS procs"
echo "they are: "
echo $PROCNAMES

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
      echo "server- advancing $nensemble members at " `date`  >> $MASTERLOG
      echo "server- advancing $nensemble members at " `date`

      # figure # batches of model runs to do, from # ensemble members and # processors
      @ nbatch = $nensemble / $NPROCS
      if ($nensemble % $NPROCS != 0 ) @ nbatch++

      # Create a directory for each member to run in for namelists
      set element = 0
      set batch = 1
      while ($batch <= $nbatch)

        ## don't wait for all the initial files, just the last one for this batch
        #@ last_num = $batch * $NPROCS 
        #set exist_secs = 1
        #while (! -e assim_model_state_ic${last_num})
        #   echo "waiting for assim_model_state_ic${last_num} to exist" >> run_job.log
        #   sleep $exist_secs
        #   if ($exist_secs < 8) @ exist_secs = 2 * $exist_secs
        #end

        ## Ensure that the last batch does not start until its files are the right size.
        ## read size of assim_model_state_ic from file written from perfect.csh or filter.csh
        ## running on a compute node.  Compare size currently in CENTRALDIR with the whole
        ## file size.

        #set size = `cat assim_size`
        #set list = `ls -l assim_model_state_ic${last_num}`
        #set size_secs = 1
        #while ($list[5] != $size)
        #   echo "waiting for assim_model_state_ic${last_num} size $size list $list[5]"  \
        #        >> $MASTERLOG 
        #   sleep $size_secs
        #   if ($size_secs < 8) @ size_secs = 2 * $size_secs
        #   set list = `ls -l assim_model_state_ic${last_num}`
        #end

         # Advance the model for each ensemble member
         # advance_model has an additional optional arg for mpi "machines"
         foreach proc ( $PROCNAMES )
            echo " " >> $MASTERLOG
               @ element++
               if ($element > $nensemble) goto all_elements_done

               set ELEMENTDIR = $SCRATCHDIR/member_$element 

         echo "ssh $proc  csh $CENTRALDIR/advance_model.csh $CENTRALDIR $element $ELEMENTDIR" >> $MASTERLOG
         echo "ssh $proc  csh $CENTRALDIR/advance_model.csh $CENTRALDIR $element $ELEMENTDIR"
               ssh $proc "csh $CENTRALDIR/advance_model.csh $CENTRALDIR $element $ELEMENTDIR" &

         end

         echo "server- waiting to finish batch $batch of $nbatch at "`date`  >> $MASTERLOG
         echo "server- waiting to finish batch $batch of $nbatch at "`date`
         wait

         @ batch++
      end
      all_elements_done:

      # Make sure all processes complete before continuing
      wait
      \rm -rf ${HOME}/lnd*.rpointer   ;# leftover cam stuff

      echo "server- Completed this advance at " `date`  >> $MASTERLOG
      echo "server- Completed this advance at " `date`
      echo "---------"                                  >> $MASTERLOG
      echo "---------"
      pwd                                               >> $MASTERLOG
      pwd
      ls -lt assim_model_state_ud*                      >> $MASTERLOG 
      ls -lt assim_model_state_ud*

      # wait till all ensemble members finish 
      # count all the ?full-size? assim_model_state* files

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
         echo "MISSING assim_model_state_ud$n and aborting" >> $MASTERLOG
         echo "MISSING assim_model_state_ud$n and aborting" >> $MASTERLOG
         echo "MISSING assim_model_state_ud$n and aborting"
         echo "MISSING assim_model_state_ud$n and aborting"
         exit
      endif

   endif

   #------------------------------------------------------------------------------------ 
   # Check to see if the go_advance_regions file exists
   #------------------------------------------------------------------------------------ 
   if( -e go_assim_regions ) then

      # First line of filter_control should have number of regions to be assimilated
      set nregions = `head -1 assim_region_control`
      echo "server- assimilating $nregions regions at " `date` >> $MASTERLOG
      echo "server- assimilating $nregions regions at " `date`

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

            set ELEMENTDIR = $SCRATCHDIR/region_$element 

      echo "ssh $proc  csh $CENTRALDIR/assim_region.csh $CENTRALDIR $element $ELEMENTDIR" >> $MASTERLOG
      echo "ssh $proc  csh $CENTRALDIR/assim_region.csh $CENTRALDIR $element $ELEMENTDIR"
            ssh $proc "csh $CENTRALDIR/assim_region.csh $CENTRALDIR $element $ELEMENTDIR" &

         end

         echo "server- waiting to finish regional batch $batch of $nbatch " `date`   >> $MASTERLOG
         echo "server- waiting to finish regional batch $batch of $nbatch " `date`
         wait

         @ batch++
      end
      all_regions_done:

      # Make sure all (background) processes complete before continuing
      wait

      # Removing 'go_assim_regions' is the signal to filter.csh to continue
      echo "server- Completed this assimilation at " `date`  >> $MASTERLOG
      echo "server- Completed this assimilation at " `date`
      echo "---------"                                       >> $MASTERLOG
      echo "---------"
      set go = true
      set n = 0
      while ($n < $nregions && $go == true)
         @ n++
         if (-z filter_assim_region_out$n || ! -e filter_assim_region_out$n) set go = false
      end

      if ($go == true) then
         rm -f go_assim_regions
      else
         echo "MISSING filter_assim_region_out$n and aborting" >> $MASTERLOG
         echo "MISSING filter_assim_region_out$n and aborting" >> $MASTERLOG
         echo "MISSING filter_assim_region_out$n and aborting"
         echo "MISSING filter_assim_region_out$n and aborting"
         exit
      endif

   endif

   # When filter is done with the obs_sequence file, 'go_end_filter' is created.
   # When filter.csh detects 'go_end_filter', 'go_end_filter_server' is created.
   # When 'go_end_filter_server' exists, we can stop.

   #if(-e go_end_filter ) then
   if(-e go_end_filter_server ) then
      echo "server- terminating normally at " `date`  >> $MASTERLOG
      echo "server- terminating normally at " `date`
      \rm -f go_end_filter_server
      rmdir -v $SCRATCHDIR
      exit
   else
      # No files found, wait and check again
      sleep $go_secs
      if ($go_secs < 8) @ go_secs = 2 * $go_secs
      # date
   endif

#--------------------------------------------------------------------------------------- 
end
