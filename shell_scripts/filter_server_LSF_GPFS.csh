#!/bin/csh

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

#----------------------------------------------------------------------
#
# Script to manage the execution of advance_model and assim_region.
# Submitted to batch queue by job.csh (the main control script).
#
# Calls advance_model.csh and assim_region.csh when signals are received
#     from filter, via job.csh.
#
# Executes on a set of compute processors; 
#    each processor advances an ensemble member and (possibly) a region.
#
#----------------------------------------------------------------------
#
# It is wise to coordinate the number of nodes requested here with the number
# of regions defined in filter (num_domains).
# Each region and model advance runs on a single processor.
# For example, 80 ensemble members runs nicely as 4 batches of 20.
# Therefore, 20 regions would be nice because 20 processors could be used.
# One processor for each region ... ergo BSUB -n 20    ... capice? 
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
### -n      number of processors  (really)
#
#BSUB -J filter_server
#BSUB -o filter_server.%J.log
#BSUB -e filter_server.%J.err
#BSUB -P 86850054
#BSUB -q regular
#BSUB -n 20

# The scratch directory will be used as the run-time directory for both the
# model advances and the regional assimilations. The SCRATCHDIR performance is
# best if ALL of the necessary directories are created ONCE before the experiment
# is run. (stooopid GPFS).
# Theoretically, job.csh could set this up (by reading input.nml) ...
#              set ELEMENTDIR = $SCRATCHDIR/member_$element 
#              set ELEMENTDIR = $SCRATCHDIR/region_$element 
# If you want the SCRATCHDIR to be a sub-directory of CENTRALDIR, there is a 
# line to uncomment about 60 lines below this. 
set SCRATCHDIR = /ptmp/${user}/filter_server

# Everything below this should not need to be changed.
# There are many useful comments, however ...

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
#
# each ensemble member or region is done entirely by one processor 
#
if ($?LSB_HOSTS) then                       ;# batch
   set NPROCS = `echo $LSB_HOSTS | wc -w`
   set PROCNAMES = ($LSB_HOSTS)
   echo "Using this block"
else if ($?NPROCS) then
   set PROCNAMES = $host
   set iproc = 2
   while($iproc <= $NPROCS)
      set PROCNAMES = ($PROCNAMES $host)
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
# set SCRATCHDIR = ${CENTRALDIR}/filter_server$$

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
               @ element++
               if ($element > $nensemble) goto some_elements_done

               set ELEMENTDIR = $SCRATCHDIR/member_$element 

         echo "ssh $proc  csh $CENTRALDIR/advance_model.csh $CENTRALDIR $element $ELEMENTDIR"
         echo "ssh $proc  csh $CENTRALDIR/advance_model.csh $CENTRALDIR $element $ELEMENTDIR" >> $MASTERLOG
               ssh $proc "csh $CENTRALDIR/advance_model.csh $CENTRALDIR $element $ELEMENTDIR" &

         end

         echo "server- waiting to finish batch $batch of $nbatch at "`date`  >> $MASTERLOG
         echo "server- waiting to finish batch $batch of $nbatch at "`date`
         wait

         @ batch++
      end
      some_elements_done:

      # Make sure all processes complete before continuing
      wait

      # wait till all ensemble members finish 
      # count all the ?full-size? assim_model_state* files

      echo " " >> $MASTERLOG
      set n = 0
      set nrerun = 0
      while ($n < $nensemble)
         @ n++
         if (-z assim_model_state_ud$n || ! -e assim_model_state_ud$n) then
            @ nrerun++
            @ procnum = $n % $NPROCS
            if ($procnum == 0) set procnum = $NPROCS

            if ($nrerun == 1) then
               set rerun = ($n)
               set badprocs = ($PROCNAMES[$procnum])
            else
               set rerun = ($rerun $n)
               set badprocs = ($badprocs $PROCNAMES[$procnum])
            endif
            echo "nrerun, rerun = $nrerun $rerun on $PROCNAMES[$procnum]" >> $MASTERLOG
         else
            echo "assim_model_state_ud$n is fine " >> $MASTERLOG
         endif
      end

      # sift out good processors
      #
      if ($nrerun > 0) then
      #  (echo "$badprocs" | mail -s "node(s) died" sghosh@ucar.edu)  &
         (echo "$badprocs" | mail -s "node(s) died" $user@ucar.edu)  &

         set ngood = 0
         set goodprocs = ' '
         foreach proc ($PROCNAMES)
            set is_good = 'true'
            foreach bad ($badprocs)
               if ($proc == $bad) set is_good = 'false'
            end
            if ($is_good == 'true') then
               @ ngood++
               if ($ngood == 1 ) then
                  set goodprocs = ($proc)
               else
                  set goodprocs = ($goodprocs $proc)
               endif
            endif
         end
         echo ' '       >> $MASTERLOG
         echo "Ensemble; nrerun rerun = $nrerun $rerun"      >> $MASTERLOG
         echo "          good processors = $goodprocs"       >> $MASTERLOG
         echo ' '       >> $MASTERLOG
      endif

      # finished with advance_model (hopefully) so remove the go_advance_model file
      # This signals 'filter' to proceed with the assimilation.

      if ($nrerun == 0) then
         rm -f go_advance_model
      else
         # number of batches and nodes limit how many bad nodes can be handled; 
         # 4 bad processors (=2 bad nodes) x 4 batches = 16 members to redo on 16 procs
         # Give up if the number to redo is greater than the number of processors.
         set element = 0
         foreach proc ( $goodprocs )
            @ element++
            if ($element > $nrerun) goto all_elements_done

            set ELEMENTDIR = $SCRATCHDIR/member_$rerun[$element] 

            echo "ssh $proc  csh $CENTRALDIR/advance_model.csh $CENTRALDIR $rerun[$element] $ELEMENTDIR" >> $MASTERLOG
            echo "ssh $proc  csh $CENTRALDIR/advance_model.csh $CENTRALDIR $rerun[$element] $ELEMENTDIR"
                  ssh $proc "csh $CENTRALDIR/advance_model.csh $CENTRALDIR $rerun[$element] $ELEMENTDIR" &

         end
         all_elements_done:

         echo "server- waiting to finish ensemble rerun "`date`  >> $MASTERLOG
         echo "server- waiting to finish ensemble rerun "`date`
         wait
      endif

      set go = true
      set n = 0
      while ($n < $nensemble)
         @ n++
         if (-z assim_model_state_ud$n || ! -e assim_model_state_ud$n) then
            set go = false
            echo "MISSING assim_model_state_ud$n and aborting" >> $MASTERLOG
            echo "MISSING assim_model_state_ud$n and aborting" >> $MASTERLOG
            echo "MISSING assim_model_state_ud$n and aborting"
            echo "MISSING assim_model_state_ud$n and aborting"
            exit
         endif
      end

      if ($go == true) rm -f go_advance_model

      echo "server- Completed this advance at " `date`  >> $MASTERLOG
      echo "server- Completed this advance at " `date`
      echo "---------"                                  >> $MASTERLOG
      echo "---------"
      pwd                                               >> $MASTERLOG
      pwd
      ls -lt assim_model_state_ud*                      >> $MASTERLOG 
      ls -lt assim_model_state_ud*

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
            if ($element > $nregions) goto some_regions_done

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
      some_regions_done:

      # Make sure all (background) processes complete before continuing
      wait
    echo ' ' >> $MASTERLOG

      echo " " >> $MASTERLOG
      set n = 0
      set nrerun = 0
      while ($n < $NPROCS)
         @ n++
         if (-z filter_assim_region_out$n || ! -e filter_assim_region_out$n) then
            @ nrerun++
            @ procnum = $n % $NPROCS
            if ($procnum == 0) set procnum = $NPROCS

            if ($nrerun == 1) then
               set rerun = ($n)
               set badprocs = ($PROCNAMES[$procnum])
            else
               set rerun = ($rerun $n)
               set badprocs = ($badprocs $PROCNAMES[$procnum])
            endif
            echo "nrerun, rerun = $nrerun $rerun on $PROCNAMES[$procnum]" >> $MASTERLOG
         else
            echo "filter_assim_region_out$n is fine " >> $MASTERLOG
         endif
      end

      # sift out good processors
      #
      if ($nrerun > 0) then
         set ngood = 0
         foreach proc ($PROCNAMES)
            set is_good = 'true'
            foreach bad ($badprocs)
               if ($proc == $bad) set is_good = 'false'
            end
            if ($is_good == 'true') then
               @ ngood++
               if ($ngood == 1 ) then
                  set goodprocs = ($proc)
               else
                  set goodprocs = ($goodprocs $proc)
               endif
            endif
         end
         echo ' ' >> $MASTERLOG
         echo "Regions; nrerun rerun = $nrerun $rerun" >> $MASTERLOG
         echo "         good processorss = $goodprocs" >> $MASTERLOG
         echo ' ' >> $MASTERLOG
      endif

      if ($nrerun == 0) then
         rm -f go_assim_regions
      else
         set element = 0
         foreach proc ( $goodprocs )
            echo " " >> $MASTERLOG
            @ element++
            if ($element > $nrerun) goto all_regions_done

            set ELEMENTDIR = $SCRATCHDIR/region_$rerun[$element] 

            echo "ssh $proc  csh $CENTRALDIR/assim_region.csh $CENTRALDIR $rerun[$element] $ELEMENTDIR"
            echo "ssh $proc  csh $CENTRALDIR/assim_region.csh $CENTRALDIR $rerun[$element] $ELEMENTDIR" >> $MASTERLOG
                  ssh $proc "csh $CENTRALDIR/assim_region.csh $CENTRALDIR $rerun[$element] $ELEMENTDIR" &

         end
         echo "server- waiting to finish regions rerun "`date`  >> $MASTERLOG
         echo "server- waiting to finish regions rerun "`date`

         all_regions_done:
         wait
      endif

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

      # Removing 'go_assim_regions' is the signal to filter.csh to continue
      echo "server- Completed this assimilation at " `date`  >> $MASTERLOG
      echo "server- Completed this assimilation at " `date`
      echo "---------"                                       >> $MASTERLOG
      echo "---------"

   endif

   # When filter is done with the obs_sequence file, 'go_end_filter' is created.
   # When filter.csh detects 'go_end_filter', 'go_end_filter_server' is created.
   # When 'go_end_filter_server' exists, we can stop.

   if(-e go_end_filter_server ) then
      echo "server- terminating normally at " `date`  >> $MASTERLOG
      echo "server- terminating normally at " `date`
      \rm -f go_end_filter_server
      exit
   else
      # No files found, wait and check again
      sleep $go_secs
      if ($go_secs < 8) @ go_secs = 2 * $go_secs
   endif

#--------------------------------------------------------------------------------------- 
end
