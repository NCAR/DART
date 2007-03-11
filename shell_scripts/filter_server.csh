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
##=============================================================================
## This block of directives constitutes the preamble for the LSF queuing system 
## LSF is used on the IBM   Linux cluster 'lightning'
## LSF is used on the IMAGe Linux cluster 'coral'
## LSF is used on the IBM   'bluevista'
## The queues on lightning and bluevista are supposed to be similar.
##
## the normal way to submit to the queue is:    bsub < filter_server.csh
##
## an explanation of the most common directives follows:
## -J Job name (master script job.csh presumes filter_server.xxxx.log)
## -o STDOUT filename
## -e STDERR filename
## -P      account
## -q queue    cheapest == [standby, economy, (regular,debug), premium] == $$$$
## -n number of processors  (really)
##=============================================================================
#BSUB -J filter_server
#BSUB -o filter_server.%J.log
#BSUB -e filter_server.%J.err
#BSUB -P 86850054
#BSUB -q regular
#BSUB -n 8
#
#
##=============================================================================
## This block of directives constitutes the preamble for the PBS queuing system 
## PBS is used on the CGD   Linux cluster 'bangkok'
## PBS is used on the CGD   Linux cluster 'calgary'
##
## the normal way to submit to the queue is:    qsub filter_server.csh
##
## an explanation of the most common directives follows:
## -N     Job name
## -r n   Declare job non-rerunable
## -e <arg>  filename for standard error 
## -o <arg>  filename for standard out 
## -q <arg>   Queue name (small, medium, long, verylong)
## -l nodes=xx:ppn=2   requests BOTH processors on the node. On both bangkok 
##                     and calgary, there is no way to 'share' the processors 
##                     on the node with another job, so you might as well use 
##                     them both.  (ppn == Processors Per Node)
##=============================================================================
#PBS -N filter_server
#PBS -r n
#PBS -e filter_server.err
#PBS -o filter_server.log
#PBS -q medium
#PBS -l nodes=4:ppn=2

# A common strategy for the beginning is to check for the existence of
# some variables that get set by the different queuing mechanisms.
# This way, we know which queuing mechanism we are working with,
# and can set 'queue-independent' variables for use for the remainder 
# of the script.

if ($?LS_SUBCWD) then                   # LSF

   # LSF has a list of processors already in a variable (LSB_HOSTS)

   set CENTRALDIR = $LS_SUBCWD
   set JOBNAME = $LSB_JOBNAME
   set PROCNAMES = ($LSB_HOSTS)
   set REMOTECMD = ssh
   set SCRATCHDIR = /ptmp/${user}/filter_server
   
else if ($?PBS_O_WORKDIR) then          # PBS

   # PBS has a list of processors in a file whose name is (PBS_NODEFILE)

   set CENTRALDIR = $PBS_O_WORKDIR
   set JOBNAME = $PBS_JOBNAME
   set PROCNAMES = `cat $PBS_NODEFILE`
   set REMOTECMD = rsh
   set SCRATCHDIR = /scratch/local/${user}/filter_server

else if ($?OCOTILLO_NODEFILE) then          # ocotillo

   # ocotillo is a 'special case'. It is the only cluster I know of with
   # no queueing system.  You must generate a list of processors in a 
   # file whose name is in $OCOTILLO_NODEFILE.  For example ... 
   # setenv OCOTILLO_NODEFILE  my_favorite_processors
   # echo "node1"  > $OCOTILLO_NODEFILE
   # echo "node5" >> $OCOTILLO_NODEFILE
   # echo "node7" >> $OCOTILLO_NODEFILE
   # echo "node3" >> $OCOTILLO_NODEFILE

   set CENTRALDIR = `pwd`
   set JOBNAME = interactive_filter_server
   set PROCNAMES = `cat $OCOTILLO_NODEFILE`
   set REMOTECMD = rsh
   set SCRATCHDIR = /var/tmp/${user}/filter_server
   
else                                    # interactive

   set CENTRALDIR = `pwd`
   set JOBNAME = interactive_filter_server
   set PROCNAMES = "$host $host $host $host"
   set REMOTECMD = csh
   set SCRATCHDIR = /tmp/${user}/filter_server
   
endif

if ( ! $?REMOVE ) then
  set REMOVE = 'rm -rf'
endif

# The scratch directory will be used as the run-time directory for both the
# model advances and the regional assimilations. The SCRATCHDIR performance is
# best if ALL of the necessary directories are created ONCE before the experiment
# is run. (stooopid GPFS).
# Theoretically, job.csh could set this up (by reading input.nml) ...
#              set ELEMENTDIR = $SCRATCHDIR/member_$element 
#              set ELEMENTDIR = $SCRATCHDIR/region_$element 
# If you want the SCRATCHDIR to be a sub-directory of CENTRALDIR, 
# uncomment the following line. 

# set SCRATCHDIR = ${CENTRALDIR}/filter_server$$

# Everything below this should not need to be changed.
# There are many useful comments, however ...

# This job's working directory; must cd to it, or it may run in /home...

cd $CENTRALDIR

set NPROCS = `echo $PROCNAMES | wc -w`

# Set Variable for a 'master' logfile
set MASTERLOG = ${CENTRALDIR}/run_job.log

# Make it here ... once ... lest everyone fight over making it
# when they each try to make their own sub-directory.
if ( ! -d ${SCRATCHDIR} ) then
   mkdir -p ${SCRATCHDIR} || exit 3
endif

### Output to confirm job characteristics
echo " "                                             >! $MASTERLOG
echo "Running $JOBNAME on host "`hostname`           >> $MASTERLOG
echo "Initialized at "`date`                         >> $MASTERLOG
echo "CENTRALDIR is "`pwd`                           >> $MASTERLOG
echo "SCRATCHDIR is ${SCRATCHDIR}"                   >> $MASTERLOG
echo "This job has allocated $NPROCS processors."    >> $MASTERLOG
echo "they are: "                                    >> $MASTERLOG
echo $PROCNAMES                                      >> $MASTERLOG

echo " "
echo "Running $JOBNAME on host "`hostname`
echo "Initialized at "`date`
echo "CENTRALDIR is "`pwd`
echo "SCRATCHDIR is ${SCRATCHDIR}"
echo "This job has allocated $NPROCS processors."
echo "they are: "
echo $PROCNAMES

# The script has three fundamental sections.
# The first  part advances the ensemble members; 
# the second part advances the regions; 
# the third  part cleans up and terminates. 

set go_secs = 1
# Hang around forever for now and wait for go_{advance,assim_regions} file to appear
while (1 == 1)

   #--------------------------------------------------------------------------------
   # Section 1. Advance the models (ensemble members).
   #--------------------------------------------------------------------------------
   # When filter wants the model to advance, it creates the file 'go_advance_model'
   # If go_advance_model exists, advance all the ensemble members.
   #--------------------------------------------------------------------------------
   if( -e go_advance_model ) then   # advance all ensemble members

      # First line of filter_control should have number of model states to be integrated
      if (-e filter_control) then
         set nensemble = `head -1 filter_control`
      else
         set nensemble = 1
      endif
      echo "$JOBNAME - advancing $nensemble members at " `date`  >> $MASTERLOG
      echo "$JOBNAME - advancing $nensemble members at " `date`

      # figure # batches of model runs to do, from # ensemble members and # processors
      @ nbatch = $nensemble / $NPROCS
      if ($nensemble % $NPROCS != 0 ) @ nbatch++

      # Create a directory for each member to run in for namelists
      set element = 0
      set batch = 1
      while ($batch <= $nbatch)

         # Advance the model for each ensemble member
         # advance_model has an additional optional arg for mpi "machines"
         foreach proc ( $PROCNAMES )
            @ element++
            if ($element > $nensemble) goto some_elements_done

            set ELEMENTDIR = $SCRATCHDIR/member_$element 

            if ($REMOTECMD == 'csh') then  # interactive
               set proc = " "
            endif

            echo "$REMOTECMD $proc  $CENTRALDIR/advance_model.csh $CENTRALDIR $element $ELEMENTDIR"
            echo "$REMOTECMD $proc  $CENTRALDIR/advance_model.csh $CENTRALDIR $element $ELEMENTDIR" >> $MASTERLOG
                  $REMOTECMD $proc  $CENTRALDIR/advance_model.csh $CENTRALDIR $element $ELEMENTDIR &
         end

         echo "$JOBNAME - waiting to finish ensemble batch $batch of $nbatch at "`date`  >> $MASTERLOG
         echo "$JOBNAME - waiting to finish ensemble batch $batch of $nbatch at "`date`
         wait

         @ batch++
      end
      some_elements_done:

      # Make sure all processes complete before continuing
      wait

      # Count all the ?full-size? assim_model_state* files to see if there
      # are any model advances that need to be rerun. If the file is zero
      # length, or not present -- we assume the model advance failed for some reason.
      # We will try again ... once ... and give up if that does not work.

      echo "$JOBNAME - Entering advance_model rerun block at "`date`
      echo "$JOBNAME - Entering advance_model rerun block at "`date` >> $MASTERLOG
      echo " " >> $MASTERLOG

      set n = 0
      set nrerun = 0
      while ($n < $nensemble)
         @ n++

         ls -l assim_model_state_ud$n

         if (-z assim_model_state_ud$n || ! -e assim_model_state_ud$n) then
            @ nrerun++
            @ procnum = $n % $NPROCS

            echo "$JOBNAME - advance_model - failed procnum is $procnum"
            echo "$JOBNAME - advance_model - failed procnum is $procnum" >> $MASTERLOG

            if ($procnum == 0) set procnum = $NPROCS

            if ($nrerun == 1) then
               set rerun = ($n)
               set badprocs = ($PROCNAMES[$procnum])
            else
               set rerun = ($rerun $n)
               set badprocs = ($badprocs $PROCNAMES[$procnum])
            endif
            echo "$JOBNAME - nrerun, rerun = $nrerun $rerun on $PROCNAMES[$procnum]" >> $MASTERLOG
            echo "$JOBNAME - nrerun, rerun = $nrerun $rerun on $PROCNAMES[$procnum]"
         else
            echo "$JOBNAME - assim_model_state_ud$n is fine" >> $MASTERLOG
            echo "$JOBNAME - assim_model_state_ud$n is fine"
         endif
      end

      # If there were some bad model advances, we must find the processors
      # that worked, and log ones that didn't.  We will resubmit to just the
      # list of good processors.
      
      if ($nrerun > 0) then
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
         echo ' '                                         >> $MASTERLOG
         echo "Ensemble; nrerun rerun = $nrerun $rerun"   >> $MASTERLOG
         echo "       good processors = $goodprocs"       >> $MASTERLOG
         echo ' '                                         >> $MASTERLOG
      endif

      # The number of batches and nodes limit how many bad nodes can be handled; 
      # 4 bad processors (=2 bad nodes) x 4 batches = 16 members to redo on 16 procs
      # Give up if the number to redo is greater than the number of processors.

      if ($nrerun > 0) then
         set element = 0
         foreach proc ( $goodprocs )
            @ element++
            if ($element > $nrerun) goto all_elements_done

            set ELEMENTDIR = $SCRATCHDIR/member_$rerun[$element] 

            if ($REMOTECMD == 'csh') then # interactive
               set proc = " "
            endif

            echo "$REMOTECMD $proc  $CENTRALDIR/advance_model.csh $CENTRALDIR $rerun[$element] $ELEMENTDIR" >> $MASTERLOG
            echo "$REMOTECMD $proc  $CENTRALDIR/advance_model.csh $CENTRALDIR $rerun[$element] $ELEMENTDIR"
                  $REMOTECMD $proc  $CENTRALDIR/advance_model.csh $CENTRALDIR $rerun[$element] $ELEMENTDIR &

         end
         all_elements_done:

         echo "$JOBNAME - waiting to finish ensemble rerun "`date`  >> $MASTERLOG
         echo "$JOBNAME - waiting to finish ensemble rerun "`date`
         wait
      endif

      # OK, we tried. If the ones we needed to rerun are not done now ... give up.

      set n = 0
      while ($n < $nensemble)
         @ n++
         if (-z assim_model_state_ud$n || ! -e assim_model_state_ud$n) then
            echo "MISSING assim_model_state_ud$n and aborting" >> $MASTERLOG
            echo "MISSING assim_model_state_ud$n and aborting" >> $MASTERLOG
            echo "MISSING assim_model_state_ud$n and aborting"
            echo "MISSING assim_model_state_ud$n and aborting"
            exit 1
         endif
      end

      # Finished with advance_model (hopefully) so remove the go_advance_model file.
      # This signals 'filter' to proceed with the assimilation.

      ${REMOVE} go_advance_model

      echo "$JOBNAME - Completed this advance at " `date`  >> $MASTERLOG
      echo "$JOBNAME - Completed this advance at " `date`
      echo "---------"                                  >> $MASTERLOG
      echo "---------"
      pwd                                               >> $MASTERLOG
      pwd
      ls -lt assim_model_state_ud*                      >> $MASTERLOG 
      ls -lt assim_model_state_ud*

   endif

   #------------------------------------------------------------------------------------ 
   # Section 2. Assimilate region-by-region.
   # TJH should just use the goodprocs ... not all the procs ... if it was bad for
   # the model advance, it is probably still bad now ... 
   #--------------------------------------------------------------------------------
   # Check to see if the go_assim_regions file exists
   #------------------------------------------------------------------------------------ 
   if( -e go_assim_regions ) then

      # First line of assim_region_control is the number of regions to be assimilated
      set nregions = `head -1 assim_region_control`
      echo "$JOBNAME - assimilating $nregions regions at " `date` >> $MASTERLOG
      echo "$JOBNAME - assimilating $nregions regions at " `date`

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

            if ($REMOTECMD == 'csh') then  # interactive
               set proc = " "
            endif

      echo "$REMOTECMD $proc  $CENTRALDIR/assim_region.csh $CENTRALDIR $element $ELEMENTDIR" >> $MASTERLOG
      echo "$REMOTECMD $proc  $CENTRALDIR/assim_region.csh $CENTRALDIR $element $ELEMENTDIR"
            $REMOTECMD $proc  $CENTRALDIR/assim_region.csh $CENTRALDIR $element $ELEMENTDIR &

         end

         echo "$JOBNAME - waiting to finish regional batch $batch of $nbatch " `date`   >> $MASTERLOG
         echo "$JOBNAME - waiting to finish regional batch $batch of $nbatch " `date`
         wait

         @ batch++
      end
      some_regions_done:

      # Make sure all (background) processes complete before continuing
      wait

      # Count all the ?full-size? filter_assim_region_out* files to see if there
      # are any regional assimilations that need to be rerun. If the file is zero
      # length, or not present -- we assume an assimilation failed for some reason.
      # We will try again ... once ... and give up if that does not work.

      echo "Entering assim_region rerun block at "`date`
      echo "Entering assim_region rerun block at "`date` >> $MASTERLOG
      echo " " >> $MASTERLOG

      set n = 0
      set nrerun = 0
      while ($n < $nregions)
         @ n++

         ls -l filter_assim_region_out$n

         if (-z filter_assim_region_out$n || ! -e filter_assim_region_out$n) then
            @ nrerun++
            @ procnum = $n % $NPROCS

            echo "$JOBNAME - assim_region - failed procnum is $procnum"
            echo "$JOBNAME - assim_region - failed procnum is $procnum" >> $MASTERLOG

            if ($procnum == 0) set procnum = $NPROCS

            if ($nrerun == 1) then
               set rerun = ($n)
               set badprocs = ($PROCNAMES[$procnum])
            else
               set rerun = ($rerun $n)
               set badprocs = ($badprocs $PROCNAMES[$procnum])
            endif
            echo "nrerun, rerun = $nrerun $rerun on $PROCNAMES[$procnum]" >> $MASTERLOG
            echo "nrerun, rerun = $nrerun $rerun on $PROCNAMES[$procnum]"
         else
            echo "filter_assim_region_out$n is fine" >> $MASTERLOG
            echo "filter_assim_region_out$n is fine"
         endif
      end

      # If there were some bad regional assimilations, we must find the processors
      # that worked, and log ones that didn't.  We will resubmit to just the
      # list of good processors.

      if ($nrerun > 0) then

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
         echo ' '                                      >> $MASTERLOG
         echo "Regions; nrerun rerun = $nrerun $rerun" >> $MASTERLOG
         echo "      good processors = $goodprocs"     >> $MASTERLOG
         echo ' '                                      >> $MASTERLOG
      endif

      # The number of batches and nodes limit how many bad nodes can be handled; 
      # 4 bad processors (=2 bad nodes) x 4 batches = 16 members to redo on 16 procs
      # Give up if the number to redo is greater than the number of processors.

      if ($nrerun > 0) then
         set element = 0
         foreach proc ( $goodprocs )
            @ element++
            if ($element > $nrerun) goto all_regions_done

            set ELEMENTDIR = $SCRATCHDIR/region_$rerun[$element] 

            if ($REMOTECMD == 'csh') then # interactive
               set proc = " "
            endif

            echo "$REMOTECMD $proc  $CENTRALDIR/assim_region.csh $CENTRALDIR $rerun[$element] $ELEMENTDIR"
            echo "$REMOTECMD $proc  $CENTRALDIR/assim_region.csh $CENTRALDIR $rerun[$element] $ELEMENTDIR" >> $MASTERLOG
                  $REMOTECMD $proc  $CENTRALDIR/assim_region.csh $CENTRALDIR $rerun[$element] $ELEMENTDIR  &

         end
         all_regions_done:

         echo "$JOBNAME - waiting to finish regions rerun "`date`  >> $MASTERLOG
         echo "$JOBNAME - waiting to finish regions rerun "`date`
         wait
      endif

      # OK, we tried. If the ones we needed to rerun are not done now ... give up.

      set n = 0
      while ($n < $nregions)
         @ n++
         if (-z filter_assim_region_out$n || ! -e filter_assim_region_out$n) then
            echo "MISSING filter_assim_region_out$n and aborting" >> $MASTERLOG
            echo "MISSING filter_assim_region_out$n and aborting" >> $MASTERLOG
            echo "MISSING filter_assim_region_out$n and aborting"
            echo "MISSING filter_assim_region_out$n and aborting"
            exit 2
         endif
      end

      # Finished with assim_region (hopefully) so remove the go_assim_regions file.
      # This signals 'filter' to proceed with the assimilation.

      ${REMOVE} go_assim_regions

      # Removing 'go_assim_regions' is the signal to filter.csh to continue
      echo "$JOBNAME - Completed this assimilation at " `date` >> $MASTERLOG
      echo "$JOBNAME - Completed this assimilation at " `date`
      echo "---------"                                      >> $MASTERLOG
      echo "---------"
      pwd                                                   >> $MASTERLOG
      pwd
      ls -lt filter_assim_region_out*                       >> $MASTERLOG
      ls -lt filter_assim_region_out*

   endif

   #--------------------------------------------------------------------------------
   # Section 3. Clean up and terminate.
   #--------------------------------------------------------------------------------
   # When filter is done with the obs_sequence file, 'go_end_filter' is created.
   # CAM filter.csh detects 'go_end_filter', 'go_end_filter_server' is created.
   # CAM When 'go_end_filter_server' exists, we can stop.
   #--------------------------------------------------------------------------------

   if(-e go_end_filter ) then
      echo "$JOBNAME - terminating normally at " `date`  >> $MASTERLOG
      echo "$JOBNAME - terminating normally at " `date`
      ${REMOVE} go_end_filter
      exit 0
   else
      # No files found, wait and check again
      sleep $go_secs
      if ($go_secs < 8) @ go_secs = 2 * $go_secs
   endif

end
