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

# DESCRIPTION:
#
# The strategy here is that this job has requested a bunch of processors
# and has many model advances (ensemble members). We are going to assign 
# a model advance to a processor. Since there are more model advances 
# (more ensemble members) than processors, we will have to advance a bunch 
# of ensemble members, wait for that to finish, then do another batch ... 
# until we have advanced all the ensemble members. We are running each 
# ensemble member using just one processor, even though some of the models 
# might be able to take advantage of multiple processors.
#
# This script also has logic in it to check to see if the model advances
# completed successfully. After all the members have been advanced once, 
# a check is made to see if the desired output files exist. If they do not,
# ONE attempt is made to gather up all the members that need to be rerun and 
# give it one more try -- this time without using the processors that failed
# in the first place. If that fails, we give up.
#
# 'Section 1' of filter_server.csh is functionally identical to this script.
# For that reason, I have kept the same syntax and formatting.
#
# This script usually gets invoked when filter_nml namelist variables
# async == 2
# adv_ens_command = "./advance_ens.csh"
#
# If you want to take advantage of the queueing systems, 
#
# adv_ens_command = "bsub < advance_ens.csh". 
#
#
# The flow outline: 
#    'filter' is running in the CENTRALDIR
#    'filter' creates a semaphor file called batchflag 
#    'filter' fires off a shell command (advance_ens.csh) to advance the models. 
#    'filter' then waits patiently for 'batchflag' to disappear.
#     advance_ens.csh (this script) advances everyone.
#     when advance_ens.csh completes, remove batchflag ... filter proceeds,
#     presumably to the assimilation. 
#
##=============================================================================
# The 'common' strategy between batch submission mechanisms is that the batch
# directives are embedded as shell comments. This provides us with the hope
# that we can use one script with several different batch mechanisms because
# the directives for one are interpreted as comments by the other.
# Both sets of directives are interpreted as comments when the script is 
# invoked interactively.
#
# The number of processors (NPROCS) is determined by one of two things.
# If this script is run in batch mode, NPROCS is defined by the batch 
# submission directives.  Run interactively, NPROCS will be 2 ... and
# multiple processes will get run on the same processor/host.
#
##=============================================================================
## This block of directives constitutes the preamble for the LSF queuing system 
## LSF is used on the IBM   Linux cluster 'lightning'
## LSF is used on the IMAGe Linux cluster 'coral'
## LSF is used on the IBM   cluster 'bluevista'
## The queues on lightning and bluevista are supposed to be similar.
##
## the normal way to submit to the queue is:    bsub < advance_ens_LSF_PBS.csh
##
## an explanation of the most common directives follows:
## -J Job name
## -o Output files
## -q queue    cheapest == [standby, economy, (regular,debug), premium] == $$$$
##=============================================================================
#BSUB -J advance_ens
#BSUB -o adv_ens.%J.log
#BSUB -q regular
#BSUB -n 8

##=============================================================================
## This block of directives constitutes the preamble for the PBS queuing system 
## PBS is used on the CGD   Linux cluster 'bangkok'
## PBS is used on the CGD   Linux cluster 'calgary'
##
## the normal way to submit to the queue is:    qsub advance_ens_LSF_PBS.csh
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
#PBS -N advance_ens
#PBS -r n
#PBS -e advance_ens.err
#PBS -o advance_ens.log
#PBS -q medium
#PBS -l nodes=4:ppn=2

# A common strategy for the beginning is to check for the existence of
# some variables that get set by the different queuing mechanisms.
# This way, we know which queuing mechanism we are working with,
# and can set 'queue-independent' variables for use for the remainder 
# of the script.

if ($?LS_SUBCWD) then

   # LSF has a list of processors already in a variable (LSB_HOSTS)

   set CENTRALDIR = $LS_SUBCWD
   set JOBNAME = $LSB_JOBNAME
   set PROCNAMES = ($LSB_HOSTS)
   set REMOTECMD = ssh
   set SCRATCHDIR = /ptmp/${user}

else if ($?PBS_O_WORKDIR) then

   # PBS has a list of processors in a file whose name is (PBS_NODEFILE)

   set CENTRALDIR = $PBS_O_WORKDIR
   set JOBNAME = $PBS_JOBNAME
   set PROCNAMES = `cat $PBS_NODEFILE`
   set REMOTECMD = rsh
   set SCRATCHDIR = /scratch/local/${user}

else if ($?OCOTILLO_NODEFILE) then

   # ocotillo is a 'special case'. It is the only cluster I know of with
   # no queueing system.  You must generate a list of processors in a 
   # file whose name is in $OCOTILLO_NODEFILE.  For example ... 
   # setenv OCOTILLO_NODEFILE  my_favorite_processors
   # echo "node1"  > $OCOTILLO_NODEFILE
   # echo "node5" >> $OCOTILLO_NODEFILE
   # echo "node7" >> $OCOTILLO_NODEFILE
   # echo "node3" >> $OCOTILLO_NODEFILE

   set CENTRALDIR = `pwd`
   set JOBNAME = advance_ens
   set PROCNAMES = `cat $OCOTILLO_NODEFILE`
   set REMOTECMD = rsh
   set SCRATCHDIR = /var/tmp/${user}

else                                    # interactive

   if ( ! $?host) then
      setenv host `uname -n`
   endif

   set CENTRALDIR = `pwd`
   set JOBNAME = interactive_advance_ens
   set PROCNAMES = "$host $host $host $host"
   set REMOTECMD = csh
   set SCRATCHDIR = `pwd`

endif

if ( ! $?REMOVE ) then
  set REMOVE = 'rm -rf'
endif

# This job's working directory; must cd to it, or it may run in /home...

cd $CENTRALDIR

set NPROCS = `echo $PROCNAMES | wc -w`


# Output to confirm job characteristics
echo " "
echo "Running $JOBNAME on host "`hostname`
echo "Time is "`date`
echo "(central) directory is $CENTRALDIR"
echo "This job has allocated $NPROCS processors."
echo "They are:"
echo $PROCNAMES
echo " "

# Set a variable for the semaphor file that indicates the ensemble advances
# have completed. This MUST be named $CENTRALDIR/batchflag ... so don't change it. 
# We can log some run-time stuff to the file with no fear. If this does not
# terminate normally, the batchflag file might contain useful information.
# If this script terminates normally, the batchflag is removed -- to indicate
# to the rest of dart that it is done.
#------------------------------------------------------------------------------
# This block is exactly the same as 'section 1' of filter_server.csh, with only
# one exception. There is no dedicated logfile for this block, so we just tack
# information onto the semaphor file 'batchflag'. If this dies prematurely,
# batchflag still exists and might contain something useful.
#------------------------------------------------------------------------------

set MASTERLOG = ${CENTRALDIR}/batchflag

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

#------------------------------------------------------------------------------
# At this point, we differ from the async = 3 scenario of filter_server.csh.
# filter_server.csh communicates by removing a file 'go_advance_model'.
#
# This script is invoked when (async = 2) and signals that it has completed
# by removing the semaphor file '${CENTRALDIR}/batchflag'
#------------------------------------------------------------------------------

echo "$JOBNAME - Completed this advance at " `date`
echo "$JOBNAME - ---------"
pwd
ls -lt assim_model_state_ud*

# signal ensemble_manager_mod:Aadvance_state() (if async=2) to continue
${REMOVE} ${CENTRALDIR}/batchflag
