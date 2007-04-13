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
#
#
#
#
#
#
# IF YOU READ ONLY ONE COMMENT -- READ THIS ONE.
# This script is designed to be as generic as possible - and still work.
# This script exploits the fact that different queueing systems set
# environment variables. If these variables exist, the values are passed
# to generic equivalents so that the vast majority of the script stays the
# same, i.e, we are not dragging around queuing-system-specific variables
# any longer than necessary.  TJH 13 Dec 2005

##=============================================================================
## This block of directives constitutes the preamble for the LSF queuing system 
## LSF is used on the IBM   Linux cluster 'lightning'
## LSF is used on the IMAGe Linux cluster 'coral'
## LSF is used on the IBM   cluster 'bluevista'
## The queues on lightning and bluevista are supposed to be similar.
##
## the normal way to submit to the queue is:    bsub < assim_filter.csh
##
## an explanation of the most common directives follows:
## -J Job name
## -o Output files
## -q queue    cheapest == [standby, economy, (regular,debug), premium] == $$$$
##=============================================================================
#BSUB -J assim_filter
#BSUB -o assim_filter.%J.o
#BSUB -q regular
#BSUB -n 8

##=============================================================================
## This block of directives constitutes the preamble for the PBS queuing system 
## PBS is used on the CGD   Linux cluster 'bangkok'
## PBS is used on the CGD   Linux cluster 'calgary'
##
## the normal way to submit to the queue is:    qsub assim_filter.csh
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
#PBS -N assim_filter
#PBS -r n
#PBS -e assim_filter.err
#PBS -o assim_filter.log
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
   set JOBNAME = assim_filter
   set PROCNAMES = `cat $OCOTILLO_NODEFILE`
   set REMOTECMD = rsh
   set SCRATCHDIR = /var/tmp/${user}

else                                    # interactive

   if ( ! $?host) then
      setenv host `uname -n`
   endif

   set CENTRALDIR = `pwd`
   set JOBNAME = interactive_assim_filter
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

#------------------------------------------------------------------------------
# This block is exactly the same as 'section 2' of filter_server.csh, with only
# one exception. There is no dedicated logfile for this block, so we just tack
# information onto the semaphor file 'batchflag'. If this dies prematurely,
# batchflag still exists and might contain something useful.
#------------------------------------------------------------------------------

set MASTERLOG = ${CENTRALDIR}/batchflag

      # First line of assim_region_control is the number of regions to be assimilated
      if ( -e assim_region_control) then
         set nregions = `head -1 assim_region_control`
      else
         set nregions = 1
      endif
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

      echo "$JOBNAME - Entering assim_region rerun block at "`date`
      echo "$JOBNAME - Entering assim_region rerun block at "`date` >> $MASTERLOG
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
            echo "$JOBNAME - nrerun, rerun = $nrerun $rerun on $PROCNAMES[$procnum]" >> $MASTERLOG
            echo "$JOBNAME - nrerun, rerun = $nrerun $rerun on $PROCNAMES[$procnum]"
         else
            echo "$JOBNAME - filter_assim_region_out$n is fine" >> $MASTERLOG
            echo "$JOBNAME - filter_assim_region_out$n is fine"
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

#------------------------------------------------------------------------------
# At this point, we differ from the async = 3 scenario of filter_server.csh.
# filter_server.csh communicates by removing a file 'go_assim_regions'.
#
# This script is invoked when (async = 2) and signals that it has completed
# by removing the semaphor file '${CENTRALDIR}/batchflag'
#------------------------------------------------------------------------------

echo "$JOBNAME - Completed this assimilation at " `date`
echo "$JOBNAME - ---------"
pwd
ls -lt filter_assim_region_out*

# signal to assim_tools_mod:filter_assim() (if async==2) to continue
${REMOVE} ${CENTRALDIR}/batchflag
