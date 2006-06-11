#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section, 
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$
#
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
#

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
#BSUB -n 9

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
#PBS -l nodes=5:ppn=2

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
   set SCRATCH_PREFIX = /ptmp/${user}

else if ($?PBS_O_WORKDIR) then          # PBS

   # PBS has a list of processors in a file whose name is (PBS_NODEFILE)

   set CENTRALDIR = $PBS_O_WORKDIR
   set JOBNAME = $PBS_JOBNAME
   set PROCNAMES = `cat $PBS_NODEFILE`
   set REMOTECMD = rsh
   set SCRATCH_PREFIX = /scratch/local/${user}

else if ($?OCOTILLO_NODEFILE) then      # ocotillo

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
   set SCRATCH_PREFIX = /var/tmp/${user}

else                                    # interactive

   set CENTRALDIR = `pwd`
   set JOBNAME = interactive_assim_filter
   set PROCNAMES = "$host $host $host $host"
   set REMOTECMD = csh
   set SCRATCH_PREFIX = /tmp/${user}

endif

if ( ! $?REMOVE ) then
  set REMOVE = 'rm -rf'
endif

# This job's working directory; must cd to it, or it may run in /home...

cd $CENTRALDIR

set NPROCS = `echo $PROCNAMES | wc -w`

# Output to confirm job characteristics

echo "Running $JOBNAME on host "`hostname`
echo "Time is "`date`
echo "Directory is $CENTRALDIR"
echo "This job has allocated $NPROCS processors."
echo "This job runs on the following nodes:"
echo $PROCNAMES
echo " "

# First line of assim_region_control should have number of regions to be assimilated.
# Given the number of ensemble members and processors, determine the number
# of 'batches' necessary using primitive shell commands.

set nregions = `head -1 assim_region_control`
@ nbatch = $nregions / $NPROCS
if ($nregions % $NPROCS != 0 ) @ nbatch++
echo "$nbatch batches of regions will be executed."

# Send jobs to nodes. Each job is done in the 'background'.
set element = 0
set batch = 1
while($batch <= $nbatch)
   foreach proc ( $PROCNAMES )
      @ element++
      if ($element > $nregions) goto all_elements_done

      if ($REMOTECMD == 'csh') then
         set proc = " "
      endif

      echo "$REMOTECMD $proc $CENTRALDIR/assim_region.csh $CENTRALDIR $element ${SCRATCH_PREFIX}/tmp$element" 
            $REMOTECMD $proc $CENTRALDIR/assim_region.csh $CENTRALDIR $element ${SCRATCH_PREFIX}/tmp$element &

   end

   # One way to monitor progress.  batchflag has other info to start,
   # so this echo can be removed and scripts will still work.
   echo "waiting to finish assim batch $batch of $nbatch"
   echo "waiting to finish assim batch $batch of $nbatch" >> $CENTRALDIR/batchflag

   wait

   @ batch++
end
all_elements_done:

# Wait for all *background* processes to finish up ... just in case.
wait

# signal to filter_assim to continue
${REMOVE} $CENTRALDIR/batchflag
