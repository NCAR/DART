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

# Script to run the filter executable.
#
# Signals filter_server.csh to advance the model or assimilate the regions.
# Can be submitted to the batch queue by job.csh (the main control script) 
# or, better yet, run locally INSIDE job.csh (preferred, since most times 
# the processor running job.csh is sitting idle).
# 
# Executes filter on a distinct node from those used by filter_server.csh
# 
# Central directory for whole filter job (CENTRALDIR)
#    batch jobs submitted from there
#    semaphor files must (dis)appear there, 
#    I/O between filter and advance_model and assim_region goes through there.
#    Final output is put there

##=============================================================================
# The 'common' strategy between batch submission mechanisms is that the batch
# directives are embedded as shell comments. This provides us with the hope
# that we can use one script with several different batch mechanisms because
# the directives for one are interpreted as comments by the other.
# Both sets of directives are interpreted as comments when the script is 
# invoked interactively.
#
##=============================================================================
## This block of directives constitutes the preamble for the LSF queuing system 
## LSF is used on the IBM   Linux cluster 'lightning'
## LSF is used on the IMAGe Linux cluster 'coral'
## LSF is used on the IBM   AIX   cluster 'bluevista'
## The queues on lightning and bluevista are supposed to be similar.
##
## the normal way to submit to the queue is:    bsub < filter.csh
##
## an explanation of the most common directives follows:
##
## -J      job name    (master script job.csh presumes filter.xxxx.log)
## -o      output listing filename 
## -P      account number
## -n      number of tasks (processors)
## -x      exclusive use of node
## -R "span[ptile=(num procs you want on each node)]"
## -q queue    cheapest == [standby, economy, (regular,debug), premium] == $$$$
##=============================================================================
#BSUB -J filter
#BSUB -o filter.%J.log
#BSUB -q economy
#BSUB -n 1
#BSUB -P 868500xx

##=============================================================================
## This block of directives constitutes the preamble for the PBS queuing system 
## PBS is used on the CGD   Linux cluster 'bangkok'
## PBS is used on the CGD   Linux cluster 'calgary'
##
## the normal way to submit to the queue is:    qsub filter.csh
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
#PBS -N filter
#PBS -r n
#PBS -e filter.err
#PBS -o filter.log
#PBS -q medium
#PBS -l nodes=1:ppn=2

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
   set SCRATCHDIR = /ptmp/${user}/filter
   
else if ($?PBS_O_WORKDIR) then

   # PBS has a list of processors in a file whose name is (PBS_NODEFILE)

   set CENTRALDIR = $PBS_O_WORKDIR
   set JOBNAME = $PBS_JOBNAME
   set PROCNAMES = `cat $PBS_NODEFILE`
   set REMOTECMD = rsh
   set SCRATCHDIR = /scratch/local/${user}/filter

else if ($?OCOTILLO_NODEFILE) then

   set CENTRALDIR = `pwd`
   set JOBNAME = interactive_filter
   set PROCNAMES = "$host"
   set REMOTECMD = rsh
   set SCRATCHDIR = /var/tmp/${user}/filter
   
else

   set CENTRALDIR = `pwd`
   set JOBNAME = interactive_filter
   set PROCNAMES = "$host"
   set REMOTECMD = csh
   set SCRATCHDIR = /tmp/${user}/filter_server
   
endif

if ( ! $?REMOVE ) then
  set REMOVE = 'rm -rf'
endif
if ( ! $?COPY ) then
  set COPY = 'cp -p'
endif
if ( ! $?MOVE ) then
  set MOVE = 'mv -f'
endif

# capture the name of this script for use in messages.
set  myname = $0

# Set Variable for a 'master' logfile

set MASTERLOG = ${CENTRALDIR}/run_job.log

# This job's working directory. 
# We are ensuring it is empty when we start, and (later) GONE when we finish.

set tempdir = ${CENTRALDIR}/filter$$
${REMOVE} $tempdir
mkdir -p $tempdir
cd       $tempdir

echo "${myname} - running    on ${PROCNAMES}"  >> $MASTERLOG
echo "${myname} - CENTRALDIR is ${CENTRALDIR}" >> $MASTERLOG
echo "${myname} - tempdir    is $tempdir"      >> $MASTERLOG
echo "${myname} - running    on ${PROCNAMES}"
echo "${myname} - CENTRALDIR is ${CENTRALDIR}"
echo "${myname} - tempdir    is $tempdir"

#-----------------------
# Get filter input files

ln -s ${CENTRALDIR}/input.nml      .
ln -s ${CENTRALDIR}/caminput.nc    .
ln -s ${CENTRALDIR}/obs_seq.out    .
ln -s ${CENTRALDIR}/filter_ic_old* .
ln -s ${CENTRALDIR}/inflate_ic_old .
ln -s ${CENTRALDIR}/final_full.*   .

# Run the filter

${CENTRALDIR}/filter | tee -a ${CENTRALDIR}/run_filter.stout &

#-----------------
# Hang around forever for now and wait for go_xxx to appear here,
# or disappear from $CENTRALDIR

set again = true
# flags to keep track of when semaphor files disappear.
set go_advance_exist = false
set go_assim_exist = false
set nsec = 1


while($again == true)

   # When go_advance_model appears, copy it to the common directory to signal
   # advance_model to go.  If it already exists, do nothing.

   if(-e go_advance_model && ${go_advance_exist} == false) then

      # remove files from previous stage to conserve disc space
      if (-e filter_assim_region_out1) then
         ${REMOVE} ${CENTRALDIR}/filter_assim_region_out* 
         ${REMOVE}               filter_assim_region_out*
      endif

      # Write assim_model_state_ic size to central directory for filter_server.csh
      # This is ultimately used to ensure all files have been completely written
      # before proceeding.
      set list = `ls -lt assim_model_state_ic1`
      echo $list[5] >! ${CENTRALDIR}/assim_size

      # Move in 'creation' order.
      # Gives the O/S more time to finish writing the last files.
      # 'while' loop stops at first non-existent file.
      #set n = 1
      #while (-e assim_model_state_ic$n)
      #   ${MOVE} assim_model_state_ic$n ${CENTRALDIR}
      #   @ n++
      #end
      ${MOVE} assim_model_state_ic* ${CENTRALDIR}

      echo "filter- go_advance_model existence at "`date` >> $MASTERLOG
      echo "filter- go_advance_model existence at "`date`
      ls -l go* ${CENTRALDIR}/go* >> $MASTERLOG
      ls -l go* ${CENTRALDIR}/go*

      set go_advance_exist = true
      ${COPY} filter_control     ${CENTRALDIR} ;# used by filter_server.csh
      ${COPY} go_advance_model   ${CENTRALDIR} ;# semaphore file for filter_server.csh

   endif
      
   # When the 'central' ${CENTRALDIR}/go_advance_model file disappears
   # the model has advanced to the next observation time and we need to assimilate.
   # Remove files needed for previous stage,
   # Move assim updated files for all ensemble members to local filespace.
   # Then, removing the local go_advance_model file signals filter to proceed.  
   if(! -e ${CENTRALDIR}/go_advance_model && ${go_advance_exist} == true) then

      ${REMOVE} ${CENTRALDIR}/assim_model_state_ic*
      ${REMOVE}               assim_model_state_ic*

      ${MOVE}   ${CENTRALDIR}/assim_model_state_ud* .

      # signal 'filter.f90' to prep for assimilation.
      set go_advance_exist = false
      ${REMOVE} go_advance_model
   endif
      
   # When filter creates 'go_assim_regions', copy to central directory to 
   # signal assim_regions to go.  If it already exists, the assimilation is
   # running, we just need to wait for it to finish.
   if(-e go_assim_regions && ${go_assim_exist} == false) then
      ${REMOVE} ${CENTRALDIR}/assim_model_state_ud*
      ${REMOVE}               assim_model_state_ud*

      ${MOVE} filter_assim_region__in*  ${CENTRALDIR}
      ${MOVE} assim_region_control      ${CENTRALDIR}
      ${COPY} filter_assim_obs_seq      ${CENTRALDIR}
      ${COPY} go_assim_regions          ${CENTRALDIR}
      set go_assim_exist = true
   endif
      
   # When the go_assim_regions disappears from the central directory, 
   # remove it from here to signal filter to proceed.
   # Copy updated files for all ensemble members to filesystem on this node.
   if(! -e ${CENTRALDIR}/go_assim_regions && ${go_assim_exist} == true) then
      ${REMOVE} ${CENTRALDIR}/filter_assim_region__in*
      ${REMOVE}               filter_assim_region__in* 

      ${MOVE} ${CENTRALDIR}/filter_assim_region_out* .
      set go_assim_exist = false
      ${REMOVE} go_assim_regions 
   endif

   #----------------------------------------------------------------------------
   # When filter writes out go_end_filter... 
   #----------------------------------------------------------------------------
   if(-e go_end_filter) then
      #wait for new filter_ics to appear,
      set msec = 1
      set go = no
      while ($go == no)
         echo "${myname} - waiting for new filter_ics to appear."
         if ( -s filter_ic_new || -s filter_ic_new.0001 ) then
            set go = yes
         else
            sleep $msec
            if ($msec < 8) @ msec = 2 * $msec
         endif
      end

      echo "${myname} - filter_ic_new has size" >> $MASTERLOG
      echo "${myname} - filter_ic_new has size"
      ls -lt filter_ic_new*                     >> $MASTERLOG
      ls -lt filter_ic_new*

      ${MOVE} filter_ic_new*                    ${CENTRALDIR}
      ${MOVE} inflate_ic_new                    ${CENTRALDIR}

      echo "${myname} - copied filter_ic_new to ${CENTRALDIR}" >> $MASTERLOG
      echo "${myname} - copied filter_ic_new to ${CENTRALDIR}"
      ls -lt ${CENTRALDIR}/filter_ic*                          >> $MASTERLOG 
      ls -lt ${CENTRALDIR}/filter_ic*

      ${MOVE} Prior_Diag.nc Posterior_Diag.nc   ${CENTRALDIR}
      ${MOVE} obs_seq.final                     ${CENTRALDIR}

      # signal job.csh that filter is done with this obs_seq.out/day.
      ${COPY} go_end_filter ${CENTRALDIR}/go_end_filter_server
      ${COPY} go_end_filter ${CENTRALDIR}/go_end_filter

      echo "filter- finished assimilating"                 >> $MASTERLOG
      echo "filter- finished assimilating"
      echo "filter- moved go_end_filter to ${CENTRALDIR}"  >> $MASTERLOG
      echo "filter- moved go_end_filter to ${CENTRALDIR}"
      ls -lt ${CENTRALDIR}/go*                             >> $MASTERLOG 
      ls -lt ${CENTRALDIR}/go*

      # setting 'again' to false exits the loop, signals we are DONE.
      set again = false
      echo "filter- filter.csh terminating normally at "`date`  >> $MASTERLOG
      echo "filter- filter.csh terminating normally at "`date`
   else
      sleep $nsec
      if ($nsec < 8) @ nsec = 2 * $nsec
      # echo "filter- waiting for go_end_filter "`date`
   endif

end

#-----------------------------------------------
# Save the run-time logfile created within filter (name specified by namelist var:
# &utilities_nml    logfilename = 'dart_log.out'    [default] 
#-----------------------------------------------
echo "filter- working directory contents at end of filter.csh" >> $MASTERLOG
echo "filter- working directory contents at end of filter.csh"
pwd                                                            >> $MASTERLOG
pwd
ls -lrt                                                        >> $MASTERLOG
ls -lrt

${MOVE} dart_log.out ${CENTRALDIR}/filter_log.out
cd ${CENTRALDIR}

# To prevent losing restart files look for signal from job.csh that it's safe
# job.csh copies the filter_ic files to 'permanent' storage. 
# Until that is done, we don't want to remove them here. 

set again = true
while($again == true)
   if (-e rm_filter_temp) then
      rm -rf $tempdir rm_filter_temp
      set again = false
   else
      echo "filter- rm_filter_temp not found yet"  >> $MASTERLOG
      sleep 10
   endif
end
