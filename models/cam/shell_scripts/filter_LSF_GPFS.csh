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

# Script to run the filter executable and signals filter_server to advance the model
#    or assimilate the regions
# Submitted to batch queue by job.csh (the main control script).
# Executes filter on a compute node, distinct from those used by filter_server.csh

#### LSF options for BSUB
### -J      job name    (master script job.csh presumes filter.xxxx.log)
### -o      output listing filename 
### -P      account number
### -q      queue
### -n      number of tasks (processors)
### -x      exclusive use of node
### -R "span[ptile=(num procs you want on each node)]"
#
#BSUB -J filter
#BSUB -o filter.%J.log
#BSUB -q economy 
#BSUB -n 1
#xxxx -x
#xxxx -R "span[ptile=1]"

# If they exist, log the values of the batch environment vars.
if ($?LSB_JOBNAME)     then
   echo "LSB_JOBNAME     is $LSB_JOBNAME"
endif
if ($?LSB_JOBFILENAME) then
   echo "LSB_JOBFILENAME is $LSB_JOBFILENAME"
endif
if ($?LSB_MCPU_HOSTS)  then
   echo "LSB_MCPU_HOSTS  is $LSB_MCPU_HOSTS"
endif
if ($?LS_SUBCWD )      then
   echo "LS_SUBCWD       is $LS_SUBCWD"
endif
if ($?LSB_HOSTS)       then
   echo "LSB_HOSTS       is $LSB_HOSTS"
endif
if ($?LSB_EXECHOSTS)   then
   echo "LSB_EXECHOSTS   is $LSB_EXECHOSTS"
endif

# Determine number of processors -- one of three ways.
# 1) Batch jobs set a variable LSB_HOSTS
# 2) Interactive jobs can have a NPROCS environment variable defined.
# 3) Interactive jobs default to 1 (one).
#
# list of hosts/machines is in $PROCNAMES
# the quoting is VERY IMPORTANT for PROCNAMES

if ($?LSB_HOSTS) then
   set NPROCS = `echo $LSB_HOSTS | wc -w`
   set PROCNAMES = "$LSB_HOSTS"
else if ($?NPROCS) then
   set PROCNAMES = $host
   set iproc = 2
   while($iproc <= $NPROCS)
      set PROCNAMES = "$PROCNAMES $host"
      @ iproc ++
   end
else
   set NPROCS = 1
   set PROCNAMES = $host
endif

# Central directory for whole filter job (CENTRALDIR)
#    batch jobs submitted from there
#    semaphor files must (dis)appear there, 
#    I/O between filter and advance_model and assim_region goes through there.
#    Final output is put there
if ( $?LS_SUBCWD ) then
   set CENTRALDIR = $LS_SUBCWD
else
   setenv CENTRALDIR `pwd`
endif

# Set Variable for a 'master' logfile
set MASTERLOG = ${CENTRALDIR}/run_job.log

# This job's working directory. 
# We are ensuring it is empty when we start, and (later) GONE when we finish.
#set tempdir = /ptmp/$user/tmp$$
set tempdir = ${CENTRALDIR}/filter$$
\rm -rf  $tempdir
mkdir -p $tempdir
cd       $tempdir

echo "filter- running tasks on ${PROCNAMES}"    >> $MASTERLOG
echo "filter- CENTRALDIR is ${CENTRALDIR}"      >> $MASTERLOG
echo "filter- tempdir is $tempdir"              >> $MASTERLOG
echo "filter- running tasks on ${PROCNAMES}"
echo "filter- CENTRALDIR is ${CENTRALDIR}"
echo "filter- tempdir is $tempdir"

#-----------------------
# Get filter input files
ln -s ${CENTRALDIR}/input.nml      .
ln -s ${CENTRALDIR}/caminput.nc    .
ln -s ${CENTRALDIR}/obs_seq.out    .
ln -s ${CENTRALDIR}/filter_ic_old* .
ln -s ${CENTRALDIR}/assim_ic_old   .

#-----------------------
# Run the filter
${CENTRALDIR}/filter >> ${CENTRALDIR}/run_filter.stout &

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

      # remove files needed for previous stage, to conserve disc space
      if (-e filter_assim_region_out1) then
         \rm -f ${CENTRALDIR}/filter_assim_region_out* 
         \rm -f               filter_assim_region_out*
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
      #   mv assim_model_state_ic$n ${CENTRALDIR}
      #   @ n++
      #end
      mv -v assim_model_state_ic* ${CENTRALDIR}

      echo "filter- go_advance_model existence at "`date` >> $MASTERLOG
      echo "filter- go_advance_model existence at "`date`
      ls -l go* ${CENTRALDIR}/go* >> $MASTERLOG
      ls -l go* ${CENTRALDIR}/go*

      set go_advance_exist = true
      cp -v filter_control        ${CENTRALDIR} ;# used by filter_server.csh
      cp -v go_advance_model      ${CENTRALDIR} ;# semaphore file for filter_server.csh

   endif
      
   # When the 'central' ${CENTRALDIR}/go_advance_model file disappears
   # the model has advanced to the next observation time and we need to assimilate.
   # Remove files needed for previous stage,
   # Move assim updated files for all ensemble members to local filespace.
   # Then, removing the local go_advance_model file signals filter to proceed.  
   if(! -e ${CENTRALDIR}/go_advance_model && ${go_advance_exist} == true) then

      \rm -f ${CENTRALDIR}/assim_model_state_ic*
      \rm -f               assim_model_state_ic*

      mv -v  ${CENTRALDIR}/assim_model_state_ud* .

      # signal 'filter.f90' to prep for assimilation.
      set go_advance_exist = false
      \rm -f go_advance_model
   endif
      
   # When filter creates 'go_assim_regions', copy to central directory to 
   # signal assim_regions to go.  If it already exists, the assimilation is
   # running, we just need to wait for it to finish.
   if(-e go_assim_regions && ${go_assim_exist} == false) then
      \rm -f ${CENTRALDIR}/assim_model_state_ud*
      \rm -f               assim_model_state_ud*

      mv -v  filter_assim_region__in*  ${CENTRALDIR}
      mv -v  assim_region_control      ${CENTRALDIR}
      cp -pv filter_assim_obs_seq      ${CENTRALDIR}
      cp -pv go_assim_regions          ${CENTRALDIR}
      set go_assim_exist = true
   endif
      
   # When the go_assim_regions disappears from the central directory, 
   # remove it from here to signal filter to proceed.
   # Copy updated files for all ensemble members to filesystem on this node.
   if(! -e ${CENTRALDIR}/go_assim_regions && ${go_assim_exist} == true) then
      \rm -f ${CENTRALDIR}/filter_assim_region__in*
      \rm -f               filter_assim_region__in* 

      mv -v ${CENTRALDIR}/filter_assim_region_out* .
      set go_assim_exist = false
      \rm -f go_assim_regions 
   endif

   #----------------------------------------------------------------------------
   # When filter writes out go_end_filter... 
   #----------------------------------------------------------------------------
   if(-e go_end_filter) then
      #wait for new filter_ics to appear,
      set msec = 1
      set go = no
      while ($go == no)
         echo "filter- waiting for new filter_ics to appear."
         if ( -s filter_ic_new || -s filter_ic_new.0001 ) then
            set go = yes
         else
            sleep $msec
            if ($msec < 8) @ msec = 2 * $msec
         endif
      end

      echo "filter- filter_ic_new has size" >> $MASTERLOG
      echo "filter- filter_ic_new has size"
      ls -lt filter_ic_new*                 >> $MASTERLOG
      ls -lt filter_ic_new*

      mv -v filter_ic_new*                     ${CENTRALDIR}
      mv -v assim_ic_new                       ${CENTRALDIR}

      echo "filter- copied filter_ic_new to ${CENTRALDIR}" >> $MASTERLOG
      echo "filter- copied filter_ic_new to ${CENTRALDIR}"
      ls -lt ${CENTRALDIR}/filter_ic*                      >> $MASTERLOG 
      ls -lt ${CENTRALDIR}/filter_ic*

      mv -v Prior_Diag.nc Posterior_Diag.nc                   ${CENTRALDIR}
      mv -v obs_seq.final                                     ${CENTRALDIR}

      # signal job.csh that filter is done with this obs_seq.out/day.
      cp -pv go_end_filter ${CENTRALDIR}/go_end_filter_server
      cp -pv go_end_filter ${CENTRALDIR}/go_end_filter

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

mv dart_log.out ${CENTRALDIR}/filter_log.out
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
