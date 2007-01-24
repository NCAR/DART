#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$

#-------------------------
# Script to run the filter executable and signals filter_server to advance the model
#    or assimilate the regions
# Submitted to batch queue by job.csh (the main control script).
# Executes filter on a compute node, distinct from those used by filter_server.csh
# 
# written 11/1/04 Kevin Raeder
# last revised 11/11/04 Kevin Raeder
# last revised  9/28/05 Tim Hoar ... verbose flags
#-------------------------

### Job name
#PBS -N run_filter
### Declare job non-rerunable
#PBS -r n
### Output files
#PBS -e run_filter.err
#PBS -o run_filter.out
### Queue name (small, medium, long, verylong)
#PBS -q long
#PBS -l nodes=1:ppn=1

### This job's working directory; must cd to it, or it will run in /home...
set tempdir = /scratch/local/tmp$$$user
\rm -f $tempdir
mkdir -p $tempdir
echo "filter- cd to $tempdir" >> $PBS_O_WORKDIR/run_job.log
echo "cd to $tempdir"          > $PBS_O_WORKDIR/run_filter.stout
cd $tempdir

if ($?PBS_O_WORKDIR) then
   # it's a batch job and the central experiment directory is automatically defined
   echo "filter- running on "`cat $PBS_NODEFILE` >> $PBS_O_WORKDIR/run_job.log
else
   echo "filter- setenv PBS_O_WORKDIR"           >> $PBS_O_WORKDIR/run_job.log
   setenv PBS_O_WORKDIR `pwd`
endif

#-----------------------
# Get filter input files
cp -v $PBS_O_WORKDIR/input.nml .
cp -v $PBS_O_WORKDIR/caminput.nc .
cp -v $PBS_O_WORKDIR/obs_seq.out .
cp -v $PBS_O_WORKDIR/filter .
cp -v $PBS_O_WORKDIR/filter_ic_old* .
cp -v $PBS_O_WORKDIR/assim_ic_old .

#-----------------------
# Run the filter
./filter >> $PBS_O_WORKDIR/run_filter.stout &

#-----------------
# Hang around forever for now and wait for go_xxx to appear here (tempdir),
# or disappear from $PBS_O_WORKDIR

set again = true
# flags to keep track of when semaphor files disappear.
set go_advance_exist = false
set go_assim_exist = false
set nsec = 1

# debug
#set first_time = true

while($again == true)

   # When filter writes out go_advance_model, copy it to the central directory 
   # to signal advance_model to go.  If it already exists, do nothing.
   if(-e go_advance_model && ${go_advance_exist} == false) then
      # remove files needed for previous stage, to conserve disc space
      if (-e filter_assim_region_out1) then
         \rm -f $PBS_O_WORKDIR/filter_assim_region_out* filter_assim_region_out*
      endif
      cp -v go_advance_model $PBS_O_WORKDIR
      cp -v filter_control $PBS_O_WORKDIR
      # copy in numerical order, not alphabetic; this permits filter_server to send
      # the first batch of model advances as soon as the necessary initial conditions
      # are there
      set n = 1
      while (-e assim_model_state_ic$n)
         cp -v assim_model_state_ic$n $PBS_O_WORKDIR
         @ n++
      end
      set go_advance_exist = true
      ## debug; do only first forecast and assim cycle. 
      #         Exit gracefully to preserve files.
      #    if (${first_time} == true) then
      #       set first_time = false
      ## end debug
      ## debug
      #    else
      #       echo done > go_end_filter
      #       echo done > $PBS_O_WORKDIR/go_end_filter
      #    endif
      ## end debug
   endif
      
  # When the go_advance_model disappears from the central directory, remove it
  # from here to signal filter to  proceed.  
  if(! -e $PBS_O_WORKDIR/go_advance_model && ${go_advance_exist} == true) then
      # remove files needed for previous stage, to conserve disc space
      \rm -f assim_model_state_ic* $PBS_O_WORKDIR/assim_model_state_ic*
      cp -v $PBS_O_WORKDIR/assim_model_state_ud* .
      set go_advance_exist = false
      \rm -f go_advance_model 
   endif
      
   # When filter writes out go_assim_regions, copy it to the central directory to 
   # signal assim_regions to go.  If it already exists, do nothing.
   if(-e go_assim_regions && ${go_assim_exist} == false) then
      # remove files needed for previous stage, to conserve disc space
      \rm -f assim_model_state_ud* $PBS_O_WORKDIR/assim_model_state_ud*
      cp -v filter_assim_region__in* $PBS_O_WORKDIR
      cp -v assim_region_control $PBS_O_WORKDIR
      cp -v filter_assim_obs_seq $PBS_O_WORKDIR
      set go_assim_exist = true
      cp -v go_assim_regions $PBS_O_WORKDIR
   endif
      
   # When the go_assim_regions disappears from the central directory, remove it
   # from here to signal filter to  proceed.
   if(! -e $PBS_O_WORKDIR/go_assim_regions && ${go_assim_exist} == true) then
      # remove files needed for previous stage, to conserve disc space
      \rm -f filter_assim_region__in* $PBS_O_WORKDIR/filter_assim_region__in*
      cp -v $PBS_O_WORKDIR/filter_assim_region_out* .
      set go_assim_exist = false
      \rm -f go_assim_regions 
   endif

   # When filter writes out go_end_filter... 
   if(-e go_end_filter) then
      #wait for new filter_ics to appear,
      set msec = 1
      set go = no
      while ($go == no)
         ls filter_ic_new* >! .garb_ic
         if ($status == 0) then
            set go = yes
         else
            sleep $msec
            if ($msec < 8) @ msec = 2 * $msec
         endif
      end

      # move output to central directory for analysis and storage
      cp -v filter_ic_new*                          $PBS_O_WORKDIR
      cp -v assim_ic_new                            $PBS_O_WORKDIR
      echo "filter- copied filter_ic_new to PBS" >> $PBS_O_WORKDIR/run_job.log
      ls -lt $PBS_O_WORKDIR/filter_ic*           >> $PBS_O_WORKDIR/run_job.log
      cp -v Prior_Diag.nc Posterior_Diag.nc         $PBS_O_WORKDIR
      cp -v obs_seq.final                           $PBS_O_WORKDIR
      cp -v filter.out                              $PBS_O_WORKDIR

      # signal job.csh that filter is done with this obs_seq.out/day.
      mv -v go_end_filter $PBS_O_WORKDIR
      echo "filter- moved go_end_filter to PBS_O_WORKDIR" >> $PBS_O_WORKDIR/run_job.log
      ls -lt $PBS_O_WORKDIR/go*                           >> $PBS_O_WORKDIR/run_job.log

      set again = false
      echo "filter- filter.csh terminating normally at "`date`  >> $PBS_O_WORKDIR/run_job.log
   else
      sleep $nsec
      if ($nsec < 8) @ nsec = 2 * $nsec
      date
   endif
end

#-----------------------------------------------
mv -v dart_out.log $PBS_O_WORKDIR

# to prevent losing restart files look for signal from job.csh that it's safe
cd $PBS_O_WORKDIR
   
set again = true
while($again == true)
   if (-e rm_filter_temp) then
      \rm -rf $tempdir rm_filter_temp
      set again = false
   else
      echo "filter- rm_filter_temp not found yet"  >> $PBS_O_WORKDIR/run_job.log
      sleep 10
   endif
end

