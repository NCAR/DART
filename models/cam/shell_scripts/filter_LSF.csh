#!/bin/csh -f
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$
#
#-------------------------
# Script to run the filter executable and signals filter_server to advance the model
#    or assimilate the regions
# Submitted to batch queue by job.csh (the main control script).
# Executes filter on a compute node, distinct from those used by filter_server.csh
# 
# written 11/1/04 Kevin Raeder
# revised 11/11/04 Kevin Raeder
# modified by committee 3/31/05 Kevin Raeder, Alain Caya, typing by Tim Hoar
#-------------------------
#
#### LSF options for BSUB
### -J      job name
### -o      output listing filename 
### -q      queue
### Queue charging    cheapest ..... most expensive   (at least on lightning)
### Queue name (standby, economy, [regular,debug], premium)
#
#BSUB -J run_filter
#BSUB -o run_filter.%J.log
#BSUB -q debug
#BSUB -n 1

# Central directory for whole filter job (DARTWORKDIR)
#    batch jobs submitted from there
#    semaphor files must (dis)appear there, 
#    I/O between filter and advance_model and assim_region goes through there.
#    Final output is put there
if ( $?LS_SUBCWD ) then           ;# batch job
   # the central experiment directory is automatically defined
   set DARTWORKDIR = $LS_SUBCWD
   echo "filter- running on $LSB_HOSTS" >> $DARTWORKDIR/run_job.log
else
   echo "filter- running locally in "   >> $DARTWORKDIR/run_job.log
   setenv DARTWORKDIR `pwd`
   echo "filter- $DARTWORKDIR"          >> $DARTWORKDIR/run_job.log
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

# This job's working directory. 
# We are ensuring it is empty when we start, and (later) GONE when we finish.
set tempdir = /ptmp/$user/tmp$$
\rm -rf $tempdir
mkdir $tempdir
echo "filter- cd to $tempdir" >> $DARTWORKDIR/run_job.log
echo "cd to $tempdir"          > $DARTWORKDIR/run_filter.stout
cd $tempdir

#-----------------------
# Get filter input files
cp $DARTWORKDIR/input.nml      .
cp $DARTWORKDIR/caminput.nc    .
cp $DARTWORKDIR/obs_seq.out    .
cp $DARTWORKDIR/filter         .
cp $DARTWORKDIR/filter_ic_old* .
cp $DARTWORKDIR/assim_ic_old   .

#-----------------------
# Run the filter
./filter >> $DARTWORKDIR/run_filter.stout &

#-----------------
# Hang around forever for now and wait for go_xxx to appear here,
# or disappear from $DARTWORKDIR

set again = true
# flags to keep track of when semaphor files disappear.
set go_advance_exist = false
set go_assim_exist = false
set nsec = 1

# debug
# set first_time = true

while($again == true)

   # When go_advance_model appears, copy it to the common directory to signal
   # advance_model to go.  If it already exists, do nothing.

   if(-e go_advance_model && ${go_advance_exist} == false) then
      # remove files needed for previous stage, to conserve disc space
      if (-e filter_assim_region_out1) then
         rm $DARTWORKDIR/filter_assim_region_out* filter_assim_region_out*
      endif

      # write size of assim_model_state_ic to Central directory for filter_server.csh to use.
      set list = `ls -lt assim_model_state_ic1`
      echo $list[5] >! $DARTWORKDIR/assim_size

      cp filter_control   $DARTWORKDIR
      cp go_advance_model $DARTWORKDIR

      echo "go_advance_model existence" >> $DARTWORKDIR/run_job.log
      ls -l go* $DARTWORKDIR/go*        >> $DARTWORKDIR/run_job.log
      # copy in numerical order, not alphabetic
      set n = 1
      while (-e assim_model_state_ic$n)
         cp assim_model_state_ic$n $DARTWORKDIR
         @ n++
      end
      set go_advance_exist = true

   endif
      
  # When the 'central' $DARTWORKDIR/go_advance_model file disappears
  # Copy assim updated files for all ensemble members to filesystem on this node.
  # Then, removing the local go_advance_model file signals filter to proceed.  
  if(! -e $DARTWORKDIR/go_advance_model && ${go_advance_exist} == true) then
      # remove files needed for previous stage, to conserve disc space
      rm assim_model_state_ic* $DARTWORKDIR/assim_model_state_ic*
      cp $DARTWORKDIR/assim_model_state_ud* .
      set go_advance_exist = false
      rm go_advance_model 
   endif
      
   # When filter writes out go_assim_regions, copy it to the central directory to 
   # signal assim_regions to go.  If it already exists, do nothing.
   if(-e go_assim_regions && ${go_assim_exist} == false) then
      rm assim_model_state_ud*    $DARTWORKDIR/assim_model_state_ud*
      cp filter_assim_region__in* $DARTWORKDIR
      cp assim_region_control     $DARTWORKDIR
      cp filter_assim_obs_seq     $DARTWORKDIR
      cp go_assim_regions         $DARTWORKDIR
      set go_assim_exist = true
   endif
      
   # When the go_assim_regions disappears from the central directory, remove it
   # from here to signal filter to  proceed.
   # Copy updated files for all ensemble members to filesystem on this node.
   if(! -e $DARTWORKDIR/go_assim_regions && ${go_assim_exist} == true) then
      # remove files needed for previous stage, to conserve disc space
      rm filter_assim_region__in* $DARTWORKDIR/filter_assim_region__in*
      cp $DARTWORKDIR/filter_assim_region_out* .
      set go_assim_exist = false
      rm go_assim_regions 
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
      echo "In filter.csh filter_ic_new has size"         >> $DARTWORKDIR/run_job.log
      ls -lt filter_ic_new                                >> $DARTWORKDIR/run_job.log
      cp filter_ic_new*                     $DARTWORKDIR
      cp assim_ic_new                       $DARTWORKDIR
      echo "filter- copied filter_ic_new to $DARTWORKDIR" >> $DARTWORKDIR/run_job.log
      ls -lt $DARTWORKDIR/filter_ic*                      >> $DARTWORKDIR/run_job.log
      cp Prior_Diag.nc Posterior_Diag.nc    $DARTWORKDIR
      cp obs_seq.final                      $DARTWORKDIR
      cp filter.out                         $DARTWORKDIR

      # signal job.csh that filter is done with this obs_seq.out/day.
      mv go_end_filter $DARTWORKDIR/go_end_filter_server
      echo "filter- finished assimilating"               >> $DARTWORKDIR/run_job.log
      echo "filter- moved go_end_filter to $DARTWORKDIR" >> $DARTWORKDIR/run_job.log
      ls -lt $DARTWORKDIR/go*                            >> $DARTWORKDIR/run_job.log

      # setting 'again' to false exits the loop, signals we are DONE.
      set again = false
      echo "filter- filter.csh terminating normally at "`date`  >> $DARTWORKDIR/run_job.log
   else
      sleep $nsec
      if ($nsec < 8) @ nsec = 2 * $nsec
      echo "filter-  "`date`
   endif
end

#-----------------------------------------------
mv dart_out.log $DARTWORKDIR
cd $DARTWORKDIR

# To prevent losing restart files look for signal from job.csh that it's safe
# job.csh copies the filter_ic files to 'permanent' storage. 
# Until that is done, we don't want to remove them here. 

#set again = true
#while($again == true)
#   if (-e rm_filter_temp) then
#      rm -rf $tempdir rm_filter_temp
#      set again = false
#   else
#      echo "filter- rm_filter_temp not found yet"  >> $DARTWORKDIR/run_job.log
#      sleep 10
#   endif
#end
