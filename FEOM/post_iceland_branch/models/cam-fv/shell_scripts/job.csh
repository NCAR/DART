#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$

#----------------------------------------------------------
# Script to run whole assimilation experiment.
# Submits filter.csh and filter_server.csh to compute nodes 
#    for each obs_seq.out file to be used.
# Runs interactively on master node in the central directory, where scripts reside 
#    and where script and program I/O are passed through.
#
# written 11/1/04 by Kevin Raeder
# revised 11/11/04   Kevin Raeder
#----------------------------------------------------------

# Run parameters to change

# Directory where output will be kept
set exp = EFG_new_deflts

# 'day'/obs_seq.out numbers to assimilate during this job
set obs_seq_1 = 1
set obs_seq_n = 7

# The "day" of the first obs_seq.out file of the experiment
set obs_seq_first = 1

# Subdirectory root name where each "day"'s output will be kept
# Note that you may want to change it for obs_seq_ > 9.
# Obs_diag restricts these subdirectory names to be 5 characters.
set output_dir = 01_0
set num_ens = 80
set obs_seq_1_ic  = /scratch/cluster/raeder/New_state/T21x80/03-01-01/DART
set obs_seq_1_cam = /scratch/cluster/raeder/New_state/T21x80/03-01-01/CAM/caminput_
set obs_seq_1_clm = /scratch/cluster/raeder/New_state/T21x80/03-01-01/CLM/clminput_
set input = input_

# Not currently used; gzip all previous days initial files at specified interval
set save_freq = 2
set mod_save = 1

#----------------------------------------------------------
echo $exp $num_ens $obs_seq_1 $obs_seq_n >! run_job.log
echo obs_seq_1_ic is $obs_seq_1_ic >> run_job.log

if (-d ${exp}) then
   echo ${exp} already exists >> run_job.log
else
   mkdir ${exp}
   cp namelistin ${exp}/namelistin
endif

# clean up old CAM inputs that may be laying around
if (-e caminput_1.nc) then
   rm clminput_[1-9]*.nc 
   rm caminput_[1-9]*.nc 
endif

# Have an overall outer loop over obs_seq.out files
set i = $obs_seq_1
while($i <= $obs_seq_n)
   echo ' ' >> run_job.log
   echo ' ' >> run_job.log
   echo starting iteration $i >> run_job.log

   @ j = $i - 1

#-----------------------
   # Get filter input files

   cp $input$i.nml input.nml

   # get rid of previous link
   rm obs_seq.out
   if ($i < 10) then
      ln -s /scratch/cluster/raeder/T21x80/obs_seq_jan0${i}_unformatted.out obs_seq.out
   else
      ln -s /scratch/cluster/raeder/T21x80/obs_seq_jan${i}_unformatted.out obs_seq.out
   endif
   echo job- obs_seq used is >> run_job.log
   ls -lt obs_seq.out >> run_job.log

   if ($i == $obs_seq_first) then
      # DART initial condition(s)
      set from_root = ${obs_seq_1_ic}

#     These are where CAM gets it's initial files, in advance_model.csh
      set cam_init = $obs_seq_1_cam
      set clm_init = $obs_seq_1_clm
   else
      set from_root = `pwd`/$exp/${output_dir}${j}/DART
      set cam_init =  `pwd`/$exp/${output_dir}${j}/CAM/caminput_
      set clm_init =  `pwd`/$exp/${output_dir}${j}/CLM/clminput_
   endif

   # transmit info to advance_model.csh, run by filter_server.csh
   # The second item echoed must be the subdirectory where CAM is kept, in the Central directory
   echo $exp cam3.0.7 $cam_init $clm_init >! casemodel

   rm assim_ic_old
   ln -s $from_root/assim_tools_ics assim_ic_old
   # link to filter_ic file(s), so that filter.csh can copy them to a compute node
   if (-e ${from_root}/filter_ic.0001) then
      set n = 1
      while($n <= ${num_ens})
           set from = ${from_root}/filter_ic*[.0]$n
           rm filter_ic_old.$from:e
           ln -s $from filter_ic_old.$from:e
           @ n++ 
      end
   else
      rm filter_ic_old
      ln -s $from_root/filter_ic filter_ic_old
   endif
   echo ' ' >> run_job.log
   echo job- filter_ic_old is/are >> run_job.log
   ls -lt filter_ic_old* >> run_job.log

#-----------------------
   # make subdirectories to store this obs_seq_s output
   mkdir ${exp}/${output_dir}$i
   mkdir ${exp}/${output_dir}$i/{CLM,CAM,DART}

#-----------------
   # Run the filter
   # This is the central directory for whole filter job
   #    qsub jobs submitted from there
   #    semaphor files must (dis)appear there, 
   #    I/O between filter and advance_model and assim_region goes through there.
   #    Final output is put there
   # It is PBS_O_WORKDIR  in filter.csh and filter_server.csh

   # advances model and assims regions
   qsub filter_server.csh 

   # runs filter, which integrates the results of model advances and region assims
   qsub filter.csh  
#-----------------

   # Hang around forever for now and wait for filter.csh to finish
   # If go_end_filter exists then stop this process

   set again = true
   set exist_end_filter = false
   set nsec = 1
   
   while($again == true)
      if(-e go_end_filter && ${exist_end_filter} == false) then
         # do this when go_end_filter first appears
         set exist_end_filter = true 
         # go_end_filter also signals filter_server to finish, but give it a second.
         sleep 1
         # Now signal this obs_seq.out loop iteration to finish
         rm go_end_filter
      else if (! -e go_end_filter && ${exist_end_filter} == true ) then
         # do this when go_end_filter first disappears
         set again = false
         echo "job.csh finishing day $i normally at " `date` >> run_job.log
         echo does PBS_O_WORKDIR contain filter_ic_new >> run_job.log
         ls -lt filter_ic_new* >> run_job.log
      else
         # do this while waiting for go_end_filter to first appear
         sleep $nsec
         if ($nsec < 8) @ nsec = 2 * $nsec
      endif
   end

#-----------------------------------------------
   # Move the output to storage after filter.csh signals 
   # that it's done moving this stuff here

   mv clminput_[1-9]*.nc                 ${exp}/${output_dir}$i/CLM
   mv caminput_[1-9]*.nc                 ${exp}/${output_dir}$i/CAM
   mv Prior_Diag.nc Posterior_Diag.nc    ${exp}/${output_dir}$i
   mv obs_seq.final                      ${exp}/${output_dir}$i

   # Move the filter restart file(s) to the storage subdirectory
   if (-e filter_ic_new) then
      mv filter_ic_new ${exp}/${output_dir}$i/DART/filter_ic
      echo moving filter_ic_new to ${exp}/${output_dir}$i/DART/filter_ic >> run_job.log
   else if (-e filter_ic_new.0001) then
      set n = 1
      while($n <= ${num_ens})
           set from = filter_ic_new*[.0]$n
           mv $from $exp/${output_dir}$i/DART/filter_ic.$from:e
           @ n++
      end
      echo moving filter_ic_newS to ${exp}/${output_dir}$i/DART/filter_icS >> run_job.log
   else
      echo NO filter_ic_new FOUND >> run_job.log
   endif
   if (-e assim_ic_new) then
      mv assim_ic_new ${exp}/${output_dir}$i/DART/assim_tools_ics
   endif

   # test whether it's safe to end this obs_seq_ by signalling filter.csh to end
   # kluge ('00num_ens') to prevent losing restart files only good for >9 ens members
   # -e doesn't care if it has size 0, so to be safer test for size 0
   # -z means if it exists AND has size 0 
   # ! -z means that it exists with size 0, OR it doesn't exist
   # We want only the first possibility, hence the additional -e test
   if (  -e $exp/${output_dir}$i/DART/filter_ic || \
       (  -e $exp/${output_dir}$i/DART/filter_ic.00${num_ens} && \
        ! -z $exp/${output_dir}$i/DART/filter_ic.00${num_ens}     )) then
       echo okay > rm_filter_temp
   else
#      This causes loop/job to exit gracefully
       @ i = $obs_seq_n + 1
       echo RETRIEVE filter_ic files from filter.csh temp directory   >> run_job.log
       echo Then remove temp and cam advance temps >> run_job.log
   endif

   # Compress older output which we won't need immediately
   #
   # if ($i % $save_freq != $mod_save) then
   # if ($i > 1) then
   #    @ j = $i - 1
   #    cd ${exp}/${output_dir}$j
   #    if (-e ${exp}/${output_dir}$i/CLM/clminput_${num_ens}.nc) then
   #       gzip -r CLM
   # #      tar cf clminput.gz.tar CLM
   # #      rm -rf CLM
   #    else
   #       echo 'NO clminput.20;  ABORTING compression' >> run_job.log
   #    endif       
   #    if (-e ${exp}/${output_dir}$i/CAM/caminput_${num_ens}.nc) then
   #       gzip -r CAM
   # #      tar cf caminput.gz.tar CAM
   # #      rm -rf CAM
   #    else
   #       echo 'NO caminput.20;  ABORTING compression' >> run_job.log
   #    endif
   #    cd ../..
   # endif


   mv cam_out_temp1       ${exp}/${output_dir}$i
   mv cam_reg_temp1       ${exp}/${output_dir}$i
   mv input.nml           ${exp}/${output_dir}$i
   mv casemodel           ${exp}/${output_dir}$i
   mv run_filter.*        ${exp}/${output_dir}$i
   mv filter_serv*.[el]*  ${exp}/${output_dir}$i
   mv filter.out          ${exp}/${output_dir}$i

# end i loop
   @ i++
end

mv run_job.* $exp
mv namelist $exp
mv dart_out.log $exp


