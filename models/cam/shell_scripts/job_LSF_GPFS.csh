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

#-------------------------------------------------------------------------------
# Script to run whole assimilation experiment.
# Submits batch jobs 'filter.csh' and 'filter_server.csh'
#    for each obs_seq.out file to be processed.
# Runs interactively on master node in the central directory, where scripts reside 
#    and where script and program I/O are passed through.
#-------------------------------------------------------------------------------
#
#### LSF options for BSUB
### -J      job name    (master script job.csh presumes filter.xxxx.log)
### -o      output listing filename 
### -P      account number
### -q      queue
### -n      number of tasks (processors)
### -x      exclusive use of node
### -R "span[ptile=(num procs you want on each node)]"
#
#BSUB -J DARTCAM
#BSUB -o DARTCAM.%J.log
#BSUB -q standby
#BSUB -n 1

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

if ( $?LS_SUBCWD ) then
   set CENTRALDIR = $LS_SUBCWD
else
   setenv CENTRALDIR `pwd`
endif

# Set Variable for a 'master' logfile 
set MASTERLOG = ${CENTRALDIR}/run_job.log

echo "job- running    on ${PROCNAMES}"     > $MASTERLOG
echo "job- CENTRALDIR is ${CENTRALDIR}"   >> $MASTERLOG
echo "job- running    on ${PROCNAMES}"
echo "job- CENTRALDIR is ${CENTRALDIR}"

# Run parameters to change

# Directory where output will be kept (relative to '.')
set exp = Experiment1

# 'day'/obs_seq.out numbers to assimilate during this job
# First and last obs seq files. 
set obs_seq_1 = 1
set obs_seq_n = 1

# number of obs_seqs / day
set obs_seq_freq = 2

# The "day" of the first obs_seq.out file of the experiment
set obs_seq_first = 1

# Subdirectory root name where each "day"'s output will be kept
# Obs_diag restricts these subdirectory names to be 5 characters.
# output_root = xxx    xxx signifies the month.
set output_root = 01_

set num_ens = 80

# Initial conditions for DART-CAM for the first obs_seq.out file
set obs_seq_1_ic  = /ncar/dart/CAM_init/T21x80/03-01-01/DART_lunes
set obs_seq_1_cam = /ncar/dart/CAM_init/T21x80/03-01-01/CAM/caminput_
set obs_seq_1_clm = /ncar/dart/CAM_init/T21x80/03-01-01/CLM/clminput_

# Get the CAM executable from a sub-directory of the whole CAM package to be used
set CAMsrc = /home/lightning/raeder/Cam3/cam3.1/models/atm/cam/bld/T21


# Each obs_seq file gets its own input.nml ... the following variable
# helps set this up. 


set input = input_


# CHANGE choice of whether forecasts on caminput_##.nc files should be replaced with
#        analyses from filter_ic_new.#### at the end of each obs_seq.  See 'stuff' below

# CHANGE names of obs_seq files below, if necessary


# Not currently used; all previous days initial files are gzipped
# Now don't even gzip those; keep available for other restarts
set save_freq = 2
set mod_save = 1

#----------------------------------------------------------
echo "exp num_ens obs_seq_1 obs_seq_n obs_seq_first"
echo "$exp $num_ens $obs_seq_1 $obs_seq_n $obs_seq_first"
echo "obs_seq_1_ic is $obs_seq_1_ic"

echo "exp num_ens obs_seq_1 obs_seq_n obs_seq_first"       >> $MASTERLOG
echo "$exp $num_ens $obs_seq_1 $obs_seq_n $obs_seq_first"  >> $MASTERLOG
echo "obs_seq_1_ic is $obs_seq_1_ic"                       >> $MASTERLOG

if (-d ${exp}) then
   echo "${exp} already exists" >> $MASTERLOG
   echo "${exp} already exists"
else
   echo "Making run-time directory $exp ..." >> $MASTERLOG
   echo "Making run-time directory $exp ..."
   mkdir -p ${exp}
   if (! -e namelistin ) then
      echo "ERROR ... need a namelistin file." >> $MASTERLOG
      echo "ERROR ... need a namelistin file."
      exit 99
   endif
   cp namelistin ${exp}/namelistin
endif

# clean up old CAM inputs that may be laying around
if (-e caminput_1.nc) then
   rm clminput_[1-9]*.nc 
   rm caminput_[1-9]*.nc 
endif

# Have an overall outer loop over obs_seq.out files
set i = $obs_seq_1
while($i <= $obs_seq_n) ;# start i loop
   echo ' '
   echo ' ' >> $MASTERLOG
   echo "starting observation sequence file $i at "`date`
   echo "starting observation sequence file $i at "`date` >> $MASTERLOG


   @ j = $i - 1
   set out_prev = ${output_root}
   if ($j < 10) set out_prev = ${output_root}0
   set out_prev = ${out_prev}$j

   set output_dir = ${output_root}
   if ($i < 10) set output_dir = ${output_root}0
   set output_dir = ${output_dir}$i

   #-----------------------
   # Get filter input files
   # The first one is different than all the rest ...

#   if ($i == 1) then
   if ($i == $obs_seq_first) then
      cp ${input}${obs_seq_first}.nml input.nml
   else
      if (-e ${input}n.nml) then
         cp ${input}n.nml input.nml
      else if (-e ${input}${i}.nml) then
         cp ${input}${i}.nml input.nml
      else
         echo "input_next is MISSING" >> $MASTERLOG
         echo "input_next is MISSING"
         exit
      endif
   endif

   # get rid of previous link
   rm obs_seq.out
   
   #----------------------
   # Get obs_seq file for this assimilation
   # set day and hour for this obs_seq
   # Normally, when starting with 12 hourly from beginning;  
   @ day = ($i + 1) / $obs_seq_freq
   @ hour = ($i % $obs_seq_freq) * 24 / $obs_seq_freq 
   if ($hour == 0) set hour = 24
   echo "obs_seq, day, hour = $i $day $hour " >> $MASTERLOG

   if ($i == 0) then
      set OBS_SEQ = /scratch/cluster/raeder/GWD_T42/obs_seq_1-1-03.out
   else if ($i < 19) then
      set OBS_SEQ = /ncar/dart/Obs_sets/Allx12/obs_seq2003010${day}${hour}
   else
      set OBS_SEQ = /ncar/dart/Obs_sets/Allx12/obs_seq200301${day}${hour}

## obs_seq_freq = 24      set OBS_SEQ = /ncar/dart/Obs_sets/All/obs_seq200301${i}

   endif

   if ( -s $OBS_SEQ ) then    ;# -s is true if xxxx has non-zero size
      ln -s $OBS_SEQ  obs_seq.out
   else
      echo "ERROR - no obs_seq.out for $i - looking for:" 
      echo $OBS_SEQ 
      exit 123
   endif

   echo "job- obs_seq $i used is $OBS_SEQ" >> $MASTERLOG
   echo "job- obs_seq $i used is $OBS_SEQ"
   ls -lt obs_seq.out

   if ($i == $obs_seq_first) then 
      # get initial conditions for DART, model from 'permanent' storage
      # from_root defines location of DART initial condition(s)
      # cam_init, clm_init are where CAM gets it's initial files (in advance_model.csh)
      set from_root = ${obs_seq_1_ic} 
      set  cam_init = $obs_seq_1_cam
      set  clm_init = $obs_seq_1_clm
   else
      # get initial files from result of previous experiment.
      set from_root = `pwd`/$exp/${out_prev}/DART
      set cam_init =  `pwd`/$exp/${out_prev}/CAM/caminput_
      set clm_init =  `pwd`/$exp/${out_prev}/CLM/clminput_
   endif

   # transmit info to advance_model.csh, run by filter_server.csh
   # The second item echoed must be the subdirectory where CAM is kept, in the Central directory
   echo "$exp $CAMsrc $cam_init $clm_init" >! casemodel

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
   else if (-e ${from_root}/filter_ic) then
      rm filter_ic_old
      ln -s ${from_root}/filter_ic filter_ic_old
   endif
   echo ' '
   echo "job- filter_ic_old is/are" >> $MASTERLOG
   echo "job- filter_ic_old is/are"
   ls -lt filter_ic_old*            >> $MASTERLOG
   ls -lt filter_ic_old*

   #-----------------------
   # make subdirectories to store this obs_seq_s output
   mkdir -p ${exp}/${output_dir}
   mkdir -p ${exp}/${output_dir}/{CLM,CAM,DART}

   #-----------------
   # Run the filter
   # This is the central directory for whole filter job
   #    qsub jobs submitted from there
   #    semaphor files must (dis)appear there, 
   #    I/O between filter and advance_model and assim_region goes through there.
   #    Final output is put there
   # It's CENTRALDIR  in filter.csh and filter_server.csh

   # advances model and assims regions (do this first to grab whole nodes)
   # We need to capture the batch job number to kill later if need be.
   bsub < filter_server.csh > batchsubmit$$
   set STRING = "1,$ s#<##g"
   sed -e "$STRING" batchsubmit$$ > bill$$
   set STRING = "1,$ s#>##g"
   sed -e "$STRING" bill$$ > batchsubmit$$
   set STRING = `cat batchsubmit$$`
   set FILTERSERVERBATCHID = $STRING[2]
   rm -f batchsubmit$$ bill$$

   # runs filter, which integrates the results of model advances and region assims
   # This only uses 1 processor, so it could fit on a node with someone else.
   bsub < filter.csh  > batchsubmit$$
   set STRING = "1,$ s#<##g"
   sed -e "$STRING" batchsubmit$$ > bill$$
   set STRING = "1,$ s#>##g"
   sed -e "$STRING" bill$$ > batchsubmit$$
   set STRING = `cat batchsubmit$$`
   set FILTERBATCHID = $STRING[2]
   rm -f batchsubmit$$ bill$$

   echo "filter_server spawned job $FILTERSERVERBATCHID at "`date`
   echo "filter        spawned job $FILTERBATCHID at "`date`
   echo "filter_server spawned job $FILTERSERVERBATCHID at "`date` >> $MASTERLOG
   echo "filter        spawned job $FILTERBATCHID at "`date`       >> $MASTERLOG

   set KILLCOMMAND = "bkill $FILTERSERVERBATCHID $FILTERBATCHID; touch BOMBED; exit"

   #-----------------
   # When filter.f90 finished, a file called 'go_end_filter' is created in the
   # filter.csh run-time directory (i.e. not CENTRALDIR).
   # Filter.csh copies that local semaphore file to CENTRALDIR to signal that 
   # it is done with this observation sequence file. 

   set again = true
   set nsec = 1
   
   while($again == true)
      if( -e go_end_filter ) then
         echo "job.csh finishing $OBS_SEQ at "`date` >> $MASTERLOG
         echo "job.csh finishing $OBS_SEQ at "`date`
         echo "does CENTRALDIR contain filter_ic_new"
         ls -lt filter_ic_new*
         echo "does CENTRALDIR contain filter_ic_new"        >> $MASTERLOG
         ls -lt filter_ic_new*                               >> $MASTERLOG

         rm -v go_end_filter   ; # do this when go_end_filter first appears
         set again = false     ; # gets us out of this perpetual loop
      else
         # do this while waiting for go_end_filter to first appear
         sleep $nsec
         if ($nsec < 8) then
            @ nsec = 2 * $nsec
            echo "job.csh waiting for go_end_filter to appear "`date` >>$MASTERLOG
            echo "job.csh waiting for go_end_filter to appear "`date`
         endif
      endif
   end

   #-----------------------------------------------
   # Move the output to storage after filter.csh signals completion.
   # At this point, all the restart,diagnostic files are in the CENTRALDIR
   # and need to be moved to the 'experiment permanent' directory.
   # We have had problems with some, but not all, files being moved
   # correctly, so we are adding bulletproofing to check to ensure the filesystem
   # has completed writing the files, etc. Sometimes we get here before
   # all the files have finished being written.
   #-----------------------------------------------
   # This was clean, but did not always work ... sigh ...
   # mv clminput_[1-9]*.nc    ${exp}/${output_dir}/CLM
   # mv caminput_[1-9]*.nc    ${exp}/${output_dir}/CAM
   #-----------------------------------------------

   echo "Listing contents of CENTRALDIR before archiving at "`date`
   ls -l
   echo "Listing contents of CENTRALDIR before archiving at "`date` >> $MASTERLOG
   ls -l >> $MASTERLOG

   foreach FILE ( Prior_Diag.nc Posterior_Diag.nc obs_seq.final )
      if ( -s $FILE ) then
         mv -v $FILE ${exp}/${output_dir} 
         if ( ! $status == 0 ) then
            echo "failed moving ${CENTRALDIR}/$FILE" >>$MASTERLOG
            echo "failed moving ${CENTRALDIR}/$FILE"
            $KILLCOMMAND
         endif
      else
         echo "... THUMP ... ${CENTRALDIR}/$FILE does not exist and should."
         echo "directory contents follows"
         ls -l
         $KILLCOMMAND
      endif
   end

   # Move the filter restart file(s) to the storage subdirectory
   echo "moving filter_ic_newS to ${exp}/${output_dir}/DART/filter_icS"

   if (-e filter_ic_new) then
      mv -v filter_ic_new ${exp}/${output_dir}/DART/filter_ic
      if (! $status == 0 ) then
         echo "failed moving filter_ic_new to ${exp}/${output_dir}/DART/filter_ic"
         $KILLCOMMAND
      endif
   else if (-e filter_ic_new.0001) then
      set n = 1
      while($n <= ${num_ens})
           set from = filter_ic_new*[.0]$n
           set dest = $exp/${output_dir}/DART/filter_ic.$from:e
           # # stuff analyses into CAM initial files using
           # echo $from >! member
           # echo caminput_${n}.nc >> member
           # ./trans_sv_pv
           # echo "stuffing analyses from $from into caminput_${n}.nc" >> $MASTERLOG
           # ls -l caminput_${n}.nc >> $MASTERLOG
           # # end stuffing
           mv -v $from $dest
           if (! $status == 0 ) then
              echo "failed moving $from to ${dest}"
              $KILLCOMMAND
           endif
           @ n++
      end
   else
      echo "NO filter_ic_new FOUND"
      echo "NO filter_ic_new FOUND"
      echo "NO filter_ic_new FOUND"
      echo "NO filter_ic_new FOUND"
      $KILLCOMMAND
   endif

   if (-e assim_ic_new) then
      mv -v assim_ic_new ${exp}/${output_dir}/DART/assim_tools_ics
      if (! $status == 0 ) then
         echo "failed moving assim_ic_new to ${exp}/${output_dir}/DART/assim_tools_ics"
         $KILLCOMMAND
      endif
   endif

   set n = 1
   while($n <= ${num_ens})    ;# loop over all ensemble members ... one-at-a-time
      set CAMINPUT = caminput_${n}.nc 
      set CLMINPUT = clminput_${n}.nc 

      if ( -s $CAMINPUT ) then
         mv -v $CAMINPUT ${exp}/${output_dir}/CAM 
         if (! $status == 0 ) then
            echo "failed moving ${CENTRALDIR}/$CAMINPUT"
            $KILLCOMMAND
         endif
      else
         echo "... aieeee ... ${CENTRALDIR}/$CAMINPUT does not exist and should."
         $KILLCOMMAND
      endif

      if ( -s $CLMINPUT ) then
         mv -v $CLMINPUT ${exp}/${output_dir}/CLM 
         if (! $status == 0 ) then
            echo "failed moving ${CENTRALDIR}/$CLMINPUT"
            $KILLCOMMAND
         endif
      else
         echo "... aieeee ... ${CENTRALDIR}/$CLMINPUT does not exist and should."
         $KILLCOMMAND
      endif

     @ n++
   end


   # test whether it's safe to end this obs_seq_ by signalling filter.csh to end
   # kluge ('00num_ens') to prevent losing restart files only good for >9 ens members
   # -e doesn't care if it has size 0, so to be safer test for size 0
   # -z means if it exists AND has size 0 
   # ! -z means that it exists with size 0, OR it doesn't exist
   # We want only the first possibility, hence the additional -e test
   if (   -s $exp/${output_dir}/DART/filter_ic || \
          -s $exp/${output_dir}/DART/filter_ic.00${num_ens} ) then
       echo "Signalling filter.csh that it is OK to exit at "`date` >> $MASTERLOG
       echo "Signalling filter.csh that it is OK to exit at "`date`
       echo okay > rm_filter_temp
   else
#      This causes loop/job to exit gracefully
       @ i = $obs_seq_n + 1
       echo "RETRIEVE filter_ic files from filter.csh temp directory"
       echo "Then remove temp and cam advance temps"
       echo "RETRIEVE filter_ic files from filter.csh temp directory" >> $MASTERLOG
       echo "Then remove temp and cam advance temps" >> $MASTERLOG
   endif

   # Compress and archive older output which we won't need immediately
   #
   # if ($i % $save_freq != $mod_save) then
   if ($i > $obs_seq_first ) then
       cd ${exp}/${out_prev}
       bsub < ~/auto_re2ms_LSF.csh                                         >>& $MASTERLOG
       cd ../..
       echo "Backing up restart $j to mass store;  in separate batch job"  >> $MASTERLOG
    endif


   mv -v cam_out_temp1       ${exp}/${output_dir}   ;# save a representative model advance
   mv -v cam_reg_temp1       ${exp}/${output_dir}   ;# ditto for assimilation region
   mv -v input.nml           ${exp}/${output_dir}
   mv -v casemodel           ${exp}/${output_dir}
   mv -v run_filter.stout    ${exp}/${output_dir}   ;# stdout   from filter.f90
   mv -v filter_log.out      ${exp}/${output_dir}   ;# DART log from filter.f90 (input.nml)
   mv -v filter_server.*.*   ${exp}/${output_dir}   ;# output from filter_server.csh
   mv -v filter.*.log        ${exp}/${output_dir}   ;# output from batch filter.csh

   cp -p filter_server.csh   ${exp}/${output_dir}
   cp -p filter.csh          ${exp}/${output_dir}
   cp -p job.csh             ${exp}/${output_dir}

   echo "completed iteration $i ($OBS_SEQ) at "`date`
# end i loop
   @ i++
end

mv -v run_job.log    ${exp}/run_job_${obs_seq_1}-${obs_seq_n}.log
mv -v namelist     ${exp}
mv -v filter.*.log ${exp}
