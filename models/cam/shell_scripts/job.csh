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

#-----------------------------------------------------------------------------
# job.csh ... Script to run whole assimilation experiment. Can easily run for 
# days, given the number of observation sequence files, the size of the model, 
# the number of observations, the number of regions, the number of ensemble
# members. Not to be taken lightly.
#
# Executes 'filter.csh' (locally) and submits 'filter_server.csh' as a batch job 
# for each obs_seq.out file to be processed.
#
# Runs interactively on master node in the central directory or as a batch
# job on a compute node (presuming the compute node has permission to submit
# jobs to the batch queue).
#
# The central directory is where the scripts reside and where script and 
# program I/O are expected to happen.
#-----------------------------------------------------------------------------
#
# BUG; when inflate_diag doesn't exist below 
#      (when do_obs_inflate=F do_single_ss_inflate=F)
#      one of the processes (filter.csh) to kill doesn't exist.
#      Then the bkill fails, and exit is not executed.
#
#=============================================================================
# This block of directives constitutes the preamble for the LSF queuing system 
# LSF is used on the IBM   Linux cluster 'lightning'
# LSF is used on the IMAGe Linux cluster 'coral'
# LSF is used on the IBM   'bluevista'
# The queues on lightning and bluevista are supposed to be similar.
#
# the normal way to submit to the queue is:    bsub < filter_server.csh
#
# an explanation of the most common directives follows:
# -J Job name (master script job.csh presumes filter_server.xxxx.log)
# -o STDOUT filename
# -e STDERR filename
# -P      account
# -q queue    cheapest == [standby, economy, (regular,debug), premium] == $$$$
# -n number of processors  (really)
##=============================================================================
#BSUB -J DARTCAM
#BSUB -o DARTCAM.%J.log
#BSUB -q standby
#BSUB -n 1
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
#PBS -N DARTCAM
#PBS -r n
#PBS -e DARTCAM.err
#PBS -o DARTCAM.log
#PBS -q verylong
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
   set SCRATCHDIR = /ptmp/${user}/filter_server
   alias submit 'bsub < \!*'
   
else if ($?PBS_O_WORKDIR) then

   # PBS has a list of processors in a file whose name is (PBS_NODEFILE)

   set CENTRALDIR = $PBS_O_WORKDIR
   set JOBNAME = $PBS_JOBNAME
   set PROCNAMES = `cat $PBS_NODEFILE`
   set REMOTECMD = rsh
   set SCRATCHDIR = /scratch/local/${user}/filter_server
   alias submit 'qsub \!*'

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
   set JOBNAME = interactive_filter_server
   set PROCNAMES = `cat $OCOTILLO_NODEFILE`
   set REMOTECMD = rsh
   set SCRATCHDIR = /var/tmp/${user}/filter_server
   alias submit 'rsh \!*'
   
else

   # interactive
   # you never really want to run filter_server.csh for a big job
   # interactively, so I am aliasing the 'submit' command to
   # use the LSF queueing system. Your mileage may vary.

   set CENTRALDIR = `pwd`
   set JOBNAME = interactive_filter_server
   set PROCNAMES = "$host $host $host $host"
   set REMOTECMD = csh
   set SCRATCHDIR = /tmp/${user}/filter_server
   alias submit 'bsub < \!*'
   
endif

set myname = $0     # this is the name of this script

setenv REMOVE 'rm -rvf'
setenv   COPY 'cp -vp'
setenv   MOVE 'mv -vf'

cd ${CENTRALDIR}

# Determine number of processors
set NPROCS = `echo $PROCNAMES | wc -w`

# Set Variable for a 'master' logfile 
set MASTERLOG = ${CENTRALDIR}/run_job.log

echo " "                                             >! $MASTERLOG
echo "Running $JOBNAME on host "`hostname`           >> $MASTERLOG
echo "Initialized at "`date`                         >> $MASTERLOG
echo "CENTRALDIR is "`pwd`                           >> $MASTERLOG
echo "This job has allocated $NPROCS processors."    >> $MASTERLOG
echo "they are: "                                    >> $MASTERLOG
echo $PROCNAMES                                      >> $MASTERLOG

echo " "
echo "Running $JOBNAME on host "`hostname`
echo "Initialized at "`date`
echo "CENTRALDIR is "`pwd`
echo "This job has allocated $NPROCS processors."
echo "they are: "
echo $PROCNAMES

# Run parameters to change

# Directory where output will be kept (relative to '.')
set exp = Experiment1

# 'day'/obs_seq.out numbers to assimilate during this job
# First and last obs seq files. 
set obs_seq_1 = 1
set obs_seq_n = 1

# The month of these obs_seq.out files (pad with 0s; they're used both as numbers,
# which does work, and as char strings in names.
set mo = 01
set mo_first = 01

# number of obs_seqs / day
set obs_seq_freq = 2

# The "day" of the first obs_seq.out file of the experiment
set obs_seq_first = 1
set obs_seq_root = /ptmp/thoar/T21/obs_seq2003

# Subdirectory root name where each "day"'s output will be kept
# Obs_diag restricts these subdirectory names to be 5 characters.
# output_root = xxx    xxx signifies the month OF THE OBS_SEQ_FIRST.
set output_root = 01_

set num_ens = 40


# T21
# inflate_1_ic is the wrong size, but I need to set it to something
# The CAMsrc directory is MORE than just the location of the executable.
# There are more support widgets expected in the directory tree.
set inflate_1_ic = ../Pre-J/Exp4/01_62/DART
set obs_seq_1_ic  = /ptmp/raeder/CAM_init/T21x80/03-01-01/DART_lunes
set obs_seq_1_cam = /ptmp/raeder/CAM_init/T21x80/03-01-01/CAM/caminput_
set obs_seq_1_clm = /ptmp/raeder/CAM_init/T21x80/03-01-01/CLM/clminput_
set CAMsrc = /home/coral/raeder/Cam3/cam3.1/models/atm/cam/bld/T21-O2

${REMOVE} caminput.nc clminput.nc namelistin
ln -s caminput_T85.nc caminput.nc
ln -s clminput_T85.nc clminput.nc
ln -s namelistin_T85 namelistin

# Each obs_seq file gets its own input.nml ... 
# the following variable helps set this up. 

set input = input_


# CHANGE choice of whether forecasts on caminput_##.nc files should be replaced with
#        analyses from filter_ic_new.#### at the end of each obs_seq.  See 'stuff' below

# CHANGE names of obs_seq files below, if necessary


# Not currently used; all previous days initial files are gzipped
# Now don't even gzip those; keep available for other restarts
set save_freq = 2
set mod_save = 1

set days_in_mo = (31 28 31 30 31 30 31 31 30 31 30 31)
# leap years (but year not defined here);   
#    if (($year % 4) == 0) @ days_in_mo[2] = $days_in_mo[2] + 1
# leap year every 4 years except for century marks, but include centuries divisible by 400
#    So, all modern years divisible by 4 are leap years.

#----------------------------------------------------------
echo "exp num_ens obs_seq_1 obs_seq_n obs_seq_first"
echo "$exp $num_ens $obs_seq_1 $obs_seq_n $obs_seq_first"
echo "obs_seq_1_ic is $obs_seq_1_ic"

echo "exp num_ens obs_seq_1 obs_seq_n obs_seq_first"       >> $MASTERLOG
echo "$exp $num_ens $obs_seq_1 $obs_seq_n $obs_seq_first"  >> $MASTERLOG
echo "obs_seq_1_ic is $obs_seq_1_ic"                       >> $MASTERLOG



# Ensure the experiment directory exists

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
   ${COPY} namelistin ${exp}/namelistin     # just for posterity
endif

# clean up old CAM inputs that may be laying around

if (-e caminput_1.nc) then
   ${REMOVE} clminput_[1-9]*.nc 
   ${REMOVE} caminput_[1-9]*.nc 
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

   if ($i == $obs_seq_first) then
      ${COPY} ${input}${obs_seq_first}.nml input.nml
   else
      if (  -e    ${input}n.nml) then
         ${COPY}  ${input}n.nml    input.nml
      else if (-e ${input}${i}.nml) then
         ${COPY}  ${input}${i}.nml input.nml
      else
         echo "input_next is MISSING" >> $MASTERLOG
         echo "input_next is MISSING"
         exit
      endif
   endif

   
   #----------------------------------------------------------------------
   # Get obs_seq file for this assimilation based on date.
   # Since observation sequence file (names) kinda sorta reflect something
   # that _might_ be a date, we find one we want and link it to this directory.
   # As such, we need to ensure any existing obs_seq.out is GONE, which is
   # a little scary.
   #
   # set day and hour for this obs_seq
   # UGLY; need more general calendar capability in here.
   # Normally, when starting with 12 hourly from beginning;  
   #----------------------------------------------------------------------

   ${REMOVE} obs_seq.out

   set day_this_mo = $i
   @ month = $mo - 1
   while ($month >= $mo_first)
       @ day_this_mo = $day_this_mo - $days_in_mo[$month] * $obs_seq_freq
       @ month = $month - 1
   end
   @ day = ($day_this_mo + 1) / $obs_seq_freq
   @ hour = ($i % $obs_seq_freq) * 24 / $obs_seq_freq 
   if ($hour == 0) set hour = 24
   echo "obs_seq, day, hour = $i $day $hour " >> $MASTERLOG

   if ($i == 0) then
      set OBS_SEQ = /scratch/cluster/raeder/GWD_T42/obs_seq_1-1-03.out
   else if ($i < 19) then
      set OBS_SEQ = /ncar/dart/Obs_sets/Allx12/obs_seq2003010${day}${hour}
   else
      set OBS_SEQ = /ncar/dart/Obs_sets/Allx12/obs_seq200301${day}${hour}
      ## obs_seq_freq = 24
      ## set OBS_SEQ = /ncar/dart/Obs_sets/All/obs_seq200301${i}
   endif

   if (  -s $OBS_SEQ ) then    ;# -s is true if xxxx has non-zero size
      ln -s $OBS_SEQ  obs_seq.out
   else
      echo "ERROR - no obs_seq.out for $i - looking for:" 
      echo $OBS_SEQ 
      exit 123
   endif

   echo "job- obs_seq $i used is $OBS_SEQ" >> $MASTERLOG
   echo "job- obs_seq $i used is $OBS_SEQ"
   ls -lt obs_seq.out

   #----------------------------------------------------------------------
   # Get initial conditions for DART, model from 'permanent' storage or
   # from the result of a previous experiment.
   #
   # from_root defines location of DART initial condition(s)
   # cam_init, clm_init are where CAM gets it's initial files (in advance_model.csh)
   #----------------------------------------------------------------------

   if ($i == $obs_seq_first) then 
      # get 'initial' initial files
      set from_root = ${obs_seq_1_ic} 
      set  cam_init = ${obs_seq_1_cam}
      set  clm_init = ${obs_seq_1_clm}
   else
      # get initial files from result of previous experiment.
      set from_root = `pwd`/$exp/${out_prev}/DART
      set cam_init =  `pwd`/$exp/${out_prev}/CAM/caminput_
      set clm_init =  `pwd`/$exp/${out_prev}/CLM/clminput_
   endif

   # transmit info to advance_model.csh, run by filter_server.csh
   # The second item echoed must be the subdirectory where CAM is kept, in the Central directory
   echo "$exp $CAMsrc $cam_init $clm_init" >! casemodel


   ${REMOVE} inflate_ic_old
   if ($i == $obs_seq_first) then 
      ln -s $inflate_1_ic/inflate_ics inflate_ic_old
   else
      ln -s    $from_root/inflate_ics inflate_ic_old
   endif

   # link to filter_ic file(s), so that filter.csh can copy them to a compute node
   if (-e ${from_root}/filter_ic.0001) then
      set n = 1
      while($n <= ${num_ens})
           set from = ${from_root}/filter_ic*[.0]$n
           ${REMOVE}   filter_ic_old.$from:e
           ln -s $from filter_ic_old.$from:e
           @ n++ 
      end
   else if (-e ${from_root}/filter_ic) then
      ${REMOVE} filter_ic_old
      ln  -s   ${from_root}/filter_ic filter_ic_old
   endif
   echo ' '
   echo "job- filter_ic_old is/are" >> $MASTERLOG
   echo "job- filter_ic_old is/are"
   ls -lt filter_ic_old*            >> $MASTERLOG
   ls -lt filter_ic_old*

   # Make subdirectories to store the output from the assimilation of this
   # observation sequence file. 

   mkdir -p ${exp}/${output_dir}
   mkdir -p ${exp}/${output_dir}/{CLM,CAM,DART}

   #======================================================================
   # Run the filter in async=3 mode.
   # This is the central directory for whole filter job
   #    qsub jobs submitted from there
   #    semaphor files must (dis)appear there, 
   #    I/O between filter and advance_model and assim_region goes through there.
   #    Final output is put there
   # It's CENTRALDIR  in filter.csh and filter_server.csh

   # advances model and assims regions (do this first to grab whole nodes)
   # We need to capture the batch job number to kill later if need be.
   #======================================================================

   submit filter_server.csh > batchsubmit$$
   set STRING = "1,$ s#<##g"
   sed -e "$STRING" batchsubmit$$ > bill$$
   set STRING = "1,$ s#>##g"
   sed -e "$STRING" bill$$ > batchsubmit$$
   set STRING = `cat batchsubmit$$`
   set FILTERSERVERBATCHID = $STRING[2]
   ${REMOVE} batchsubmit$$ bill$$

   # runs filter, which integrates the results of model advances and region assims
   # This only uses 1 processor, so it could fit on a node with someone else.
   submit filter.csh  > batchsubmit$$
   set STRING = "1,$ s#<##g"
   sed -e "$STRING" batchsubmit$$ > bill$$
   set STRING = "1,$ s#>##g"
   sed -e "$STRING" bill$$ > batchsubmit$$
   set STRING = `cat batchsubmit$$`
   set FILTERBATCHID = $STRING[2]
   ${REMOVE} batchsubmit$$ bill$$

   echo "filter_server spawned as job $FILTERSERVERBATCHID at "`date`
   echo "filter        spawned as job $FILTERBATCHID at "`date`
   echo "filter_server spawned as job $FILTERSERVERBATCHID at "`date` >> $MASTERLOG
   echo "filter        spawned as job $FILTERBATCHID at "`date`       >> $MASTERLOG

   set KILLCOMMAND = "bkill $FILTERSERVERBATCHID $FILTERBATCHID; touch BOMBED; exit"

# BUG; when inflate_diag doesn't exist below, one of the processes to kill doesn't exist.
#      Then the bkill fails, and exit is not executed.

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

         ${REMOVE} go_end_filter   ; # do this when go_end_filter first appears
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
   # ${MOVE} clminput_[1-9]*.nc    ${exp}/${output_dir}/CLM
   # ${MOVE} caminput_[1-9]*.nc    ${exp}/${output_dir}/CAM
   #-----------------------------------------------

   echo "Listing contents of CENTRALDIR before archiving at "`date`
   ls -l
   echo "Listing contents of CENTRALDIR before archiving at "`date` >> $MASTERLOG
   ls -l >> $MASTERLOG

   foreach FILE ( Prior_Diag.nc Posterior_Diag.nc obs_seq.final )
      if ( -s $FILE ) then
         ${MOVE} $FILE ${exp}/${output_dir} 
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

   # inflate_diag may or may not exist (when do_obs_inflate=F do_single_ss_inflate=F)
   # so don't die if it's missing.  We'd have to query input.nml to learn if it should exist.

   if ( -s inflate_diag ) then
      ${MOVE} inflate_diag ${exp}/${output_dir} 
      if ( ! $status == 0 ) then
         echo "failed moving ${CENTRALDIR}/inflate_diag" >>$MASTERLOG
         echo "failed moving ${CENTRALDIR}/inflate_diag"
         $KILLCOMMAND
      endif
   endif

   # Move the filter restart file(s) to the storage subdirectory
   echo "moving filter_ic_newS to ${exp}/${output_dir}/DART/filter_icS"

   if (-e filter_ic_new) then
      ${MOVE} filter_ic_new ${exp}/${output_dir}/DART/filter_ic
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
           ${MOVE} $from $dest
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

#   if (-e assim_ic_new) then
#      ${MOVE} assim_ic_new ${exp}/${output_dir}/DART/assim_tools_ics
#      if (! $status == 0 ) then
#         echo "failed moving assim_ic_new to ${exp}/${output_dir}/DART/assim_tools_ics"
#         $KILLCOMMAND
#      endif
#   endif

   if (-e inflate_ic_new) then
      ${MOVE} inflate_ic_new ${exp}/${output_dir}/DART/inflate_ics
      if (! $status == 0 ) then
         echo "failed moving inflate_ic_new to ${exp}/${output_dir}/DART/inflate_ics"
         $KILLCOMMAND
      endif
   endif

   set n = 1
   while($n <= ${num_ens})    ;# loop over all ensemble members ... one-at-a-time
      set CAMINPUT = caminput_${n}.nc 
      set CLMINPUT = clminput_${n}.nc 

      if (  -s   $CAMINPUT ) then
         ${MOVE} $CAMINPUT ${exp}/${output_dir}/CAM 
         if (! $status == 0 ) then
            echo "failed moving ${CENTRALDIR}/$CAMINPUT"
            $KILLCOMMAND
         endif
      else
         echo "... aieeee ... ${CENTRALDIR}/$CAMINPUT does not exist and should."
         $KILLCOMMAND
      endif

      if (  -s   $CLMINPUT ) then
         ${MOVE} $CLMINPUT ${exp}/${output_dir}/CLM 
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
       echo okay >! rm_filter_temp
   else
#      This causes loop/job to exit gracefully
       @ i = $obs_seq_n + 1
       echo "RETRIEVE filter_ic files from filter.csh temp directory"
       echo "Then remove temp and cam advance temps"
       echo "RETRIEVE filter_ic files from filter.csh temp directory" >> $MASTERLOG
       echo "Then remove temp and cam advance temps" >> $MASTERLOG
   endif

   # Compress and archive older output which we won't need immediately
   # This relies on having auto_re2ms_LSF.csh in your $HOME directory ... TJH
   # if ($i % $save_freq != $mod_save) then
   if ($i > $obs_seq_first ) then
      cd ${exp}/${out_prev}
      eval ${SUBMIT} ~/auto_re2ms_LSF.csh                                 >>& $MASTERLOG
      cd ../..
      echo "Backing up restart $j to mass store;  in separate batch job"  >> $MASTERLOG
# Try this for Ave;  works
        echo "Backing up restart $j to mass store;  in separate batch job"  >> $MASTERLOG
   else
        rm CAM/* CLM/* DART/filter*
   endif
   cd ../..
endif



   ${MOVE} cam_out_temp1       ${exp}/${output_dir}   ;# save a representative model advance
   ${MOVE} cam_reg_temp1       ${exp}/${output_dir}   ;# ditto for assimilation region
   ${MOVE} input.nml           ${exp}/${output_dir}
   ${MOVE} casemodel           ${exp}/${output_dir}
   ${MOVE} run_filter.stout    ${exp}/${output_dir}   ;# stdout   from filter.f90
   ${MOVE} filter_log.out      ${exp}/${output_dir}   ;# DART log from filter.f90 (input.nml)
   ${MOVE} filter_server.*.*   ${exp}/${output_dir}   ;# output from filter_server.csh
   ${MOVE} filter.*.log        ${exp}/${output_dir}   ;# output from batch filter.csh

   ${COPY} filter_server.csh   ${exp}/${output_dir}
   ${COPY} filter.csh          ${exp}/${output_dir}
   ${COPY} job.csh             ${exp}/${output_dir}

   echo "completed iteration $i ($OBS_SEQ) at "`date`
   @ i++
end # end of the huge "i" loop

${MOVE} run_job.log  ${exp}/run_job_${obs_seq_1}-${obs_seq_n}.log
${MOVE} namelist     ${exp}
${MOVE} filter.*.log ${exp}
${REMOVE} ~/lnd.*.rpointer

