#!/bin/csh

# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$

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
#BSUB -q debug
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
   alias submit 'bsub < \!*'
   
else if ($?PBS_O_WORKDIR) then

   # PBS has a list of processors in a file whose name is (PBS_NODEFILE)

   set CENTRALDIR = $PBS_O_WORKDIR
   set JOBNAME = $PBS_JOBNAME
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
   set JOBNAME = DARTCAM
   alias submit 'csh \!*'
   
else

   # interactive
   # YOU need to know if you are using the PBS or LSF queuing
   # system ... and set 'submit' accordingly.

   set CENTRALDIR = `pwd`
   set JOBNAME = DARTCAM
   alias submit 'bsub < \!*'
   
endif

set myname = $0     # this is the name of this script

# Set the experiment name.

set experiment = CAM1X

cd ${CENTRALDIR}

# some systems don't like the -v option to any of the following 

set OSTYPE = `uname -s` 
switch ( ${OSTYPE} )
   case IRIX64:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      breaksw
   case AIX:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      breaksw
   default:
      setenv REMOVE 'rm -rvf'
      setenv   COPY 'cp -vp'
      setenv   MOVE 'mv -fv'
      breaksw
endsw

echo " "
echo "Running $JOBNAME on host "`hostname`
echo "Initialized at "`date`
echo "CENTRALDIR is "`pwd`

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------

set DARTDIR = /image/home/${user}/DART
set DARTCAMDIR = ${DARTDIR}/models/cam
set CAMDATADIR = /image/home/${user}/CAMDATA

#-----------------------------------------------------------------------------
# Get the queueing-system (independent?) scripts
#-----------------------------------------------------------------------------

${COPY} ${DARTDIR}/shell_scripts/advance_ens.csh      .
${COPY} ${DARTDIR}/shell_scripts/assim_filter.csh     .
${COPY} ${DARTDIR}/shell_scripts/filter_server.csh    .

#-----------------------------------------------------------------------------
# Get the model-specific scripts
#-----------------------------------------------------------------------------

${COPY} ${DARTCAMDIR}/shell_scripts/assim_region.csh  .
${COPY} ${DARTCAMDIR}/shell_scripts/advance_model.csh .

#-----------------------------------------------------------------------------
# Get the DARTCAM executables
#-----------------------------------------------------------------------------

${COPY} ${DARTCAMDIR}/work/filter                     .
${COPY} ${DARTCAMDIR}/work/assim_region               .
${COPY} ${DARTCAMDIR}/work/trans_date_to_dart         .
${COPY} ${DARTCAMDIR}/work/trans_pv_sv                .
${COPY} ${DARTCAMDIR}/work/trans_pv_sv_time0          .
${COPY} ${DARTCAMDIR}/work/trans_sv_pv                .
${COPY} ${DARTCAMDIR}/work/trans_time                 .

#-----------------------------------------------------------------------------
# Get the necessary data files -- this is the hard part.
# This script does not involve 'cold starting' CAM, nor spinning up DART.
# The DARTics directory has one initial conditions file for
# each ensemble member. We need one for each ...
# The input.nml has a restart_in_file_name of 'filter_ic_old'
# which must match the filename here. 
# Because that same namelist has 'single_restart_file_in' as .false.,
# the restart_in_file_name gets an ensemble member number appended to it.
#-----------------------------------------------------------------------------

${COPY} ${CAMDATADIR}/input.nml                       .
${COPY} ${CAMDATADIR}/obs_seq.out                     .

# try to discover the ensemble size from the input.nml
# this is some gory shell programming ... all to do 'something simple'

grep ens_size input.nml >! ensstring.$$
set  STRING = "1,$ s#,##g"
set ensstring = `sed -e "$STRING" ensstring.$$`
set num_ens = $ensstring[3]

${REMOVE} ensstring.$$

echo "There are ${num_ens} ensemble members."

# This just copies just the initial conditions for the correct number
# of ensemble members.

set DARTics  = /ptmp/raeder/CAM_init/T21x80/03-01-01/DART_lunes

set n = 1
while($n <= ${num_ens})
   set from = ${DARTics}/filter_ic*[.0]$n
   ${COPY} $from filter_ic_old.$from:e
   @ n++ 
end

${COPY} ${CAMDATADIR}/namelistin                      .
${COPY} ${CAMDATADIR}/caminput.nc                     .
${COPY} ${CAMDATADIR}/clminput.nc                     .
set CAMics = /ptmp/raeder/CAM_init/T21x80/03-01-01/CAM/caminput_
set CLMics = /ptmp/raeder/CAM_init/T21x80/03-01-01/CLM/clminput_

#-----------------------------------------------------------------------------
# T21
# inflate_1_ic is the wrong size, but I need to set it to something
# The CAMsrc directory is MORE than just the location of the executable.
# There are more support widgets expected in the directory tree.
#-----------------------------------------------------------------------------

# set inflate_1_ic = ../Pre-J/Exp4/01_62/DART

set CAMsrc = /home/coral/raeder/Cam3/cam3.1/models/atm/cam/bld/T21-O2

#-----------------------------------------------------------------------------
# Ensure the (output) experiment directory exists
# All the  CAM-related files will get put in ${experiment}/CAM
# All the  CLM-related files will get put in ${experiment}/CLM
# All the DART-related files will get put in ${experiment}/DART
#-----------------------------------------------------------------------------

if (-d ${experiment}) then
   echo "${experiment} already exists"
else
   echo "Making run-time directory ${experiment} ..."
   mkdir -p ${experiment}
endif
mkdir -p ${experiment}/{CLM,CAM,DART}

#-----------------------------------------------------------------------------
# This is where I should check to make sure all the required files exist.
#-----------------------------------------------------------------------------

if (! -e namelistin ) then
   echo "ERROR ... need a namelistin file."
   exit 99
endif

#-----------------------------------------------------------------------------
# Some information about CAM must be made available to advance_model.csh
# filter_server.csh spawns advance_ens.csh which spawns advance_model.csh
# 'casemodel' is required (by advance_model.csh) to be in the Central directory
#-----------------------------------------------------------------------------

echo "${experiment} ${CAMsrc} ${CAMics} ${CLMics}" >! casemodel

#-----------------------------------------------------------------------------
# Run the filter in async=3 mode.
# This is the central directory for whole filter job
#    qsub jobs submitted from there
#    semaphor files must (dis)appear there, 
#    I/O between filter and advance_model and assim_region goes through there.
#    Final output is put there
# It's CENTRALDIR  in filter.csh and filter_server.csh
# advances model and assims regions (do this first to grab whole nodes)
# We need to capture the batch job number to kill later if need be.
#-----------------------------------------------------------------------------
# A 20 member ensemble @ T21 can take anywhere between 10-30 minutes.
#-----------------------------------------------------------------------------

submit filter_server.csh

#-----------------------------------------------------------------------------
# Runs filter - IN THE FOREGROUND - which integrates the results of model 
# advances and region assimilations.
# This only uses 1 processor, so it easily runs on this node.
#-----------------------------------------------------------------------------

./filter || exit 23

#-----------------------------------------------------------------------------
# When filter.f90 finished, it creates a file  called 'go_end_filter' in this
# runtime directory (i.e. CENTRALDIR). The existence of 'go_end_filter' is
# enough to signal filter_server.csh 
#
# time to end. 
# filter_server.csh is in an infinite loop looking for the existence
# of any of three files: 
#      go_advance_model   (time to advance the ensemble members)
#      go_assim_regions   (time to assimilate the observations)
#      go_end_filter      (time to end)
#-----------------------------------------------------------------------------

echo "Finished at "`date`

#-----------------------------------------------------------------------------
# Move the output to storage after filter completes.
# At this point, all the restart,diagnostic files are in the CENTRALDIR
# and need to be moved to the 'experiment permanent' directory.
# We have had problems with some, but not all, files being moved
# correctly, so we are adding bulletproofing to check to ensure the filesystem
# has completed writing the files, etc. Sometimes we get here before
# all the files have finished being written.
#-----------------------------------------------------------------------------

echo "Listing contents of CENTRALDIR before archiving"
ls -l

${MOVE} clminput_[1-9]*.nc         ${experiment}/CLM

${MOVE} cam_out_temp[1-9]*         ${experiment}/CAM
${MOVE} caminput_[1-9]*.nc         ${experiment}/CAM

${MOVE} filter_ic_old*             ${experiment}/DART
${MOVE} filter_ic_new*             ${experiment}/DART
${MOVE} assim_model_state_ud[1-9]* ${experiment}/DART
${MOVE} assim_model_state_ic[1-9]* ${experiment}/DART
${MOVE} inflate_ic_new             ${experiment}/DART
${MOVE} filter_control             ${experiment}/DART
${MOVE} Posterior_Diag.nc          ${experiment}/DART
${MOVE} Prior_Diag.nc              ${experiment}/DART
${MOVE} obs_seq.final              ${experiment}/DART
${MOVE} dart_log.out               ${experiment}/DART
${MOVE} run_job.log                ${experiment}/DART   # filter_server runtime log

${COPY} namelistin                 ${experiment}
${MOVE} namelist                   ${experiment}
${MOVE} casemodel                  ${experiment}

# Good style dictates that you save the scripts so you can see what worked.

${COPY} input.nml                  ${experiment}
${COPY} *.csh                      ${experiment}
${COPY} $myname                    ${experiment}

# CAM leaves a bunch of remnants in your $HOME directory.
# I have not figured out how to use them ... so I clean up.

${REMOVE} ~/lnd.*.rpointer

ls -lrt

echo "Depending on when filter_server.csh finishes, you may wind up"
echo "with a couple files called filter_server.xxxx.[log,err]"
echo "You could/should move them to ${experiment}/DART"
echo "filter_server.csh will also remove the semaphor file go_end_filter,"
echo "so do not remove it. If it still exists after filter_server has completed"
echo "something is wrong ..."
echo "Cheers."
