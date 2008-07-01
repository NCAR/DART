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
# job.simple.csh ... Top level script to run a single assimilation experiment.
#
#  Unlike the more complex job.csh, this script only processes a single 
#  observation file.  Still fairly complex; requires a raft of
#  data files and most of them are in hardcoded locations.
#
# You need to know which of several batch systems you are using.  The most
# common one is LSF.   PBS is also common.  (POE is another but is
# not supported directly by this script.  It is not recommended that you have a
# parallel cluster without a batch system (it schedules which nodes are assigned
# to which processes) but it is possible to run that way -- you have to do
# more work to get the information about which nodes are involved to the 
# parallel tasks -- but anyway, there is a section below that uses ssh and no
# batch.
#
# How to submit this job:
#  1. Look at the #BSUB or #PBS sections below and adjust any of the parameters
#     on your cluster.  Queue names are very system specific; some systems 
#     require wall-clock limits; some require an explicit charge code.
#  2. Submit this script to the queue:
#        LSF:   bsub < job.simple.csh
#        PBS:   qsub job.simple.csh
#       NONE:   job.simple.csh
#
# The script moves the necessary files to the current directory and then
# starts 'filter' as a parallel job on all nodes; each of these tasks will 
# call some a separate model_advance.csh when necessary.
#
# The central directory is where the scripts reside and where script and 
# program I/O are expected to happen.
#-----------------------------------------------------------------------------

set CENTRALDIR = `pwd`
set experiment = DARTMIT
alias submit 'bsub < \!*'

set myname = $0     # this is the name of this script

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
echo "Running $experiment on host "`hostname`
echo "Initialized at "`date`
echo "CENTRALDIR is "`pwd`

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------

set DARTDIR = /fs/image/home/${user}/SVN/DART
set DARTMITDIR = ${DARTDIR}/models/MITgcm_ocean
set MITDATADIR = ${DARTMITDIR}/inputs

#-----------------------------------------------------------------------------
# Get the DART/MIT executables and scripts
#-----------------------------------------------------------------------------

${COPY} ${DARTMITDIR}/work/filter                     .
${COPY} ${DARTMITDIR}/work/wakeup_filter              .
${COPY} ${DARTMITDIR}/work/trans_pv_sv                .
${COPY} ${DARTMITDIR}/work/trans_sv_pv                .
${COPY} ${DARTMITDIR}/work/mitgcmuv_20p               .
${COPY} ${DARTMITDIR}/shell_scripts/advance_model.csh .
${COPY} ${DARTMITDIR}/shell_scripts/runme_filter      .

#-----------------------------------------------------------------------------
# Get the necessary data files -- this is the hard part.
# This script does not 'cold start' the ocean model, nor spin up DART.
#-----------------------------------------------------------------------------

if (-d inputs) then
   echo "using existing 'inputs' directory"
else
   echo "Making 'inputs' directory"
   mkdir inputs
endif
${COPY} ${MITDATADIR}/*                               inputs
${COPY} ${DARTMITDIR}/work/obs_seq.out                .
${COPY} ${DARTMITDIR}/work/filter_ics                 .
${COPY} ${DARTMITDIR}/work/input.nml                  .
${COPY} inputs/data                                   .
${COPY} inputs/data.cal                               .

#-----------------------------------------------------------------------------
# Ensure the (output) experiment directory exists
# All the  MIT-related files will get put in ${experiment}/MIT
# All the DART-related files will get put in ${experiment}/DART
#-----------------------------------------------------------------------------

if (-d ${experiment}) then
   echo "${experiment} already exists"
else
   echo "Making run-time directory ${experiment} ..."
   mkdir -p ${experiment}
endif
mkdir -p ${experiment}/{MIT,DART}

#-----------------------------------------------------------------------------
# Runs filter which integrates the results of model advances  (async=4).
#-----------------------------------------------------------------------------

submit runme_filter

echo "Finished at "`date`

exit

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

${MOVE} *.data *.meta         ${experiment}/MIT
${MOVE} STD*                  ${experiment}/MIT

${MOVE} filter_restart*            ${experiment}/DART
${MOVE} assim_model_state_ud[1-9]* ${experiment}/DART
${MOVE} assim_model_state_ic[1-9]* ${experiment}/DART
${MOVE} Posterior_Diag.nc          ${experiment}/DART
${MOVE} Prior_Diag.nc              ${experiment}/DART
${MOVE} obs_seq.final              ${experiment}/DART
${MOVE} dart_log.out               ${experiment}/DART

# Good style dictates that you save the scripts so you can see what worked.

${COPY} input.nml                  ${experiment}/DART
${COPY} *.csh                      ${experiment}/DART
${COPY} $myname                    ${experiment}/DART

ls -lrt
