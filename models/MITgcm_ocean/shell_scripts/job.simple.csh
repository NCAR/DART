#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# job.simple.csh ... Top level script to run a single assimilation experiment.
#
# Unlike the more complex job.csh, this script only processes a single 
# observation file.  Still fairly complex; requires a raft of
# data files and most of them are in hardcoded locations.
#
# This script is designed to be run from the command line (as a single thread)
# and should only take a few seconds to a minute to complete, depending on
# the filesystem performance and data file size.
#
# The script moves the necessary files to the current directory - in DART
# nomenclature, this will be called CENTRALDIR. 
# After everything is confirmed to have been assembled, it is possible
# to edit the data, data.cal, and input.nml files for the specifics of 
# the experiment; as well as allow final configuration of a 'nodelist' file.
#
# Once the 'table is set', all that remains is to start/submit the 
# 'runme_filter' script. That script will spawn 'filter' as a 
# parallel job on the appropriate nodes; each of these tasks will 
# call a separate model_advance.csh when necessary.
#
# The central directory is where the scripts reside and where script and 
# program I/O are expected to happen.
#-----------------------------------------------------------------------------

set CENTRALDIR = `pwd`
set experiment = DARTMIT
alias submit 'bsub < \!*'

set myname = $0     # this is the name of this script
set nonomatch       # suppress "rm" warnings if wildcard does not match anything

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
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
${COPY} ${DARTMITDIR}/work/mitgcmuv_20p               mitgcmuv
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
# Early exit just to check that everything moved OK, and the table is set.
# Gives us a chance to edit the local input.nml, data.cal, etc. if needed.
#-----------------------------------------------------------------------------

exit

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
${MOVE} data data.cal         ${experiment}/MIT
${MOVE} STD*                  ${experiment}/MIT

${MOVE} filter_restart*            ${experiment}/DART
${MOVE} assim_model_state_ud[1-9]* ${experiment}/DART
${MOVE} assim_model_state_ic[1-9]* ${experiment}/DART
${MOVE} analysis.nc                ${experiment}/DART
${MOVE} preassim.nc                ${experiment}/DART
${MOVE} obs_seq.final              ${experiment}/DART
${MOVE} dart_log.out               ${experiment}/DART

# Good style dictates that you save the scripts so you can see what worked.

${COPY} input.nml                  ${experiment}/DART
${COPY} *.csh                      ${experiment}/DART
${COPY} $myname                    ${experiment}/DART

ls -lrt

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

