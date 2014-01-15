#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Top level script to generate observations and a TRUE state.
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
#
#BSUB -J perfect
#BSUB -o perfect.%J.log
#BSUB -P P3507xxxx
#BSUB -q economy
#BSUB -n 1
#BSUB -W 1:00
#BSUB -N -u ${USER}@ucar.edu

#----------------------------------------------------------------------
# Turns out the scripts are a lot more flexible if you don't rely on 
# the queuing-system-specific variables -- so I am converting them to
# 'generic' names and using the generics throughout the remainder.
#----------------------------------------------------------------------

if ($?LSB_HOSTS) then

   setenv ORIGINALDIR $LS_SUBCWD
   setenv JOBNAME     $LSB_JOBNAME
   setenv JOBID       $LSB_JOBID
   setenv MYQUEUE     $LSB_QUEUE
   setenv MYHOST      $LSB_SUB_HOST

else

   #-------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     tiegcm
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $host

endif

#----------------------------------------------------------------------
# Just an echo of job attributes
#----------------------------------------------------------------------

echo
echo "${JOBNAME} ($JOBID) submitted   from $ORIGINALDIR"
echo "${JOBNAME} ($JOBID) submitted   from $MYHOST"
echo "${JOBNAME} ($JOBID) running in queue $MYQUEUE"
echo "${JOBNAME} ($JOBID) running       on $MYHOST"
echo "${JOBNAME} ($JOBID) started   at "`date`
echo

#----------------------------------------------------------------------
# Make a unique, (empty, clean) temporary directory.
#----------------------------------------------------------------------

setenv TMPDIR /glade/scratch/${user}/DART/${JOBNAME}/job_${JOBID}

mkdir -p ${TMPDIR}
cd ${TMPDIR}

set CENTRALDIR = `pwd`
set myname = $0          # this is the name of this script

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
      setenv   COPY 'cp -v'
      setenv   MOVE 'mv -fv'
      breaksw
endsw

echo "${JOBNAME} ($JOBID) CENTRALDIR == $CENTRALDIR"

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------

set    DARTDIR = /glade/u/home/tmatsuo/DART/models/tiegcm
set  TIEGCMDIR = /glade/u/home/tmatsuo/DART/models/tiegcm/tiegcm_files
set EXPERIMENT = /glade/scratch/tmatsuo/2002_03_28_tiegcm

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
#-----------------------------------------------------------------------------

# executables

 ${COPY} ${DARTDIR}/work/perfect_model_obs          .
 ${COPY} ${DARTDIR}/work/dart_to_model              .
 ${COPY} ${DARTDIR}/work/model_to_dart              .

# shell scripts
 ${COPY} ${DARTDIR}/shell_scripts/advance_model.csh .

# data files
 ${COPY} ${EXPERIMENT}/initial/obs_seq.in           .
 ${COPY} ${DARTDIR}/work/input.nml                  .

#-----------------------------------------------------------------------------
# Get the tiegcm executable, control files, and data files.
#-----------------------------------------------------------------------------

 ${COPY} ${TIEGCMDIR}/tiegcm-nompi                  tiegcm
#${COPY} ${TIEGCMDIR}/tiegcm                        .

 ${COPY} ${EXPERIMENT}/initial/tiegcm_restart_p.nc  .
 ${COPY} ${EXPERIMENT}/initial/tiegcm_s.nc          .
 ${COPY} ${EXPERIMENT}/initial/tiegcm.nml           tiegcm.nml.original

#-----------------------------------------------------------------------------
# Remove all the comments that follow (;) symbol from tiegcm.nml namelist file
#-----------------------------------------------------------------------------

grep -v "^;" tiegcm.nml.original >! tiegcm.nml

#-----------------------------------------------------------------------------
# Check that everything moved OK, and the table is set.
# Convert a TIEGCM file 'tiegcm_restart.nc' to a DART ics file 'perfect_ics'
# 'model_to_dart' has a hardwired output filename of 'temp_ud' ...
#-----------------------------------------------------------------------------

./model_to_dart || exit 1
mv temp_ud perfect_ics

#-----------------------------------------------------------------------------
# Run perfect_model_obs ... harvest the observations to populate obs_seq.out
# model_mod expects a generic name // advance_model.csh expects a filename
# with the ensemble member ID tacked on - must provide both.
#-----------------------------------------------------------------------------

ln -sf tiegcm_restart_p.nc tiegcm_restart_p.nc.0001
ln -sf tiegcm_s.nc         tiegcm_s.nc.0001
ln -sf tiegcm.nml          tiegcm.nml.0001

./perfect_model_obs || exit 2

echo "${JOBNAME} ($JOBID) finished at "`date`

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

exit

${MOVE} tiegcm_s.nc.0001           ${EXPERIMENT}/perfect/tiegcm_s.nc
${MOVE} tiegcm_restart_p.nc.0001   ${EXPERIMENT}/perfect/tiegcm_restart_p.nc
${MOVE} tiegcm.nml                 ${EXPERIMENT}/perfect
${MOVE} obs_seq.out                ${EXPERIMENT}/perfect
${MOVE} True_State.nc              ${EXPERIMENT}/perfect
${MOVE} perfect_restart            ${EXPERIMENT}/perfect

${MOVE} tiegcm_out_1               ${EXPERIMENT}/perfect/tiegcm_out
${MOVE} dart_log.out               ${EXPERIMENT}/perfect
${MOVE} dart_log.nml               ${EXPERIMENT}/perfect
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

