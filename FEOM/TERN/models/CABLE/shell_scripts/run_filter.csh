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
#BSUB -J filter
#BSUB -o filter.%J.log
#BSUB -P P3507xxxx
#BSUB -q economy
#BSUB -n 64
#BSUB -R "span[ptile=64]"
#BSUB -W 3:00
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

 ${COPY} ${DARTDIR}/work/filter                     .
 ${COPY} ${DARTDIR}/work/dart_to_model              .
 ${COPY} ${DARTDIR}/work/model_to_dart              .

# shell scripts
 ${COPY} ${DARTDIR}/shell_scripts/advance_model.csh .

# data files
 ${COPY} ${EXPERIMENT}/initial/obs_seq.out          .
 ${COPY} ${DARTDIR}/work/input.nml                  .

#-----------------------------------------------------------------------------
# Get the tiegcm executable, control files, and data files.
#-----------------------------------------------------------------------------

 ${COPY} ${TIEGCMDIR}/tiegcm-nompi                  tiegcm
#${COPY} ${TIEGCMDIR}/tiegcm                        .

#-----------------------------------------------------------------------------
# Get the tiegcm input state ... for this experiment, we generated the ensemble by: 
#
# ${COPY} ${TIEGCMDIR}/tiegcm_s.nc                   .
# ${COPY} ${TIEGCMDIR}/tiegcm_restart_p.nc           .
# ./model_to_dart || exit 1
# mv temp_ud filter_ics
#
# REQUIREMENT for the case where we have an initial ensemble:
# input.nml:filter_nml:start_from_restart = .TRUE.
# input.nml:filter_nml:restart_in_file    = 'filter_ics'
#-----------------------------------------------------------------------------
# Put all of the DART initial conditions files and all of the TIEGCM files
# in the CENTRALDIR - preserving the ensemble member ID for each filename.
# The advance_model.csh script will copy the appropriate files for each 
# ensemble member into the model advance directory.
# These files may be linked to CENTRALDIR since they get copied to the
# model advance directory. 
#-----------------------------------------------------------------------------

set ENSEMBLESTRING = `/usr/local/bin/grep -A 42 filter_nml input.nml | grep ens_size`
set NUM_ENS = `echo $ENSEMBLESTRING[3] | sed -e "s#,##"`

@ i = 1
while ( $i <= $NUM_ENS )

  set darticname  = `printf "filter_ics.%04d"          $i`
  set tiesecond   = `printf "tiegcm_s.nc.%04d"         $i`
  set tierestart  = `printf "tiegcm_restart_p.nc.%04d" $i`
  set tieinp      = `printf "tiegcm.nml.%04d"          $i`

  ln -sf ${EXPERIMENT}/initial/$darticname .
  ln -sf ${EXPERIMENT}/initial/$tiesecond  .
  ln -sf ${EXPERIMENT}/initial/$tierestart .
  ln -sf ${EXPERIMENT}/initial/$tieinp     tiegcm.nml.original 

  sed -e 's/;.*//' -e '/^$/ d' tiegcm.nml.original >! $tieinp

  @ i += 1
end

#-----------------------------------------------------------------------------
# Run filter ... 
#-----------------------------------------------------------------------------

ln -sf tiegcm_restart_p.nc.0001 tiegcm_restart_p.nc
ln -sf tiegcm_s.nc.0001         tiegcm_s.nc
ln -sf tiegcm.nml.0001          tiegcm.nml

mpirun.lsf ./filter || exit 2

echo "${JOBNAME} ($JOBID) finished at "`date`

#-----------------------------------------------------------------------------
# Move the output to storage after filter completes.
# At this point, all the DART restart,diagnostic files are in the CENTRALDIR
# and need to be moved to the 'experiment permanent' directory.
#
# TJH: At this point, the output files have pretty 'generic' names.
# The files should be archived with the assimilation date in their name.
#-----------------------------------------------------------------------------

exit

${MOVE} tiegcm_s.nc*               ${experiment}/tiegcm
${MOVE} tiegcm_restart_p.nc*       ${experiment}/tiegcm
${MOVE} tiegcm_out_*               ${experiment}/tiegcm

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

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

