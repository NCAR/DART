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
# This script is designed to be submitted as a batch job but may be run from 
# the command line (as a single thread) to check for file motion, etc.
# If running interactively, please comment out the part that actually runs filter.
#
#-----------------------------------------------------------------------------
#
#BSUB -J tiegcm_perfect
#BSUB -o tiegcm_perfect.%J.log
#BSUB -P P3507xxxx
#BSUB -q premium
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
   setenv JOBNAME     tiegcm_perfect
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $HOST

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

setenv CENTRALDIR /glade/scratch/${user}/DART/${JOBNAME}/job_${JOBID}

mkdir -p ${CENTRALDIR}
cd ${CENTRALDIR}

set myname = $0          # this is the name of this script

# some systems don't like the -v option to any of the following 

set OSTYPE = `uname -s` 
switch ( ${OSTYPE} )
   case IRIX64:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      setenv   LINK 'ln -s'
      breaksw
   case AIX:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      setenv   LINK 'ln -s'
      breaksw
   default:
      setenv REMOVE 'rm -rvf'
      setenv   COPY 'cp -v'
      setenv   MOVE 'mv -fv'
      setenv   LINK 'ln -s'
      breaksw
endsw

echo "${JOBNAME} ($JOBID) CENTRALDIR == $CENTRALDIR"

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
# DARTDIR      The location of the DART tiegcm model directory
# TIEGCMDIR    The location of the TIEGCM executable
# ENSEMBLEDIR  The location of the initial ensemble of TIEGCM files
# EXPERIMENT   The (safe) location for the results of this run.
#-----------------------------------------------------------------------------

set     DARTDIR = /glade/u/home/${USER}/DART/tiegcm/models/tiegcm
set   TIEGCMDIR = /glade/p/image/RDA_strawman/TIEGCM_files
set ENSEMBLEDIR = /glade/p/image/RDA_strawman/TIEGCM_files/ensembles
set  EXPERIMENT = /glade/p/work/${USER}/${JOBNAME}

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
# Get the tiegcm executable, control files, and data files.
# The tiegcm initial conditions are in the next block.
#-----------------------------------------------------------------------------

${COPY} ${DARTDIR}/work/perfect_model_obs            . || exit 1
${COPY} ${DARTDIR}/work/dart_to_model                . || exit 1
${COPY} ${DARTDIR}/work/model_to_dart                . || exit 1
${COPY} ${DARTDIR}/work/input.nml   input.nml.original || exit 1
${COPY} ${DARTDIR}/shell_scripts/advance_model.csh   . || exit 1
${COPY} ${DARTDIR}/work/obs_seq.in                   . || exit 1

${COPY} ${TIEGCMDIR}/tiegcm-nompi               tiegcm || exit 1

${COPY} ${ENSEMBLEDIR}/tiegcm_restart_p.nc           . || exit 1
${COPY} ${ENSEMBLEDIR}/tiegcm_s.nc                   . || exit 1
${COPY} ${ENSEMBLEDIR}/tiegcm.nml  tiegcm.nml.original || exit 1

#-----------------------------------------------------------------------------
# Remove all the comments that follow (;) symbol from tiegcm.nml namelist file
# That is a non-standard syntax for fortran namelists.
#
# Ensure that the tiegcm.nml for all the ensemble members is identical
# in all the ways that matter. This will result in a miniumum of changes
# in the advance_model.csh script. This script REQUIRES that there is a  
# SINGLE tiegcm_restart_p.nc. Just keep appending all the timesteps to
# the same file. If you need to subset the large file, use the NCO
# operators. for example    ncks -d time,20,30 tiegcm_restart_p.nc bob.nc 
# If you need more than 300 timesteps in the file, increase it here.
#-----------------------------------------------------------------------------

sed -e 's/;.*//' -e '/^$/ d' \
    -e "/ MXHIST_PRIM /c\ MXHIST_PRIM = 300" \
    -e "/ MXHIST_SECH /c\ MXHIST_SECH = 300" \
    -e "/ SOURCE /c\ SOURCE = 'tiegcm_restart_p.nc'" \
    -e "/ OUTPUT /c\ OUTPUT = 'tiegcm_restart_p.nc'" \
    -e "/ SECOUT /c\ SECOUT = 'tiegcm_s.nc'"         \
    tiegcm.nml.original >! tiegcm.nml  || exit 2

#-----------------------------------------------------------------------------
# Convert a TIEGCM file 'tiegcm_restart.nc' to a DART ics file 'dart_ics'
# There are some requirements for this script and advance_model.csh.
# The requirements for this script are enforced here, the requirements for
# advance_model.csh are enforced there.
# 
# REQUIREMENTS: for input.nml
# model_nml            : tiegcm_restart_file_name   = 'tiegcm_restart_p.nc'
# model_nml            : tiegcm_secondary_file_name = 'tiegcm_s.nc'
# model_nml            : tiegcm_namelist_file_name  = 'tiegcm.nml'
# model_to_dart_nml    : file_out                   = 'dart_ics'
#-----------------------------------------------------------------------------
# perfect_model_obs_nml: async                      = 2
# perfect_model_obs_nml: adv_ens_command            = 'advance_model.csh'
# perfect_model_obs_nml: start_from_restart         = .TRUE.
# perfect_model_obs_nml: restart_in_file_name       = 'dart_ics'
#-----------------------------------------------------------------------------
# dart_to_model_nml    : file_in                    = 'dart_restart'
# dart_to_model_nml    : file_namelist_out          = 'namelist_update'

sed -e "/ tiegcm_restart_file_name /c\ tiegcm_restart_file_name = 'tiegcm_restart_p.nc'" \
    -e "/ tiegcm_secondary_file_name /c\ tiegcm_secondary_file_name = 'tiegcm_s.nc'" \
    -e "/ tiegcm_namelist_file_name /c\ tiegcm_namelist_file_name = 'tiegcm.nml'" \
    -e "/ file_out /c\ file_out = 'dart_ics'" \
    -e "/ async /c\ async = 2" \
    -e "/ adv_ens_command /c\ adv_ens_command = './advance_model.csh'" \
    -e "/ start_from_restart /c\ start_from_restart = .TRUE." \
    -e "/ restart_in_file_name /c\ restart_in_file_name = 'dart_ics'" \
    -e "/ file_in /c\ file_in = 'dart_restart'" \
    -e "/ file_namelist_out /c\ file_namelist_out = 'namelist_update'" \
    input.nml.original >! input.nml  || exit -3

./model_to_dart || exit 2

#-----------------------------------------------------------------------------
# Run perfect_model_obs ... harvest the observations to populate obs_seq.out
# model_mod expects a generic name // advance_model.csh expects a filename
# with the ensemble member ID tacked on - must provide both.
#-----------------------------------------------------------------------------

${REMOVE} tiegcm_restart_p.nc.0001 tiegcm_s.nc.0001 tiegcm.nml.0001 

${LINK} tiegcm_restart_p.nc tiegcm_restart_p.nc.0001 || exit 3
${LINK} tiegcm_s.nc         tiegcm_s.nc.0001         || exit 3
${LINK} tiegcm.nml          tiegcm.nml.0001          || exit 3


./perfect_model_obs || exit 3

#-----------------------------------------------------------------------------
# At this point, all the restart,diagnostic files are in the run/CENTRALDIR.
# You may want to move them to someplace more 'permanent'.
#
# TJH: At this point, the output files have pretty 'generic' names.
# The files could be archived with the assimilation date in their name.
#-----------------------------------------------------------------------------

# ${MOVE} tiegcm_s.nc.0001           ${EXPERIMENT}/perfect/tiegcm_s.nc
# ${MOVE} tiegcm_restart_p.nc.0001   ${EXPERIMENT}/perfect/tiegcm_restart_p.nc
# ${MOVE} tiegcm.nml                 ${EXPERIMENT}/perfect
# ${MOVE} obs_seq.out                ${EXPERIMENT}/perfect
# ${MOVE} True_State.nc              ${EXPERIMENT}/perfect

# ${MOVE} tiegcm_out_1               ${EXPERIMENT}/perfect/tiegcm_out
# ${MOVE} dart_log.out               ${EXPERIMENT}/perfect
# ${MOVE} dart_log.nml               ${EXPERIMENT}/perfect
# Good style dictates that you save the scripts so you can see what worked.

# ${COPY} input.nml                  ${EXPERIMENT}/DART
# ${COPY} *.csh                      ${EXPERIMENT}/DART
# ${COPY} $myname                    ${EXPERIMENT}/DART

echo "${JOBNAME} ($JOBID) finished at "`date`
echo "These are the files in the run directory at completion:"
ls -lrt

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

