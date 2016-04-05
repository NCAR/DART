#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Top level script to perform an assimilation.
#
# This script is designed to be submitted as a batch job but may be run from 
# the command line (as a single thread) to check for file motion, etc.
# If running interactively, please comment out the part that actually runs filter.
#
# PLEASE READ THE FOLLOWING: 
#    Setting the number of tasks and choosing the right ptile requires work.
# The number of tasks (-n) can be be as big as the ensemble size for
# a single-threaded tiegcm (i.e. async == 2) so that all ensemble members can 
# run simultaneously. The setting of ptile specifies the number of tasks on each 
# node, which usually depends on the model resolution and subsequent memory use 
# of each ensemble member. Think of ptile as the number of ensemble members you 
# can run on one node and not run out of the shared memory on that node.
#    If you specify more tasks than ensemble members, there are tasks that have
# nothing to do during the model advance. If the model advance step takes longer
# than the MPI timeout on your machine, you may need to disable the MPI timeout.
#-----------------------------------------------------------------------------
#
#BSUB -J tiegcm_filter
#BSUB -o tiegcm_filter.%J.log
#BSUB -P P3507xxxx
#BSUB -q regular
#BSUB -n 80
#BSUB -R "span[ptile=16]"
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
   setenv MPI_RUN_CMD mpirun.lsf

   # MP_DEBUG_NOTIMEOUT may alleviate MPI timeouts that may occur under
   # certain task geometries. It is NOT a good idea to use it in general. 
   # setenv MP_DEBUG_NOTIMEOUT yes

else

   #-------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     tiegcm_filter
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $HOST
   setenv MPI_RUN_CMD ''

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

${COPY} ${DARTDIR}/work/filter                       . || exit 1
${COPY} ${DARTDIR}/work/dart_to_model                . || exit 1
${COPY} ${DARTDIR}/work/model_to_dart                . || exit 1
${COPY} ${DARTDIR}/work/input.nml   input.nml.original || exit 1
${COPY} ${DARTDIR}/shell_scripts/advance_model.csh   . || exit 1
${COPY} ${EXPERIMENT}/observation/obs_seq.out        . || exit 1

${COPY}  ${TIEGCMDIR}/tiegcm-nompi              tiegcm || exit 1

#-----------------------------------------------------------------------------
# Put all of the DART initial conditions files and all of the TIEGCM files
# in the CENTRALDIR - preserving the ensemble member ID for each filename.
# The advance_model.csh script will copy the appropriate files for each 
# ensemble member into the model advance directory.
# These files may be linked to CENTRALDIR since they get copied to the
# model advance directory. 
#
# REQUIREMENTS: for input.nml
# model_nml            : tiegcm_restart_file_name   = 'tiegcm_restart_p.nc'
# model_nml            : tiegcm_secondary_file_name = 'tiegcm_s.nc'
# model_nml            : tiegcm_namelist_file_name  = 'tiegcm.nml'
# model_to_dart_nml    : file_out                   = 'dart_ics'
#-----------------------------------------------------------------------------
# ensemble_manager_nml : single_restart_file_in     = .false.
# filter_nml           : async                      = 2
# filter_nml           : adv_ens_command            = './advance_model.csh'
# filter_nml           : start_from_restart         = .TRUE.
# filter_nml           : restart_in_file_name       = 'filter_ics'
# filter_nml           : restart_out_file_name      = 'filter_restart'
#-----------------------------------------------------------------------------
# dart_to_model_nml    : file_in                    = 'dart_restart'
# dart_to_model_nml    : file_namelist_out          = 'namelist_update'

sed -e "/ tiegcm_restart_file_name /c\ tiegcm_restart_file_name = 'tiegcm_restart_p.nc'" \
    -e "/ tiegcm_secondary_file_name /c\ tiegcm_secondary_file_name = 'tiegcm_s.nc'" \
    -e "/ tiegcm_namelist_file_name /c\ tiegcm_namelist_file_name = 'tiegcm.nml'" \
    -e "/ file_out /c\ file_out = 'dart_ics'" \
    -e "/ single_restart_file_in /c\ single_restart_file_in = .FALSE." \
    -e "/ async /c\ async = 2" \
    -e "/ adv_ens_command /c\ adv_ens_command = './advance_model.csh'" \
    -e "/ start_from_restart /c\ start_from_restart = .TRUE." \
    -e "/ restart_in_file_name /c\ restart_in_file_name = 'filter_ics'" \
    -e "/ restart_out_file_name /c\ restart_out_file_name = 'filter_restart'" \
    -e "/ file_in /c\ file_in = 'dart_restart'" \
    -e "/ file_namelist_out /c\ file_namelist_out = 'namelist_update'" \
    input.nml.original >! input.nml  || exit 2

set ENSEMBLESTRING = `grep -A 42 filter_nml input.nml | grep ens_size`
set NUM_ENS = `echo $ENSEMBLESTRING[3] | sed -e "s#,##"`

@ instance = 1
while ( $instance <= $NUM_ENS )

  set darticname  = `printf "filter_ics.%04d"          $instance`
  set tiesecond   = `printf "tiegcm_s.nc.%04d"         $instance`
  set tierestart  = `printf "tiegcm_restart_p.nc.%04d" $instance`
  set tieinp      = `printf "tiegcm.nml.%04d"          $instance`

  ${COPY} ${ENSEMBLEDIR}/$tiesecond  .                   || exit 2
  ${COPY} ${ENSEMBLEDIR}/$tierestart .                   || exit 2
  ${COPY} ${ENSEMBLEDIR}/$tieinp     tiegcm.nml.original || exit 2

  # Ensure that the tiegcm.nml for all the ensemble members is identical
  # in all the ways that matter. This will result in a miniumum of changes
  # in the advance_model.csh script. This script REQUIRES that there is a  
  # SINGLE tiegcm_restart_p.nc. Just keep appending all the timesteps to
  # the same file. If you need to subset the large file, use the NCO
  # operators. for example    ncks -d time,20,30 tiegcm_restart_p.nc bob.nc 
  # If you need more than 300 timesteps in the file, increase it here.
  
  sed -e 's/;.*//' -e '/^$/ d' \
      -e "/ MXHIST_PRIM /c\ MXHIST_PRIM = 300" \
      -e "/ MXHIST_SECH /c\ MXHIST_SECH = 300" \
      -e "/ SOURCE /c\ SOURCE = 'tiegcm_restart_p.nc'" \
      -e "/ OUTPUT /c\ OUTPUT = 'tiegcm_restart_p.nc'" \
      -e "/ SECOUT /c\ SECOUT = 'tiegcm_s.nc'"         \
      tiegcm.nml.original >! $tieinp  || exit 2

  # If an existing ensemble of filter_ics.#### exist, use it.
  # If not, generate one. Be aware - even if they exist, they may
  # not have the same variable set as your current input.nml
  # If that is the case, you will have to generate your own set anyway.
  # If you get an error from aread_state_restart(), this is likely the case.

  if (  -e  ${ENSEMBLEDIR}/initial/$darticname.GENERATE ) then
     ${REMOVE} $darticname
     ${LINK} ${ENSEMBLEDIR}/initial/$darticname . || exit 2
  else
     # We must convert a tiegcm_restart_p.nc file to a dart_ics file
     # for each ensemble member. So - momentarily, we must
     # create links to the static filenames expected by model_to_dart

     ${REMOVE} tiegcm_restart_p.nc tiegcm_s.nc tiegcm.nml 

     ${LINK} $tierestart tiegcm_restart_p.nc   || exit 2
     ${LINK} $tiesecond  tiegcm_s.nc           || exit 2
     ${LINK} $tieinp     tiegcm.nml            || exit 2

     ./model_to_dart || exit 2

     if (-e dart_ics ) then
        ${MOVE} dart_ics $darticname
     else
        echo "ERROR: File conversion from $tierestart to $darticname failed."
        echo "ERROR: File conversion from $tierestart to $darticname failed."
        echo "ERROR: File conversion from $tierestart to $darticname failed."
        exit 2
     endif
  endif

  @ instance++
end

#-----------------------------------------------------------------------------
# Run filter ... 
#-----------------------------------------------------------------------------

${REMOVE} tiegcm_restart_p.nc tiegcm_s.nc tiegcm.nml 

${LINK} tiegcm_restart_p.nc.0001 tiegcm_restart_p.nc   || exit 3
${LINK} tiegcm_s.nc.0001         tiegcm_s.nc           || exit 3
${LINK} tiegcm.nml.0001          tiegcm.nml            || exit 3

${MPI_RUN_CMD} ./filter || exit 3

#-----------------------------------------------------------------------------
# At this point, all the restart,diagnostic files are in the run/CENTRALDIR.
# You may want to move them to someplace more 'permanent'.
#
# TJH: At this point, the output files have pretty 'generic' names.
# The files could be archived with the assimilation date in their name.
#-----------------------------------------------------------------------------

# ${COPY} tiegcm.nml                 ${EXPERIMENT}/tiegcm
# ${MOVE} tiegcm_s.nc*               ${EXPERIMENT}/tiegcm
# ${MOVE} tiegcm_restart_p.nc*       ${EXPERIMENT}/tiegcm
# ${MOVE} tiegcm_out_*               ${EXPERIMENT}/tiegcm

# ${MOVE} Posterior_Diag.nc          ${EXPERIMENT}/DART
# ${MOVE} Prior_Diag.nc              ${EXPERIMENT}/DART
# ${MOVE} obs_seq.final              ${EXPERIMENT}/DART
# ${MOVE} dart_log.out               ${EXPERIMENT}/DART

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

