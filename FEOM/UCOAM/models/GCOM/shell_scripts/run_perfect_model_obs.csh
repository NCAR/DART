#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Top level script to generate observations and a TRUE state.
# This presumes two directories exists that contain all the required bits
# for GCOM and for DART.
#
# This script only processes a single observation file.
# Still fairly complex; requires a raft of
# data files and most of them are in hardcoded locations.
#
# This script is designed to be executed as a batch job.
# However, if you comment out the model advance part - you can run this
# interactively to check for logic, file motion, syntax errors and the like.
# It is entirely fine to have directives for both PBS and LSF in the same file.
# I guarantee that as soon as you delete one set, you will change machines
# and wish you had not deleted the directives.
#
# The script moves the necessary files to a temporary directory that is the
# basis for the DART experiment; this will be called CENTRALDIR.
#
##==============================================================================
## This block of directives constitutes the preamble for the PBS queuing system
##
## the normal way to submit to the queue is:    qsub run_perfect_model_obs.csh
##
## an explanation of the most common directives follows:
## -N <arg>   Job name
## -r n       Declare job non-rerunable
## -e <arg>   filename for standard error AFTER job completes
## -o <arg>   filename for standard out AFTER job completes
## -q <arg>   Queue name
## -l nodes=xx:ppn=16   request xx nodes and 16 processors on each node.
## -l walltime=hh:mm:ss request hh wallclock hours of runtime ..
##==============================================================================
#
#PBS -N gcom_pmo
#PBS -r n
#PBS -e gcom_pmo.err
#PBS -o gcom_pmo.out
#PBS -q batch
#PBS -l nodes=1:ppn=1:reserved
#PBS -l walltime=24:00:00
#
##==============================================================================
## This block of directives constitutes the preamble for the LSF queuing system
##
## the normal way to submit to the queue is:    bsub < run_perfect_model_obs.csh
##
## an explanation of the most common directives follows:
## -J <arg>      Job name
## -o <arg>      output listing filename
## -P <arg>      account
## -q <arg>      queue
## -W <arg>      wall-clock hours:minutes required
## -n <arg>      number of MPI tasks (or processors, usually)
## -R <arg>      number of MPI tasks per node
## -N -u <arg>   mail this user when job finishes
##==============================================================================
#
#BSUB -J gcom_pmo
#BSUB -o gcom_pmo.%J.log
#BSUB -P 8685xxxx
#BSUB -q regular
#BSUB -W 2:00
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -N -u ${USER}@ucar.edu

#-------------------------------------------------------------------------------
# Turns out the scripts are a lot more flexible if you don't rely on
# the queuing-system-specific variables -- so I am converting them to
# 'generic' names and using the generics throughout the remainder.
#-------------------------------------------------------------------------------

if ($?PBS_QUEUE) then

   #----------------------------------------------------------------------------
   # This is used by PBS
   #----------------------------------------------------------------------------

   setenv ORIGINALDIR $PBS_O_WORKDIR
   setenv JOBNAME     $PBS_JOBNAME
   setenv JOBID       $PBS_JOBID
   setenv MYQUEUE     $PBS_QUEUE
   setenv MYHOST      $PBS_O_HOST
  setenv RUN_CMD     "mpirun -np 1 -machinefile $PBS_NODEFILE"

else if ($?LSB_QUEUE) then

   #----------------------------------------------------------------------------
   # This is used by LSF
   #----------------------------------------------------------------------------

   setenv ORIGINALDIR $LS_SUBCWD
   setenv JOBNAME     $LSB_JOBNAME
   setenv JOBID       $LSB_JOBID
   setenv MYQUEUE     $LSB_QUEUE
   setenv MYHOST      $LSB_SUB_HOST
   setenv RUN_CMD     mpirun.lsf

else

   #----------------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   #----------------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     pmo
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $HOST
   setenv RUN_CMD     csh

endif

#-------------------------------------------------------------------------------
# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.
#-------------------------------------------------------------------------------

set nonomatch  # suppress "rm" warnings if wildcard does not match anything

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following

switch ("`hostname`")

   case ys*:
      # NCAR "yellowstone"
      setenv   MOVE 'mv -fv'
      setenv   COPY 'cp -fv --preserve=timestamps'
      setenv   LINK 'ln -fvs'
      setenv REMOVE 'rm -fr'

      setenv EXPERIMENT /glade/p/work/${USER}/${JOBNAME}
      setenv CENTRALDIR /glade/scratch/${USER}/${JOBNAME}/job_${JOBID}
      setenv    DARTDIR ${HOME}/work/DART/UCOAM/models/GCOM
      setenv   SERUCOAM ${HOME}/work/DART/UCOAM/models/GCOM/serucoam
      setenv BASEOBSDIR ${HOME}/work/DART/UCOAM/models/GCOM/work
   breaksw

   case node*:
      # SDSU cluster "Cincinatti" has hosts 'node??'
      setenv   MOVE 'mv -fv'
      setenv   COPY 'cp -fv --preserve=timestamps'
      setenv   LINK 'ln -fvs'
      setenv REMOVE 'rm -fr'

      setenv  EXPERIMENT /usr/scratch/${USER}/deep_storage
      setenv  CENTRALDIR /usr/scratch/${USER}/${JOBNAME}/job_${JOBID}
      setenv     DARTDIR /home/${USER}/DART-UCOAM/models/GCOM
      setenv    SERUCOAM /home/${USER}/serucoamATMG-UVsmall
#     setenv    SERUCOAM /home/${USER}/serucoamsvn
      setenv  BASEOBSDIR /home/${USER}/DART-UCOAM/models/GCOM/work
   breaksw

   default:
      # SDSU "dulcinea"
      setenv   MOVE 'mv -fv'
      setenv   COPY 'cp -fv --preserve=timestamps'
      setenv   LINK 'ln -fvs'
      setenv REMOVE 'rm -fr'

     # setenv EXPERIMENT /gcemproject/${USER}/${JOBNAME}
      setenv CENTRALDIR /cinci/${USER}/${JOBNAME}/job_${JOBID}
      setenv    DARTDIR /home/${USER}/DART-UCOAM/models/GCOM
      setenv   SERUCOAM /home/${USER}/serucoamAT
      setenv BASEOBSDIR /home/${USER}/DART-UCOAM/models/GCOM/work
   breaksw

endsw

#-------------------------------------------------------------------------------
# Make a unique, (empty, clean) temporary directory.
#-------------------------------------------------------------------------------

echo "`date` -- BEGIN PERFECT_OBS_GENERATION"

mkdir -p ${CENTRALDIR}
cd ${CENTRALDIR}

set myname = $0          # this is the name of this script

#-------------------------------------------------------------------------------
# Just an echo of the job attributes
#-------------------------------------------------------------------------------

echo
echo "${JOBNAME} ($JOBID) submitted   from $ORIGINALDIR"
echo "${JOBNAME} ($JOBID) submitted   from $MYHOST"
echo "${JOBNAME} ($JOBID) running in queue $MYQUEUE"
echo "${JOBNAME} ($JOBID) running       on $HOST"
echo "${JOBNAME} ($JOBID) started       at "`date`
echo "${JOBNAME} ($JOBID) CENTRALDIR    is $CENTRALDIR"
echo

#===============================================================================
# Block 1: Build all the GCOM executables we will need for this run.
# Since the compute nodes cannot execute things compiled on the head node,
# you have to compile what you need on the compute node. Really annoying.
# This script requires that gcom is serial code - no MPI at this point.

echo "`date` -- Assembling the GCOM pieces."

# cd ${SERUCOAM}/src
# make clean || exit -1
# make       || exit -1
# ${REMOVE} *.o *.mod
# cd ${CENTRALDIR}

${COPY} ${SERUCOAM}/Main.exe          gcom.serial.exe || exit 1
${COPY} ${SERUCOAM}/Grid.dat          Grid.dat        || exit 1
${COPY} ${SERUCOAM}/Gridll.dat          Gridll.dat        || exit 1
${COPY} ${SERUCOAM}/ProbSize.dat      ProbSize.dat    || exit 1
${COPY} ${SERUCOAM}/param.dat         param.dat       || exit 1
${COPY} ${SERUCOAM}/gcom_restart.nc   gcom_restart.nc || exit 1

#===============================================================================
# Block 2: Populate CENTRALDIR with everything needed to run DART and GCOM.
#===============================================================================

# Get the DART executables, scripts, and input files
# The input.nml will be copied from the DART directory and modified appropriately.
# Again, if the compute nodes can execute code compiled on the head node,
# it should not be necessary to compile the code here.

echo "`date` -- Assembling the DART pieces"

# cd ${DARTDIR}/work
# csh quickbuild.csh -nompi || exit 2
# cd ${CENTRALDIR}

${COPY} ${DARTDIR}/work/perfect_model_obs          . || exit 2
${COPY} ${DARTDIR}/work/dart_to_gcom               . || exit 2
${COPY} ${DARTDIR}/work/gcom_to_dart               . || exit 2
${COPY} ${BASEOBSDIR}/obs_seq.in                   . || exit 2
${COPY} ${DARTDIR}/shell_scripts/advance_model.csh . || exit 2

#===============================================================================
# Block 3: Convert 1 GCOM restart file to a DART initial conditions file.
# At the end of the block, we have a DART initial condition file  dart_ics
#===============================================================================
#
# DART namelist settings required (to make advance_model.csh easy):
#
# &perfect_model_obs_nml:  start_from_restart       = .true.
# &perfect_model_obs_nml:  async                    = 2
# &perfect_model_obs_nml:  restart_in_file_name     = 'dart_ics'
# &perfect_model_obs_nml:  obs_sequence_in_name     = 'obs_seq.in'
# &perfect_model_obs_nml:  obs_sequence_out_name    = 'obs_seq.perfect'
# &perfect_model_obs_nml:  init_time_days           = -1,
# &perfect_model_obs_nml:  init_time_seconds        = -1,
# &perfect_model_obs_nml:  first_obs_days           = -1,
# &perfect_model_obs_nml:  first_obs_seconds        = -1,
# &perfect_model_obs_nml:  last_obs_days            = -1,
# &perfect_model_obs_nml:  last_obs_seconds         = -1,
# &model_nml:              gcom_restart_file        = 'gcom_restart.nc'
# &model_nml:              gcom_geometry_file       = 'gcom_geometry.nc'
# &gcom_to_dart_nml:       gcom_to_dart_output_file = 'dart_ics'
# &dart_to_gcom_nml:       dart_to_gcom_input_file  = 'dart_restart'
#===============================================================================

if ( ! -e ${DARTDIR}/work/input.nml ) then
   echo "ERROR ... DART required file ${DARTDIR}/work/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${DARTDIR}/work/input.nml not found ... ERROR"
   exit 3
endif

# Ensure the namelist has the values required by this script.

sed -e "/ start_from_restart /c\ start_from_restart = .true." \
    -e "/ async /c\ async = 2" \
    -e "/ restart_in_file_name /c\ restart_in_file_name = 'dart_ics'" \
    -e "/ obs_sequence_in_name /c\ obs_sequence_in_name = 'obs_seq.in'" \
    -e "/ obs_sequence_out_name /c\ obs_sequence_out_name = 'obs_seq.perfect'" \
    -e "/ gcom_restart_file /c\ gcom_restart_file = 'gcom_restart.nc'" \
    -e "/ gcom_geometry_file /c\ gcom_geometry_file = 'gcom_geometry.nc'" \
    -e "/ gcom_to_dart_output_file /c\ gcom_to_dart_output_file = 'dart_ics'" \
    -e "/ dart_to_gcom_input_file /c\ dart_to_gcom_input_file = 'dart_restart'" \
       ${DARTDIR}/work/input.nml >! input.nml || exit 3

echo "`date` -- BEGIN GCOM-TO-DART"

${LINK} gcom_restart.nc gcom_geometry.nc || exit 3

${RUN_CMD} ./gcom_to_dart

if ($status != 0) then
   echo "ERROR ... DART died in 'gcom_to_dart' ... ERROR"
   echo "ERROR ... DART died in 'gcom_to_dart' ... ERROR"
   exit 3
endif

echo "`date` -- END GCOM-TO-DART"

#===============================================================================
# Block 4: Advance the model and harvest the synthetic observations.
# output files are:
# True_state.nc   ...... the DART state
# obs_seq.perfect ...... the synthetic observations
# dart_log.out    ...... run-time output of all DART routines
# perfect_restart ...... which we don't need
#===============================================================================

# advance_model.csh needs a 'gcom_restart_nnnn.nc' in this directory

${LINK} gcom_restart.nc gcom_restart_0001.nc || exit 4

echo "`date` -- BEGIN GCOM PERFECT_MODEL_OBS"

${RUN_CMD} ./perfect_model_obs

if ($status != 0) then
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   exit 4
endif

echo "`date` -- END   GCOM PERFECT_MODEL_OBS"

#===============================================================================
# Block 5: Copy/Move the good stuff back to someplace safe.
# CENTRALDIR is usually on a volatile or temporary filesystem.
# EXPERIMENT is usually someplace long-term.

mkdir -p ${EXPERIMENT}
${MOVE} True_State.nc    ${EXPERIMENT}
${MOVE} obs_seq.perfect  ${EXPERIMENT}
${MOVE} dart_log.out     ${EXPERIMENT}
${MOVE} input.nml        ${EXPERIMENT}
${MOVE} *.csh            ${EXPERIMENT}
${MOVE} $myname          ${EXPERIMENT}

echo "${JOBNAME} ($JOBID) finished at "`date`
echo "Listing contents of CENTRALDIR after archiving (miss anything?)"
ls -l

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

