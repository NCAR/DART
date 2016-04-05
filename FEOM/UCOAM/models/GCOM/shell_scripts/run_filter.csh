#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Script to assimilate observations using DART and the GCOM ocean model.
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
## the normal way to submit to the queue is:    qsub run_filter
##
## an explanation of the most common directives follows:
## -N <arg>   Job name
## -r n       Declare job non-rerunable
## -e <arg>   filename for standard error AFTER job completes
## -o <arg>   filename for standard out AFTER job completes
## -q <arg>   Queue name
## -l nodes=xx:ppn=16   request xx nodes and 16 processors on each node.
## -l walltime=hh:mm:ss request hh wallclock hours of runtime ..
##PBS -l nodes=node12+node13
##PBS -l nodes=reserved
##PBS -l nodes=node13:ppn=15+node9:ppn=15
##PBS -l nodes=2:ppn=15:reserved
##PBS -l walltime=24:00:00
##==============================================================================
#
#PBS -N gcom_filter
#PBS -r n
#PBS -e gcom_filter.err
#PBS -o gcom_filter.log
#PBS -q batch
#PBS -l nodes=2:ppn=15
#
##==============================================================================
## This block of directives constitutes the preamble for the LSF queuing system
##
## the normal way to submit to the queue is:    bsub < run_filter
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
#BSUB -J gcom_filter
#BSUB -o gcom_filter.%J.log
#BSUB -P NIMGxxxx
#BSUB -q premium
#BSUB -W 1:00
#BSUB -n 16
#BSUB -R "span[ptile=16]"
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
   setenv ORIGINALDIR   $PBS_O_WORKDIR
   setenv JOBNAME       $PBS_JOBNAME
   setenv JOBID         $PBS_JOBID
   setenv MYQUEUE       $PBS_QUEUE
   setenv MYHOST        $PBS_O_HOST
   setenv NCORES        `cat $PBS_NODEFILE | wc -l`
   setenv RUN_CMD       "mpirun -np $NCORES -machinefile $PBS_NODEFILE"

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
   setenv JOBNAME     gcom
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $HOST
   setenv RUN_CMD     'no_interactive_mpi'

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
      setenv   MOVE 'mv -v'
      setenv   COPY 'cp -v --preserve=timestamps'
      setenv   LINK 'ln -vs'
      setenv REMOVE 'rm -fr'

      setenv  EXPERIMENT /glade/p/work/${USER}/${JOBNAME}
      setenv  CENTRALDIR /glade/scratch/${USER}/${JOBNAME}/job_${JOBID}
      setenv     DARTDIR ${HOME}/work/DART/UCOAM/models/GCOM
      setenv    SERUCOAM ${HOME}/work/DART/UCOAM/models/GCOM/serucoam
      setenv  BASEOBSDIR ${HOME}/work/DART/UCOAM/models/GCOM/work
      setenv ENSEMBLEDIR ${HOME}/work/DART/UCOAM/models/GCOM/serucoam
   breaksw

   case node*:
      # SDSU cluster "Cincinatti" has hosts 'node??'
      setenv   MOVE 'mv -v'
      setenv   COPY 'cp -v --preserve=timestamps'
      setenv   LINK 'ln -vs'
      setenv REMOVE 'rm -fr'

      setenv  EXPERIMENT /cinci/${USER}/deep_storage
      setenv  CENTRALDIR /usr/scratch/${USER}/${JOBNAME}/job_${JOBID}
      setenv     DARTDIR /home/${USER}/DART-UCOAM/models/GCOM
      setenv    SERUCOAM /home/${USER}/serucoamATMG-UVsmall
      setenv  BASEOBSDIR /home/${USER}/DART-UCOAM/models/GCOM/work
      setenv ENSEMBLEDIR /home/${USER}/serucoamATMG-UVsmall
   breaksw

   default:
      # SDSU "dulcinea"
      setenv   MOVE 'mv -v'
      setenv   COPY 'cp -v --preserve=timestamps'
      setenv   LINK 'ln -vs'
      setenv REMOVE 'rm -fr'
      setenv  EXPERIMENT /cinci/${USER}/deep_storage
      setenv  CENTRALDIR /cinci/${USER}/${JOBNAME}/job_${JOBID}
      setenv     DARTDIR /home/${USER}/DART-UCOAM/models/GCOM
      setenv    SERUCOAM /home/${USER}/serucoamAT
      setenv  BASEOBSDIR /home/${USER}/DART-UCOAM/models/GCOM/work
      setenv ENSEMBLEDIR /home/${USER}/DART-UCOAM/models/GCOM/serucoam
   breaksw

endsw

#-------------------------------------------------------------------------------
# Make a unique, (empty, clean) temporary directory.
#-------------------------------------------------------------------------------

echo "`date` -- BEGIN FILTER"

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
# This script requires that gcom is serial code - but DART's filter
# and wakeup_filter must be compiled WITH mpi.

echo "`date` -- Assembling the GCOM pieces."

# cd ${SERUCOAM}/src
# make clean || exit -1
# make       || exit -1
# ${REMOVE} *.o *.mod
# cd ${CENTRALDIR}

${COPY} ${SERUCOAM}/Main.exe          gcom.serial.exe || exit 1
${COPY} ${SERUCOAM}/Grid.dat          Grid.dat        || exit 1
${COPY} ${SERUCOAM}/Gridll.dat        Gridll.dat      || exit 1
${COPY} ${SERUCOAM}/ProbSize.dat      ProbSize.dat    || exit 1
${COPY} ${SERUCOAM}/param.dat         param.dat       || exit 1

#===============================================================================
# Block 2: Populate CENTRALDIR with everything needed to run DART and GCOM.
#===============================================================================

# Get the DART executables, scripts, and input files
# The input.nml will be copied from the DART directory and modified appropriately.
# Again, if the compute nodes can execute code compiled on the head node,
# it should not be necessary to compile the code here.

echo "`date` -- Assembling the DART pieces"

# cd ${DARTDIR}/work
# csh quickbuild.csh -mpi || exit 2
# cd ${CENTRALDIR}

${COPY} ${DARTDIR}/work/advance_time               . || exit 2
${COPY} ${DARTDIR}/work/filter                     . || exit 2
${COPY} ${DARTDIR}/work/dart_to_gcom               . || exit 2
${COPY} ${DARTDIR}/work/gcom_to_dart               . || exit 2
${COPY} ${BASEOBSDIR}/obs_seq.out                  . || exit 2
${COPY} ${DARTDIR}/shell_scripts/advance_model.csh . || exit 2

#===============================================================================
# Block 3: convert N GCOM restart files to DART initial conditions file(s).
# At the end of this block, we have DART restart files   filter_ics.[1-N]
#===============================================================================
#
# DART namelist settings required (to make advance_model.csh easy):
#
# &filter_nml:           start_from_restart       = .true.
# &filter_nml:           async                    = 2
# &filter_nml:           restart_in_file_name     = 'dart_ics'
# &filter_nml:           obs_sequence_in_name     = 'obs_seq.in'
# &filter_nml:           obs_sequence_out_name    = 'obs_seq.perfect'
# &filter_nml:           init_time_days           = -1,
# &filter_nml:           init_time_seconds        = -1,
# &filter_nml:           first_obs_days           = -1,
# &filter_nml:           first_obs_seconds        = -1,
# &filter_nml:           last_obs_days            = -1,
# &filter_nml:           last_obs_seconds         = -1,
# &ensemble_manager_nml: single_restart_file_in   = .false.
# &model_nml:            gcom_restart_file        = 'gcom_restart.nc'
# &model_nml:            gcom_geometry_file       = 'gcom_geometry.nc'
# &gcom_to_dart_nml:     gcom_to_dart_output_file = 'dart_ics'
# &dart_to_gcom_nml:     dart_to_gcom_input_file  = 'dart_restart'
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
    -e "/ obs_sequence_in_name /c\ obs_sequence_in_name = 'obs_seq.out'" \
    -e "/ obs_sequence_out_name /c\ obs_sequence_out_name = 'obs_seq.final'" \
    -e "/ single_restart_file_in /c\ single_restart_file_in = .false." \
    -e "/ gcom_restart_file /c\ gcom_restart_file = 'gcom_restart.nc'" \
    -e "/ gcom_geometry_file /c\ gcom_geometry_file = 'gcom_geometry.nc'" \
    -e "/ gcom_to_dart_output_file /c\ gcom_to_dart_output_file = 'dart_ics'" \
    -e "/ dart_to_gcom_input_file /c\ dart_to_gcom_input_file = 'dart_restart'" \
       ${DARTDIR}/work/input.nml >! input.nml || exit 3

#-------------------------------------------------------------------------------
# Grab the start time from the param.dat file and convert it to a DART time.
# Use that as the start time for all the initial conditions files.
#-------------------------------------------------------------------------------

set starttime = `grep Start_Time param.dat | sed -e "s/[:=,a-zA-Z_'\-]/ /g"`
set   startyear = $starttime[1]
set  startmonth = $starttime[2]
set    startday = $starttime[3]
set   starthour = $starttime[4]
@ startminute = `echo $starttime[5] | bc`
@ startsecond = `echo $starttime[6] | bc`
@ seconds = ${startminute} * 60 + ${startsecond}

set modeltime = `echo "${startyear}${startmonth}${startday}${starthour} +${seconds}s -g" | ./advance_time`

set modeldays = $modeltime[1]
set modelseconds = $modeltime[2]

echo "$starttime is DART time $modeldays $modelseconds"

#-------------------------------------------------------------------------------
# Determine the number of ensemble members from input.nml,
# It may exist in more than one place - we will use the first instance.
# Parse out the filter_nml string and use the next hunk of lines.
# ditto for the advance command
#-------------------------------------------------------------------------------

set ENSEMBLESTRING = `grep -A 42 filter_nml input.nml | grep ens_size`
set  ADVANCESTRING = `grep -A 42 filter_nml input.nml | grep adv_ens_command`
set    ASYNCSTRING = `grep -A 42 filter_nml input.nml | grep async`
set  ensemble_size = `echo $ENSEMBLESTRING[3] | sed -e "s#,##"`
set        ADV_CMD = `echo  $ADVANCESTRING[3] | sed -e 's#,##' -e 's#"##g'`
set     ASYNC_TYPE = `echo    $ASYNCSTRING[3] | sed -e 's#,##'`

echo "The model advance command is ${ADV_CMD}"

#-------------------------------------------------------------------------------
# Loop over all the ensemble members and create a DART ics file for each.
# This requires an ensemble of gcom restart files "gcom_restart_nnnn.nc"
# GCOM_POINTER will always point to the most current GCOM file for each member.
#-------------------------------------------------------------------------------

echo "`date` -- BEGIN GCOM-TO-DART"

set member = 1
while ($member <= $ensemble_size)

   set  GCOM_POINTER = `printf gcom_restart_%04d.nc $member`
   set  rest_file_in = `printf   dart_input.%04d    $member`
   set  DART_IC_FILE = `printf     dart_ics.%04d    $member`

   ${REMOVE} gcom_restart.nc gcom_geometry.nc

   ${COPY} ${ENSEMBLEDIR}/${GCOM_POINTER}  .  || exit 3
   ${LINK} ${GCOM_POINTER} gcom_restart.nc    || exit 3
   ${LINK} ${GCOM_POINTER} gcom_geometry.nc   || exit 3

   ./gcom_to_dart                       || exit 3
   ${MOVE} dart_ics ${DART_IC_FILE}     || exit 3

   @ member++

end

#===============================================================================
# Block 4: run filter in the appropriate manner.
#===============================================================================

# detect whether the model is supposed to run as an MPI job or not
# by reading the "async = " from the &filter_nml namelist in input.nml.

if ( "${ASYNC_TYPE}" == "0" || "${ASYNC_TYPE}" == "2") then
  set parallel_model = "false"
  echo "The model is believed to be single-threaded."
else if ( "${ASYNC_TYPE}" == "4") then
  set parallel_model = "true"
  echo "The model is believed to be MPI-aware."
else
  echo 'ERROR - Cannot autodetect async value in the filter_nml namelist in input.nml.'
  echo 'ERROR - hardcode the parallel_model shell variable and comment out these lines.'
  exit -1
  set parallel_model = "false"
endif

#-------------------------------------------------------------------------------
# A common strategy for the beginning is to check for the existence of
# some variables that get set by the different queuing mechanisms.
# This way, we know which queuing mechanism we are working with,
# and can set 'queue-independent' variables for use for the remainder
# of the script.

if ($?LSB_QUEUE || $?PBS_QUEUE) then

    # Must be using LSF or PBS as the queueing system.
    echo "Using ${RUN_CMD} for execution"

    if ( "$parallel_model" == "false" ) then
       # each filter task advances the ensembles, each running on 1 proc.
       ${RUN_CMD} ./filter || exit 4

    else
       # 1) filter runs in parallel until time to do a (parallel) model advance.
       # 2) advance_model.csh successively runs N GCOM instances,
       #    each using the entire processor set - one after another.
       # 3) wakeup_filter wakes up filter so it can continue.

      \rm -f model_to_filter.lock filter_to_model.lock
      mkfifo model_to_filter.lock filter_to_model.lock

      set filterhome = ~/.filter$$
      if ( ! -e $filterhome) mkdir $filterhome

      # this starts filter but also returns control back to
      # this script immediately.

      ( setenv HOME $filterhome; ${RUN_CMD} ./filter ) &

      while ( -e filter_to_model.lock )

        set todo=`cat < filter_to_model.lock`
        echo "todo received, value = ${todo}"

        if ( "${todo}" == "finished" ) then
          echo "main script: filter done."
          wait
          break

        else if ( "${todo}" == "advance" ) then

          echo "calling model advance now:"
          ./advance_model.csh 0 ${ensemble_size} filter_control00000 || exit 9

          echo "restarting filter."
          ${RUN_CMD} ./wakeup_filter

        else

          echo "main script: unexpected value received."
          break

        endif

      end

      echo "filter finished, removing pipes."
      \rm -f model_to_filter.lock filter_to_model.lock

      if ( -d $filterhome) rmdir $filterhome
    endif

else

    # If you have a linux cluster with no queuing software, use this
    # section. The list of computational nodes is given to the mpirun
    # command and it assigns them as they appear in the file. In some
    # cases it seems to be necessary to wrap the command in a small
    # script that changes to the current directory before running.

    echo "running with no queueing system"
    echo "This is untested for GCOM -- ending now."
    exit 9

    # before running this script, do this once. the syntax is
    # node name : how many tasks you can run on it
    #setenv MYNODEFILE ~/nodelist
    #echo "node7:2" >! $MYNODEFILE
    #echo "node5:2" >> $MYNODEFILE
    #echo "node3:2" >> $MYNODEFILE
    #echo "node1:2" >> $MYNODEFILE

#   for compas
    setenv NUM_PROCS `cat nodelist-pgi | wc -l`
    set MPIRUN = /opt/mpich/myrinet/pgi/bin/mpirun
    setenv RUN_CMD "$MPIRUN -np $NUM_PROCS -nolocal -machinefile nodelist-pgi"

#   for atlas-pgi
    setenv NUM_PROCS `cat nodelist-pgi | wc -l`
    set MPIRUN = /share/apps/mpich1/pgi/bin/mpirun
    setenv RUN_CMD "$MPIRUN -np $NUM_PROCS -nolocal -machinefile nodelist-pgi"

#   for atlas-gfortran
    set MPIRUN = /share/apps/openmpi/gfortran/bin/mpirun
    setenv RUN_CMD "$MPIRUN --hostfile nodelist-gfortran --mca mtl mx --mca pml cm -np 72"

    echo "RUN_CMD = ${RUN_CMD}"

    # filter runs in parallel until time to do a model advance,
    # and then this script starts up the mitgcmuv jobs, each one
    # running in parallel. then it runs wakeup_filter to wake
    # up filter so it can continue.

    \rm -f model_to_filter.lock filter_to_model.lock
    mkfifo model_to_filter.lock filter_to_model.lock

    set filterhome = ~/.filter$$
    if ( ! -e $filterhome) mkdir $filterhome

    # this starts filter but also returns control back to
    # this script immediately.

    (setenv HOME $filterhome; ${RUN_CMD} ./filter) &

    while ( -e filter_to_model.lock )

        set todo=`cat < filter_to_model.lock`
        echo "todo received, value = ${todo}"

        if ( "${todo}" == "finished" ) then
          echo "main script: filter done."
          wait
          break

        else if ( "${todo}" == "advance" ) then

          # the second number below must match the number
          # of ensembles. Also, in input.nml, the advance model
          # command must have -np N with N equal to the number
          # of processors this job is using.

          echo "calling model advance now:"
          ./advance_model.csh 0 ${ensemble_size} filter_control00000 || exit 9

          echo "restarting filter."
          ${RUN_CMD} ./wakeup_filter

        else

          echo "main script: unexpected value received."
          break

        endif

    end

    echo "filter finished, removing pipes."
    \rm -f model_to_filter.lock filter_to_model.lock

    if ( -d $filterhome) rmdir $filterhome

endif

#===============================================================================
# Block 5: Think about it. The last GCOM restart file is not modified at this point.
# The assimilated state is in the DART filter_restart.nnnnn.
# If you want to create a GCOM posterior state from this file, you will
# have to manually run dart_to_gcom on this last set of DART and GCOM files.
#===============================================================================

# FIXME ... should update the last gcom_restart files with the assimilated state.

#===============================================================================
# Block 6: Move the output to storage after filter completes.
# At this point, all the restart,diagnostic files are in the CENTRALDIR
# and need to be moved to the 'experiment permanent' directory.
# We have had problems with some, but not all, files being moved
# correctly, so we are adding bulletproofing to check to ensure the filesystem
# has completed writing the files, etc. Sometimes we get here before
# all the files have finished being written.
#===============================================================================

echo "${JOBNAME} ($JOBID) finished at "`date`
echo "Listing contents of CENTRALDIR before archiving"
ls -l

exit 0

#-------------------------------------------------------------------------------
# everything after here is just bits and pieces that may be useful.
# the preceeding exit means that we never get here.
#
# ${MOVE} *.data *.meta              ${experiment}/GCOM
# ${MOVE} data data.cal              ${experiment}/GCOM
# ${MOVE} STD*                       ${experiment}/GCOM

# ${MOVE} filter_restart*            ${experiment}/DART
# ${MOVE} assim_model_state_ud[1-9]* ${experiment}/DART
# ${MOVE} assim_model_state_ic[1-9]* ${experiment}/DART
# ${MOVE} Posterior_Diag.nc          ${experiment}/DART
# ${MOVE} Prior_Diag.nc              ${experiment}/DART
# ${MOVE} obs_seq.final              ${experiment}/DART
# ${MOVE} dart_log.out               ${experiment}/DART

# Good style dictates that you save the scripts so you can see what worked.

# ${COPY} input.nml                  ${experiment}/DART
# ${COPY} *.csh                      ${experiment}/DART
# ${COPY} $myname                    ${experiment}/DART

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

