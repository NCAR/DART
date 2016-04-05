#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# NOGAPS version.
#
# This script assumes that the experiment directory already has
# a full setup of files that nogaps needs to advance.  anything
# that is specific to a particular ensemble will be copied by
# the advance_model.csh script.
#
#=============================================================================
#
# This block of directives constitutes the preamble for the LSF queuing system
#
# the normal way to submit to the queue is:    bsub < run_filter.csh
#
# an explanation of the most common directives follows:
# -J	Job name 
# -o	STDOUT filename
# -e	STDERR filename
# -P	account
# -q	queue
# -n	number of processors  (really)
#
#=============================================================================
#
#BSUB -J filter
#BSUB -o filter.%J.log
#BSUB -q economy 
#BSUB -n 8 
#BSUB -W 12:00
#BSUB -R "span[ptile=1]"
#
#=============================================================================
#
# This block of directives constitutes the preamble for the PBS queuing system
#
# the normal way to submit to the queue is:    qsub run_filter.csh
#
# an explanation of the most common directives follows:
# -N     Job name
# -r n   Declare job non-rerunable
# -e <arg>  filename for standard error
# -o <arg>  filename for standard out
# -q <arg>   Queue name (small, medium, long, verylong)
# -l nodes=xx:ppn=2   requests two processors on the node.
#
#=============================================================================
#
#PBS -N filter
#PBS -r n
#PBS -e filter.err
#PBS -o filter.log
#PBS -q dedicated
#PBS -l nodes=16

if ($?LSB_QUEUE) then

   #-------------------------------------------------------------------
   # This is used by LSF
   #-------------------------------------------------------------------

   setenv ORIGINALDIR $LS_SUBCWD
   setenv JOBNAME     $LSB_JOBNAME
   setenv JOBID       $LSB_JOBID
   setenv MYQUEUE     $LSB_QUEUE
   setenv MYHOST      $LSB_SUB_HOST

else if ($?PBS_QUEUE) then

   #-------------------------------------------------------------------
   # This is used by PBS
   #-------------------------------------------------------------------

   setenv ORIGINALDIR $PBS_O_WORKDIR
   setenv JOBNAME     $PBS_JOBNAME
   setenv JOBID       $PBS_JOBID
   setenv MYQUEUE     $PBS_QUEUE
   setenv MYHOST      $PBS_O_HOST

else

   #-------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     filter
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $host

endif

# The experiment configuration is contained in these two files:
# config.csh
# input.nml
# 
# this is shared with the advance_model.csh script, and should set
# where the nogaps files, executables, etc are located.

source ${ORIGINALDIR}/config.csh

# If the directory doesn't exist yet, we don't know where to get
# the initial ensemble, so we should FAIL here. 
if ( ! -d ${experiment_dir}) then
    echo "ERROR There is no viable experiment directory."
    echo "ERROR ${ORIGINALDIR}/config.csh indicates it should be"
    echo ${experiment_dir}
    exit 1
else
    echo "Running the experiment in ${experiment_dir}"
endif

# Move into the experiment directory - DART calls this 'CENTRALDIR'

cd ${experiment_dir}

echo "copying all required files"

cp ${ORIGINALDIR}/input.nml                   . || exit 1
cp ${ORIGINALDIR}/config.csh                  . || exit 1
cp $FILTER_exec_name                          . || exit 1
cp $WAKEUPFILTER_exec_name                    . || exit 1
cp $NOGAPS2DART_exec_name                     . || exit 1
cp $DART2NOGAPS_exec_name                     . || exit 1
cp $TRANSTIME_exec_name                       . || exit 1
cp ${DART_script_dir}/advance_model.csh       . || exit 1
cp ${NOGAPS_exec_dir}/${NOGAPS_exec_name}     . || exit 1

# FIXME:
#  which observations are we going to use?  for now,
#  copy a single obs file into the same dir.
# set seqdate = `echo $dtg | cut -c1-8`
# cp $DART_obs/obs_seq${seqdate} obs_seq.out

# FIXME: for now, just use the obs from a perfect model
ln -fs /ptmp/hansenj/Exp1/test3/obs_seq.out .

# put in some bulletproofing here to be sure the files that
# the NOGAPS model_mod.f90 needs are indeed in the right place.
if ( ! -f noggeom.txt ) then
 cp -v ${climo}/noggeom${resolution}l${n_levels}_full.txt noggeom.txt
endif

# Determine the ensemble size from the filter_nml portion of input.nml
# Locate "filter_nml" and see which "ens_size" is immediately after it ...
# ditto for the advance command

set ENSEMBLESTRING = `grep -A 42 filter_nml input.nml | grep ens_size`
set        NUM_ENS = `echo $ENSEMBLESTRING[3] | sed -e "s#,##"`
set  ADVANCESTRING = `grep -A 42 filter_nml input.nml | grep adv_ens_command`
set        ADV_CMD = `echo $ADVANCESTRING[3] | sed -e 's#,##' -e 's#"##g'`

echo "The model advance command is ${ADV_CMD}"

# detect whether the model is supposed to run as an MPI job or not
# by reading the "async = " from the &filter_nml namelist in input.nml.
# if async=2, e.g. you are going to run './modelxxx', single process
# (or possibly 'mpirun -np 1 ./modelxxx'), so each processor advances
# one ensemble independently of the others, leave this as false.
#
# if async=4, e.g. all the processors advance each modelxxx in turn with
# mpirun -np 64 modelxxx (or whatever) for as many ensembles as you have,
# set this to "true"

# if async=4, also check that the call to advance_model.csh
# has the right number of ensemble members below; it must match
# the input.nml number.

set ASYNCSTRING = `grep -A 42 filter_nml input.nml | grep async`
set  ASYNC_TYPE = `echo $ASYNCSTRING[3] | sed -e 's#,##'`

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

# From here down, this script is pretty much like all the other
# models which work with dart.  above here, the changes are specific
# to NOGAPS.

# A common strategy for the beginning is to check for the existence of
# some variables that get set by the different queuing mechanisms.
# This way, we know which queuing mechanism we are working with,
# and can set 'queue-independent' variables for use for the remainder
# of the script.

if ($?LSB_QUEUE || $?PBS_QUEUE) then

    # Must be using LSF or PBS as the queueing system.
    echo "Using ${MPI} for execution"

    # each filter task advances the ensembles, each running on 1 proc.
    if ( "$parallel_model" == "false" ) then

       ${MPI} ./${FILTER_exec_name:t}

    else

    # filter runs in parallel until time to do a model advance,
    # and then this script starts up the modelxxx jobs, each one
    # running in parallel. then it runs wakeup_filter to wake
    # up filter so it can continue.

      \rm -f model_to_filter.lock filter_to_model.lock
      mkfifo model_to_filter.lock filter_to_model.lock

      set filterhome = ~/.filter$$
      if ( ! -e $filterhome) mkdir $filterhome

      # this starts filter but also returns control back to
      # this script immediately.

      ( setenv HOME $filterhome; ${MPI} ./${FILTER_exec_name:t} ) &

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
          ${ADV_CMD} 0 ${NUM_ENS} filter_control00000 || exit 9

          echo "restarting filter."
          ${MPI} ./${WAKEUPFILTER_exec_name:t} || exit 9

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

    # before running this script, do this once. the syntax is
    # node name : how many tasks you can run on it
    #setenv MYNODEFILE ~/nodelist
    #echo "node7:2" >! $MYNODEFILE
    #echo "node5:2" >> $MYNODEFILE
    #echo "node3:2" >> $MYNODEFILE
    #echo "node1:2" >> $MYNODEFILE

#   one possibility
    setenv NUM_PROCS `cat nodelist-pgi | wc -l`
    set MPIRUN = /opt/mpich/myrinet/pgi/bin/mpirun
    set MPICMD = $MPIRUN -np $NUM_PROCS -nolocal -machinefile nodelist-pgi

#   another possibility - note hardwired NP ...
    set MPIRUN = /share/apps/openmpi/gfortran/bin/mpirun
    set MPICMD = $MPIRUN --hostfile nodelist-gfortran --mca mtl mx --mca pml cm -np 72

    echo "MPICMD = ${MPICMD}"

    # filter runs in parallel until time to do a model advance,
    # and then this script starts up the modelxxx jobs, each one
    # running in parallel. then it runs wakeup_filter to wake
    # up filter so it can continue.

    \rm -f model_to_filter.lock filter_to_model.lock
    mkfifo model_to_filter.lock filter_to_model.lock

    set filterhome = ~/.filter$$
    if ( ! -e $filterhome) mkdir $filterhome

    # this starts filter but also returns control back to
    # this script immediately.

    (setenv HOME $filterhome; ${MPICMD} $FILTER_exec_name:t) &

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
          ${ADV_CMD} 0 ${NUM_ENS} filter_control00000 || exit 9

          echo "restarting filter."
          ${MPICMD} $WAKEUPFILTER_exec_name:t

        else

          echo "main script: unexpected value received."
          break

        endif

    end

    echo "filter finished, removing pipes."
    \rm -f model_to_filter.lock filter_to_model.lock

    if ( -d $filterhome) rmdir $filterhome

endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

# original collaboration information - do not edit
# $orgURL: https://svn2.assembla.com/svn/ngdart/shell_scripts/run_filter.csh $
# $orgId: run_filter.csh 107 2010-06-09 20:42:36Z thoar $
# $orgRevision: 107 $
# $orgDate: 2010-06-09 14:42:36 -0600 (Wed, 09 Jun 2010) $
