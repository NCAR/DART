#!/bin/tcsh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# Script to start an MPI version of filter, and then optionally
# run the model advance if &filter_nml has async=4 (parallel filter
# AND parallel model).  This version gets the number of ensemble members
# and advance command out of the input.nml namelist file automatically.
# It also gets the async setting and sets serial vs parallel model
# automatically.  The theory is that once you get this script working on
# your system, you will not need to change anything here as you change the
# number of ensemble members, async setting, or model advance command.
#
#=============================================================================
# This block of directives constitutes the preamble for the LSF queuing system
# LSF is used on the IMAGe Linux cluster 'coral'
# LSF is used on the IBM   'bluefire'
#
# the normal way to submit to the queue is:    bsub < run_filter.csh
#
# an explanation of the most common directives follows:
# -J Job_name
# -o STDOUT_filename
# -e STDERR_filename
# -P account_code_number
# -q queue    cheapest == [standby, economy, (regular,debug), premium] == $$$$
# -n number of MPI processes (not nodes)
# -W hh:mm  wallclock time (required on some systems)
##=============================================================================
#BSUB -J filter
#BSUB -o filter.%J.log
#BSUB -q standby
#BSUB -n 20
#BSUB -W 1:00
#
##=============================================================================
## This block of directives constitutes the preamble for the PBS queuing system
## PBS is used on the CGD Linux cluster 'bangkok'
## PBS is used on the CGD Linux cluster 'calgary'
##
## the normal way to submit to the queue is:    qsub run_filter.csh
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
##                     them both. (ppn == Processors Per Node)
##=============================================================================
#PBS -N filter
#PBS -r n
#PBS -e filter.err
#PBS -o filter.log
#PBS -q dedicated
#PBS -l nodes=10:ppn=2

# if async=2, e.g. you are going to run './modelxxx', single process
# (or possibly 'mpirun -np 1 ./modelxxx'), so each processor advances
# one ensemble independently of the others, leave this as false.
#
# if async=4, e.g. all the processors advance each modelxxx in turn with
# mpirun -np 64 modelxxx (or whatever) for as many ensembles as you have,
# set this to "true"

# this script is going to determine several things by reading the input.nml
# file which contains the &filter_nml namelist.  make sure it exists first.
if ( ! -e input.nml ) then
   echo "ERROR - input.nml does not exist in local directory."
   echo "ERROR - input.nml needed to determine several settings for this script."
   exit 1
endif

# detect whether the model is supposed to run as an MPI job or not
# by reading the "async = " from the &filter_nml namelist in input.nml.
# some namelists contain the same string - be sure to get the filter_nml one
# by grepping for lines which follow it.

set ASYNCSTRING = `grep -A 42 filter_nml input.nml | grep async`
set ASYNC_TYPE = `echo $ASYNCSTRING[3] | sed -e "s#,##"`

if ( "${ASYNC_TYPE}" == "0" || "${ASYNC_TYPE}" == "2") then
  set parallel_model = "false"
else if ( "${ASYNC_TYPE}" == "4") then
set parallel_model = "true"
else 
  echo 'cannot autodetect async value in the filter_nml namelist in input.nml file.'
  echo 'hardcode the parallel_model shell variable and comment out these lines.'
  exit -1
  set parallel_model = "false"
endif

# Determine the number of ensemble members from input.nml,
# as well as the command for advancing the model.

set ENSEMBLESTRING = `grep -A 42 filter_nml input.nml | grep ens_size`
set NUM_ENS = `echo $ENSEMBLESTRING[3] | sed -e "s#,##"`
set ADVANCESTRING = `grep -A 42 filter_nml input.nml | grep adv_ens_command`
set ADV_CMD  = `echo $ADVANCESTRING[3] | sed -e "s#,##"`

# A common strategy for the beginning is to check for the existence of
# some variables that get set by the different queuing mechanisms.
# This way, we know which queuing mechanism we are working with,
# and can set 'queue-independent' variables for use for the remainder
# of the script.

if ($?LS_SUBCWD) then

    # LSF has a list of processors already in a variable (LSB_HOSTS)
    echo "LSF - using mpirun.lsf for execution"

    # each filter task advances the ensembles, each running on 1 proc.
    if ( "$parallel_model" == "false" ) then

       mpirun.lsf ./filter

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

      (setenv HOME $filterhome; mpirun.lsf ./filter) &

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
          mpirun.lsf ./wakeup_filter

        else

          echo "main script: unexpected value received."
          break

        endif

      end

      echo "filter finished, removing pipes."
      \rm -f model_to_filter.lock filter_to_model.lock

      if ( -d $filterhome) rmdir $filterhome
    endif


else if ($?PBS_O_WORKDIR) then

    # PBS has a list of processors in a file whose name is (PBS_NODEFILE)
    echo "PBS - using mpirun for execution"

    # each filter task advances the ensembles, each running on 1 proc.
    if ( "$parallel_model" == "false" ) then

      mpirun ./filter

    else

    # filter runs in parallel until time to do a model advance,
    # and then this script starts up the modelxxx jobs, each one
    # running in parallel. then it runs wakeup_filter to wake
    # up filter so it can continue.

      \rm -f model_to_filter.lock filter_to_model.lock
      mkfifo model_to_filter.lock filter_to_model.lock

      set filterhome = ~/.filter
      if ( ! -e $filterhome) mkdir $filterhome

      # this starts filter but also returns control back to
      # this script immediately.

      (setenv HOME $filterhome; mpirun ./filter) &

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
          mpirun ./wakeup_filter

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

    (setenv HOME $filterhome; ${MPICMD} ./filter) &

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
          ${MPICMD} ./wakeup_filter

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

