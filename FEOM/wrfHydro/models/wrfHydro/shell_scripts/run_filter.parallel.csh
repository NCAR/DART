#!/bin/csh


#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
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
# jlm - there were previously BSUB directives and code, which I removed.
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
## -l nodes=1:ppn=1   (ppn == Processors Per Node)
##=============================================================================
#PBS -N run_filter
#PBS -r n
#PBS -e filter.err
#PBS -o filter.log
#PBS -l nodes=3:ppn=11,walltime=168:00:00

####################################################################
####################################################################

# in batch mode you have to cd to the work dir.  if run interactively
# don't mess with where we are.
if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

echo `pwd`

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

# if async=2, e.g. you are going to run './modelxxx', single process
# (or possibly 'mpirun -np 1 ./modelxxx'), so each processor advances
# one ensemble independently of the others, leave this as false.
#
# if async=4, e.g. all the processors advance each modelxxx in turn with
# mpirun -np 64 modelxxx (or whatever) for as many ensembles as you have,
# set this to "true"

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

set ADVANCESTRING = `grep -A 42 filter_nml input.nml | grep adv_ens_command | grep -v '!'`
set ADV_CMD  = `echo $ADVANCESTRING | cut -d= -f2 | tr -d '"'`
echo $ADV_CMD

# this variable is optional. helps manage writing files to node disks rather than 
# across the network. 
# nodeDisk should be a location available on all nodes for writing. 
set nodeDisk = '/c1'
if ($?nodeDisk) then 
    set nodeDir  = ${nodeDisk}/DART.RUN.`date | tr ' ' '-'`
    echo $nodeDir > nodeDir.control
    if ( ! -e nodeDir.control ) exit 1
endif

echo parallel_model: $parallel_model

# A common strategy for the beginning is to check for the existence of
# some variables that get set by the different queuing mechanisms.
# This way, we know which queuing mechanism we are working with,
# and can set 'queue-independent' variables for use for the remainder
# of the script.
if ($?PBS_O_WORKDIR) then

    #####################################################################
    #####################################################################
    # PBS has a list of processors in a file whose name is (PBS_NODEFILE)
    echo "PBS - using mpirun for execution"
    # each filter task advances the ensembles, each running on 1 proc.
    if ( "$parallel_model" == "false" ) then  
    ## async==2 ##

      mpirun ./filter

    else  
    ## async==4 ##
   
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

    #####################################################################
    #####################################################################
    # If you have a linux cluster with no queuing software, use this
    # section. The list of computational nodes is given to the mpirun
    # command and it assigns them as they appear in the file. In some
    # cases it seems to be necessary to wrap the command in a small
    # script that changes to the current directory before running.

    echo "running with no queueing system from current directory"
    
    # USER CONFIG
    # before running this script
    # node name : how many tasks you can run on it
    #setenv MYNODEFILE `pwd`/nodelist
    #echo "node5:5" >! $MYNODEFILE
    #echo "node5:2" >> $MYNODEFILE
    #echo "node3:2" >> $MYNODEFILE
    #echo "node1:2" >> $MYNODEFILE
    setenv NUM_PROCS 1
    set MPIRUN = '/opt/openmpi/bin/mpirun'
    set MPICMD = "$MPIRUN -np $NUM_PROCS"
    # -machinefile $MYNODEFILE"
    
    echo "MPICMD = ${MPICMD}"

    if ( "$parallel_model" == "false" ) then  
    ## async==2 ##

      $MPICMD ./filter

      if (! $?) exit 8

    else  
    ## async==4 ##

        # filter runs in parallel until time to do a model advance,
        # and then this script starts up the modelxxx jobs, each one
        # running in parallel. then it runs wakeup_filter to wake
        # up filter so it can continue.

        \rm -f model_to_filter.lock filter_to_model.lock
        mkfifo model_to_filter.lock filter_to_model.lock

        set filterhome = $HOME/.filter$$
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
            # jlm somewhat unclear. 
            # and also giving an example calling seq in input.nml
            # e.g.: mpirun -np $NUM_PROCS ?
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

    endif  ## async if 2 else 4

endif  ## ($?PBS_O_WORKDIR) 

## bring back desired output on the nodes.  
if ($?nodeDisk) then 
    ./getNodeFiles.sh `pwd`
    rm nodeDir.control
endif


exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

