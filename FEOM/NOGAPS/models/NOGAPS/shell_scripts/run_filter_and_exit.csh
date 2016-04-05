#!/bin/tcsh -v
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# NOGAPS version - this one tries to cut down the obs_seq file so it 
# only has obs for this assimilation window.  this version is intended
# to run filter, assimilate without trying to advance the model, and then
# exit.  a separate script will advance the model.
#
#=============================================================================
# This block of directives constitutes the preamble for the LSF queuing system
# LSF is used on the IMAGe Linux cluster 'coral'
# LSF is used on the IBM   'bluefire'
#
# the normal way to submit to the queue is:    bsub < run_filter.csh
#
# an explanation of the most common directives follows:
# -J Job name (master script job.csh presumes filter_server.xxxx.log)
# -o STDOUT filename
# -e STDERR filename
# -P      account
# -q queue    cheapest == [standby, economy, (regular,debug), premium] == $$$$
# -n number of processors  (really)
##=============================================================================
#BSUB -J filter
#BSUB -o filter.%J.log
#BSUB -q economy 
#BSUB -n 8 
#BSUB -W 1:00
#
##=============================================================================
## This block of directives constitutes the preamble for the PBS queuing system
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
## ? is ppn a general flag?
##=============================================================================
#PBS -N filter
#PBS -r n
#PBS -e filter.err
#PBS -o filter.log
#PBS -q dedicated
#PBS -l nodes=16

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

# if you are running and exiting without trying to advance the
# model, set this to false.
set parallel_model = "false"

# this is shared with the advance_model.csh script, and should set
# where the nogaps files, executables, etc are located.  see the
# script for more details of what it needs to set.
source ./config.csh

# if the directory doesn't exist yet, make it
if ( ! -d $scratch_dir/$experiment_name) mkdir -p $scratch_dir/$experiment_name

# copy any files needed over to $scratch_dir/$experiment_name 
# before cd'ing there and starting to run.
cp ./input.nml             $scratch_dir/$experiment_name
cp ./config.csh            $scratch_dir/$experiment_name
cp $FILTER_exec_name       $scratch_dir/$experiment_name
cp $WAKEUPFILTER_exec_name $scratch_dir/$experiment_name
cp $NOGAPS2DART_exec_name  $scratch_dir/$experiment_name
cp $DART2NOGAPS_exec_name  $scratch_dir/$experiment_name
cp $TRANSTIME_exec_name    $scratch_dir/$experiment_name
cp $OBSTOOL_exec_name      $scratch_dir/$experiment_name
cp $RESTART_exec_name      $scratch_dir/$experiment_name
# FIXME: this will eventually copy directly from the shell_scripts dir
cp ${DART_script_dir}/advance_model.csh       $scratch_dir/$experiment_name
cp ${DART_script_dir}/advance_model_batch.csh $scratch_dir/$experiment_name

# FIXME:
#  which observations are we going to use?  for now,
#  copy a single obs file into the same dir.
# set seqdate = `echo $dtg | cut -c1-8`
# cp $DART_obs/obs_seq${seqdate} $scratch_dir/$experiment_name/obs_seq.out

# FIXME: for now, just use the obs from perfect model
# because there is 1 observation per 6 hours which makes
# this fast to test.  eventually we want to use the real
# observations (~150,000 obs per 6 hours)

# do this in a subdirectory so we don't disturb the input.nml
# and the obs_sequence files here.
mkdir -p fixobs
cd fixobs

# this is the full observation sequence file to chop up
# it will be the input to the obs_sequence_tool; and obs_seq.out
# will be the output with only a 6 hour chunk of obs in it
cp -fv ../obs_seq.out obs_seq.full

cp -fv ${DART_script_dir}/input.nml.obstool input.nml.obstool
cp -fv $OBSTOOL_exec_name    .
cp -fv $ADVTIME_exec_name    .
set otool = ./${OBSTOOL_exec_name:t}   # just the filename, no path
set ttool = ./${ADVTIME_exec_name:t}   # just the filename, no path

set START=(`echo $dtg -3h+1s -g | $ttool`)
set STARTDAY=$START[1]
set STARTSEC=$START[2]
set END=(`echo $dtg +3h -g | $ttool`)
set ENDDAY=$END[1]
set ENDSEC=$END[2]
sed -e "s/STARTDAY/$STARTDAY/" \
    -e "s/STARTSEC/$STARTSEC/" \
    -e "s/ENDDAY/$ENDDAY/"     \
    -e "s/ENDSEC/$ENDSEC/"  input.nml.obstool > input.nml

# run it - this creates obs_seq.out in this directory
# (the name is set in the input.nml, in the namelist)
$otool

# and copy it into the execution directory
cp obs_seq.out       $scratch_dir/$experiment_name/obs_seq.out

cd ..

# now change into the scratch dir - this is sometimes
# called the 'centraldir' in dart parlance.
cd $scratch_dir/$experiment_name
echo current directory now `pwd`

# put in some bulletproofing here to be sure the files that
# the NOGAPS model_mod.f90 needs are indeed in the right place.
if ( ! -f noggeom.txt ) then
 cp -v $climo_dir/climo${resolution}/noggeom${resolution}l${n_levels}_full.txt noggeom.txt
endif

# Determine the number of ensemble members from input.nml,
# it may exist in more than one place.
# Parse out the filter_nml string and see which 
# one is immediately after it ...

if ( ! -e input.nml ) then
   echo "ERROR - input.nml does not exist in local directory."
   echo "ERROR - input.nml needed to determine number of ensemble members."
   exit 1
endif

set ENSEMBLESTRING = `grep -A 42 filter_nml input.nml | grep ens_size`
set NUM_ENS = `echo $ENSEMBLESTRING[3] | sed -e "s#,##"`

set ADVANCESTRING = `grep -A 42 filter_nml input.nml | grep advance`
set ADVANCECMD = `echo $ADVANCESTRING[3] | sed -e "s#,##"`

# the run_nogapsIC_to_dart.csh script leaves each of the dart state
# vectors in a subdirectory mbrNNN/dart_new_vector.  filter expects
# to read in filter_ic.NNNN in the current directory, so loop over
# the members making symbolic links.

@ i = 1
while ( $i <= $NUM_ENS )
  echo i is $i

  set index   = `printf "%03d" $i`
  set dirname = "mbr$index"
  set index   = `printf "%04d" $i`
  set icname  = "filter_ic.$index"

  ln -fsv $dirname/dart_new_vector $icname

  @ i += 1
end

# this script assumes that the current directory already has
# a full setup of files that nogaps needs to advance.  anything
# that is specific to a particular ensemble will be copied by
# the advance_model.csh script.

# from here down, this script is pretty much like all the other
# models which work with dart.  above here, the changes are specific
# to NOGAPS.


# A common strategy for the beginning is to check for the existence of
# some variables that get set by the different queuing mechanisms.
# This way, we know which queuing mechanism we are working with,
# and can set 'queue-independent' variables for use for the remainder
# of the script.

if ($?LS_SUBCWD) then

    # LSF has a list of processors already in a variable (LSB_HOSTS)
    # alias submit 'bsub < \!*'
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
          ./advance_model_batch.csh 0 ${NUM_ENS} filter_control00000 || exit 9

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
    # alias submit 'qsub \!*'
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
          ./advance_model.csh 0 ${NUM_ENS} filter_control00000 || exit 9

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
          ./advance_model.csh 0 ${NUM_ENS} filter_control00000 || exit 9

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

# original collaboration information - do not edit
# $orgURL: https://svn2.assembla.com/svn/ngdart/shell_scripts/run_filter_and_exit.csh $
# $orgId: run_filter_and_exit.csh 111 2010-06-09 21:55:44Z thoar $
# $orgRevision: 111 $
# $orgDate: 2010-06-09 15:55:44 -0600 (Wed, 09 Jun 2010) $

