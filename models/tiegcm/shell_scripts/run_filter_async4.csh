#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
##------------------------------------------------------------------------------
## This block of directives constitutes the preamble for the LSF queuing system
## LSF is used on the IBM   Linux cluster 'lightning'
## LSF is used on the IMAGe Linux cluster 'coral'
## LSF is used on the IBM   'bluevista'
## The queues on lightning and bluevista are supposed to be similar.
##
## the normal way to submit to the queue is:    bsub < run_filter.csh
##
## an explanation of the most common directives follows:
## -J Job name (master script job.csh presumes filter_server.xxxx.log)
## -o STDOUT filename
## -e STDERR filename
## -P      account
## -q queue    cheapest == [standby, economy, (regular,debug), premium] == $$$$
## -n number of processors  (really)
##
#BSUB -J filter
#BSUB -o filter.%J.log
#BSUB -P P3507xxxx
#BSUB -q debug
#BSUB -n 4
#BSUB -W 0:30
#BSUB -N -u ${USER}@ucar.edu
#
##------------------------------------------------------------------------------
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
##
#PBS -N filter
#PBS -r n
#PBS -e filter.err
#PBS -o filter.log
#PBS -q dedicated
#PBS -l nodes=10:ppn=2

#===============================================================================

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

# FIXME ... read the async value from input.nml and set parallel_model accordingly.
# if async=2, e.g. you are going to run './modelxxx', single process
# (or possibly 'mpirun -np 1 ./modelxxx'), so each processor advances
# one ensemble independently of the others, leave this as false.
#
# if async=4, e.g. all the processors advance each modelxxx in turn with
# mpirun -np 64 modelxxx (or whatever) for as many ensembles as you have,
# set this to "true"

set parallel_model = "true"

# A common strategy for the beginning is to check for the existence of
# some variables that get set by the different queuing mechanisms.
# This way, we know which queuing mechanism we are working with,
# and can set 'queue-independent' variables for use for the remainder
# of the script.

if ($?LS_SUBCWD) then

    # LSF has a list of processors already in a variable (LSB_HOSTS)
    # alias submit 'bsub < \!*'
    echo "LSF - using mpirun.lsf for execution"
    setenv MPICMD mpirun.lsf

else if ($?PBS_O_WORKDIR) then

    # PBS has a list of processors in a file whose name is (PBS_NODEFILE)
    # alias submit 'qsub \!*'
    echo "PBS - using mpirun for execution"
    setenv MPICMD mpirun

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

endif

#-------------------------------------------------------------------------------
# Everything below this separator should not need to be modified if everything
# above the separator is set correctly.
#-------------------------------------------------------------------------------

if ( "$parallel_model" == "false" ) then

   # each filter task advances the ensembles, each running on 1 proc.

   ${MPICMD} ./filter

else

   # filter runs in parallel until time to do a model advance,
   # and then this script starts up the modelxxx jobs, each one
   # running in parallel. then it runs wakeup_filter to wake
   # up filter so it can continue. The communication happens through
   # 'named pipes' created by the mkfifo command.

   \rm -f model_to_filter.lock filter_to_model.lock
   mkfifo model_to_filter.lock filter_to_model.lock

   set filterhome = ~/.filter$$
   if ( ! -e $filterhome) mkdir $filterhome

   # start filter and immediately return control back to this script

   (setenv HOME $filterhome; ${MPICMD} ./filter) &

   while ( -e filter_to_model.lock )

      set todo=`cat < filter_to_model.lock`
      echo "todo received, value = ${todo}"

      if ( "${todo}" == "finished" ) then
         echo "main script: filter done."
         wait
         break

      else if ( "${todo}" == "advance" ) then

         # FIXME : in input.nml, the advance model command must
         # have -np N with N equal to the number of processors this job is using.

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

   if ( -d $filterhome) \rmdir $filterhome

endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

