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
# a (set of?) DART initial condition (i.e. restart file) that was created by
# run_nogapsIC_to_dart.csh.
#
# run_nogapsIC_to_dart.csh creates dart state vectors 'filter_ic.NNNN'
# in the experiment directory. Use any of these as input to perfect_model_obs:
# &perfect_model_obs_nml:restart_in_file_name = "filter_ic.0001"
#
# perfect_model_obs is inherently single-threaded, but the model advance can
# be anything.
#
#=============================================================================
#
# This block of directives constitutes the preamble for the LSF queuing system
#
# the normal way to submit to the queue is:    bsub < run_perfect.csh
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
#BSUB -J perfect
#BSUB -o perfect.%J.log
#BSUB -q standby
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
#PBS -N perfect
#PBS -r n
#PBS -e perfect.err
#PBS -o perfect.log
#PBS -q dedicated
#PBS -l nodes=4:ppn=2

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
   setenv JOBNAME     perfect
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $host

endif

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
cp $PERFECT_exec_name                         . || exit 1
cp $NOGAPS2DART_exec_name                     . || exit 1
cp $DART2NOGAPS_exec_name                     . || exit 1
cp $TRANSTIME_exec_name                       . || exit 1
cp ${DART_script_dir}/advance_model.csh       . || exit 1
cp ${NOGAPS_exec_dir}/${NOGAPS_exec_name}     . || exit 1

# for perfect model, this is the template of observations
# to fill in.
cp ${ORIGINALDIR}/obs_seq.in                  . || exit 1

# put in some bulletproofing here to be sure the files that
# the NOGAPS model_mod.f90 needs are indeed in the right place.
if ( ! -f noggeom.txt ) then
 cp -v ${climo}/noggeom${resolution}l${n_levels}_full.txt noggeom.txt
endif

# From here down, this script is pretty much like all the other
# models which work with DART.

# A common strategy for the beginning is to check for the existence of
# some variables that get set by the different queuing mechanisms.
# This way, we know which queuing mechanism we are working with,
# and can set 'queue-independent' variables for use for the remainder
# of the script.

if ($?LSB_QUEUE || $?PBS_QUEUE) then

    # Must be using LSF or PBS as the queueing system.
    # The model advance is launched using the MPI, but 

    echo "using ${MPI} for model executable advance"

    # start perfect_model_obs here, not an MPI job

    ./${PERFECT_exec_name:t}

else

    # WARNING: This block is untested with NOGAPS. It is unlikely to work.

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
    set    MPI = $MPIRUN -np $NUM_PROCS -nolocal -machinefile nodelist-pgi

#   another possibility - note hardwired NP ...
    set MPIRUN = /share/apps/openmpi/gfortran/bin/mpirun
    set    MPI = $MPIRUN --hostfile nodelist-gfortran --mca mtl mx --mca pml cm -np 72

    echo "MPI = ${MPI}"

    # start perfect_model_obs here, not an MPI job

    ./${PERFECT_exec_name:t}

endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

# original collaboration information - do not edit
# $orgURL: https://svn2.assembla.com/svn/ngdart/shell_scripts/run_perfect.csh $
# $orgId: run_perfect.csh 111 2010-06-09 21:55:44Z thoar $
# $orgRevision: 111 $
# $orgDate: 2010-06-09 15:55:44 -0600 (Wed, 09 Jun 2010) $

