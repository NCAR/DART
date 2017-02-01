#!/usr/bin/env python
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# CREDIT: This script was donated to DART by Luke Madaus during his time
# at the University of Washington. Thanks Luke!
#
# Modified for Python (Oct. 2015) Luke Madaus, University of Washington
#-----------------------------------------------------------------------------
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
# Adapted for Python (Luke Madaus, Oct 2015)
# Python imports go here
from __future__ import print_function, division
import os, time
from datetime import datetime, timedelta
from netCDF4 import Dataset
from namelist_utils import read_namelist, write_namelist



# if async=2, e.g. you are going to run './modelxxx', single process
# (or possibly 'mpirun -np 1 ./modelxxx'), so each processor advances
# one ensemble independently of the others, leave this as false.
#
# if async=4, e.g. all the processors advance each modelxxx in turn with
# mpirun -np 64 modelxxx (or whatever) for as many ensembles as you have,
# set this to "true"

# this script is going to determine several things by reading the input.nml
# file which contains the &filter_nml namelist.  make sure it exists first.
#if ( ! -e input.nml ) then
if not os.path.exists('input.nml'):
    print("ERROR - input.nml does not exist in local directory.")
    print("ERROR - input.nml needed to determine several settings for this script.")
    exit(1)


# detect whether the model is supposed to run as an MPI job or not
# by reading the "async = " from the &filter_nml namelist in input.nml.
# some namelists contain the same string - be sure to get the filter_nml one
# by grepping for lines which follow it.
nmld = read_namelist('input.nml')
async_type = nmld['filter_nml']['async']

if (async_type == 0) or (async_type == 2):
    parallel_model = False
elif (async_type == 4):
    parallel_model = True
else: 
    print("cannot autodetect async value in the filter_nml namelist in input.nml file.")
    print("hardcode the parallel_model shell variable and comment out these lines.")
    exit(-1)
    parallel_model = False

# Determine the number of ensemble members from input.nml,
# as well as the command for advancing the model.
num_ens = nmld['filter_nml']['ens_size']
adv_cmd = nmld['filter_nml']['adv_ens_command']


# A common strategy for the beginning is to check for the existence of
# some variables that get set by the different queuing mechanisms.
# This way, we know which queuing mechanism we are working with,
# and can set 'queue-independent' variables for use for the remainder
# of the script.

if (os.environ.get('LS_SUBCWD') not in ['', None]):

    # LSF has a list of processors already in a variable (LSB_HOSTS)
    print("LSF - using mpirun.lsf for execution")

    # each filter task advances the ensembles, each running on 1 proc.
    if not parallel_model:
        os.system('mpirun.lsf ./filter')

    else:

        # filter runs in parallel until time to do a model advance,
        # and then this script starts up the modelxxx jobs, each one
        # running in parallel. then it runs wakeup_filter to wake
        # up filter so it can continue.

        os.system('rm -f model_to_filter.lock filter_to_model.lock')
        os.system('mkfifo model_to_filter.lock filter_to_model.lock')

        curpid = os.getpid()
        filterhome = '~/.filter{:d}'.format(curpid)
        if not os.path.exists(filterhome): os.system('mkdir {:s}'.format(filterhome))

        # this starts filter but also returns control back to
        # this script immediately.
        os.envrion.put("HOME", filterhome)
        os.system('mpirun.lsf ./filter &')

        while os.path.exists('filter_to_model.lock'):
            with open('filter_to_model.lock','r') as todofile:
                todo = todofile.read()
                print("todo received, value ={:s}".format(todo))

                if ("finished" in todo):
                    print("main script: filter done.")
                    time.sleep(5)
                    break

                elif ("advance" in todo):

                    # the second number below must match the number
                    # of ensembles. Also, in input.nml, the advance model
                    # command must have -np N with N equal to the number
                    # of processors this job is using.

                    print("calling model advance now:")

                    os.system('{adv_cmd} 0 {num_ens} filter_control00000 || exit 9'.format(**locals()))

                    print("restarting filter.")
                    os.system('mpirun.lsf ./wakeup_filter')

                else:

                    print("main script: unexpected value received.")
                    break



      print("filter finished, removing pipes.")
      os.system('rm -f model_to_filter.lock filter_to_model.lock')

      if os.path.exists(filterhome): os.system('rmdir {:s}'.format(filterhome))


elif (os.environ.get(PBS_O_WORKDIR) not in ['', None]):

    # PBS has a list of processors in a file whose name is (PBS_NODEFILE)
    print("PBS - using mpirun for execution")

    # each filter task advances the ensembles, each running on 1 proc.
    if not parallel_model:
        os.system('mpirun ./filter')

    else:

        # filter runs in parallel until time to do a model advance,
        # and then this script starts up the modelxxx jobs, each one
        # running in parallel. then it runs wakeup_filter to wake
        # up filter so it can continue.

        os.system('rm -f model_to_filter.lock filter_to_model.lock')
        os.system('mkfifo model_to_filter.lock filter_to_model.lock')

        filterhome = '~/.filter'
        if not os.path.exists(filterhome): os.system('mkdir {:s}'.format(filterhome))

        # this starts filter but also returns control back to
        # this script immediately.
        os.environ.put('HOME', filterhome)  
        os.system('mpirun ./filter &')

        while os.path.exists('filter_to_model.lock')
            with open('filter_to_model.lock','r') as todofile:
                todo = todofile.read()
                print("todo received, value = {:s}".format(todo))

            if "finished" in todo:
                print("main script: filter done.")
                time.sleep(5)
                break

            elif "advance" in todo:

                # the second number below must match the number
                # of ensembles. Also, in input.nml, the advance model
                # command must have -np N with N equal to the number
                # of processors this job is using.

                print("calling model advance now:")
                os.system('{adv_cmd} 0 {num_ens} filter_control00000 || exit 9'.format(**locals()))

                print("restarting filter.")
                os.system('mpirun ./wakeup_filter')

            else:

                print("main script: unexpected value received.")
                break


        print("filter finished, removing pipes.")
        os.system('rm -f model_to_filter.lock filter_to_model.lock')

        if os.path.exists(filterhome): os.system('rmdir {:s}'.format(filterhome))

else:

    # If you have a linux cluster with no queuing software, use this
    # section. The list of computational nodes is given to the mpirun
    # command and it assigns them as they appear in the file. In some
    # cases it seems to be necessary to wrap the command in a small
    # script that changes to the current directory before running.

    print("running with no queueing system")

    # before running this script, do this once. the syntax is
    # node name : how many tasks you can run on it
    #setenv MYNODEFILE ~/nodelist
    #print("node7:2" >! $MYNODEFILE
    #print("node5:2" >> $MYNODEFILE
    #print("node3:2" >> $MYNODEFILE
    #print("node1:2" >> $MYNODEFILE

#   another possibility - note hardwired NP ...
    mpirun = '/share/apps/openmpi/gfortran/bin/mpirun'
    mpicmd = '$MPIRUN --hostfile nodelist-gfortran --mca mtl mx --mca pml cm -np 72'

    print("mpicmd = {:s}".format(mpicmd))

    # filter runs in parallel until time to do a model advance,
    # and then this script starts up the modelxxx jobs, each one
    # running in parallel. then it runs wakeup_filter to wake
    # up filter so it can continue.

    os.system('rm -f model_to_filter.lock filter_to_model.lock')
    os.system('mkfifo model_to_filter.lock filter_to_model.lock')

    filterhome = '~/.filter{:d}'.format(os.getpid())
    if not os.path.exists(filterhome): os.system('mkdir {:s}'.format(filterhome))

    # this starts filter but also returns control back to
    # this script immediately.

    os.environ.put('HOME', filterhome)
    os.system('{:s} ./filter) &'.format(mpicmd))

    while os.path.exists('filter_to_model.lock'):

        with open('filter_to_model.lock','r') as todofile:
            todo = todofile.read()
            print("todo received, value = {:s}".format(todo))

        if "finished" in todo:
            print("main script: filter done.")
            time.sleep(5)
            break
        elif "advance" in todo:

            # the second number below must match the number
            # of ensembles. Also, in input.nml, the advance model
            # command must have -np N with N equal to the number
            # of processors this job is using.

            print("calling model advance now:")
            os.system('{adv_cmd} 0 {num_ens} filter_control00000 || exit 9'.format(**locals()))

            print("restarting filter.")
            os.system('mpirun ./wakeup_filter')

        else:

            print("main script: unexpected value received.")
            break



    print("filter finished, removing pipes.")
    os.system('rm -f model_to_filter.lock filter_to_model.lock')

    if os.path.exists(filterhome): os.system('rm -rf {:s}'.format(filterhome))


exit(0)

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

