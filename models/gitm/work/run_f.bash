#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

noe=$1 #number of ensemble members
nop=$2 #number of processors the filter can run on


rm hf #if there are hf's from old runs
head -n $nop $PBS_NODEFILE > hf #create the hostfile for filter, assumes we are in work dir
cat hf

    # filter runs in parallel until time to do a model advance,
    # and then this script starts up the modelxxx jobs, each one
    # running in parallel. then it runs wakeup_filter to wake
    # up filter so it can continue.

    \rm -f model_to_filter.lock filter_to_model.lock
    mkfifo model_to_filter.lock filter_to_model.lock

    export filterhome=~/.filter$$
    if [ ! -e $filterhome ] ; then
	mkdir $filterhome
    fi
    
    # this starts filter but also returns control back to
    # this script immediately.

    (export HOME=$filterhome; mpiexec -hostfile hf -n $nop ./filter || exit 41 ) & #run $nop instances of filter on hf
    while [ -e filter_to_model.lock ]
    do
        export todo=`cat < filter_to_model.lock`
        echo "todo received, value = ${todo}"

        if [ "${todo}" == "finished" ] ; then
          echo "main script: filter done."
          wait
          break

        elif [ "${todo}" == "advance" ] ; then

          # the second number below must match the number
          # of ensembles. Also, in input.nml, the advance model
          # command must have -np N with N equal to the number
          # of processors this job is using.

          echo "calling model advance now:"
          ../shell_scripts/advance_model.csh 0 $noe filter_control00000 || exit 9

          echo "restarting filter."
          mpiexec -hostfile hf -n $nop ./wakeup_filter 

        else

          echo "main script: unexpected value received."
          break

        fi

    done

    echo "filter finished, removing pipes."
    \rm -f model_to_filter.lock filter_to_model.lock

    if [ -d $filterhome ] ; then
	rm -r $filterhome
    fi
#!/bin/bash

noe=$1 #number of ensemble members
nop=$2 #number of processors the filter can run on




rm hf #if there are hf's from old runs
head -n $nop $PBS_NODEFILE > hf #create the hostfile for filter, assumes we are in work dir
cat hf

    # filter runs in parallel until time to do a model advance,
    # and then this script starts up the modelxxx jobs, each one
    # running in parallel. then it runs wakeup_filter to wake
    # up filter so it can continue.

    \rm -f model_to_filter.lock filter_to_model.lock
    mkfifo model_to_filter.lock filter_to_model.lock

    export filterhome=~/.filter$$
    if [ ! -e $filterhome ] ; then
	mkdir $filterhome
    fi
    
    # this starts filter but also returns control back to
    # this script immediately.

    (export HOME=$filterhome; mpiexec -hostfile hf -n $nop ./filter || exit 41 ) & #run $nop instances of filter on hf
    while [ -e filter_to_model.lock ]
    do
        export todo=`cat < filter_to_model.lock`
        echo "todo received, value = ${todo}"

        if [ "${todo}" == "finished" ] ; then
          echo "main script: filter done."
          wait
          break

        elif [ "${todo}" == "advance" ] ; then

          # the second number below must match the number
          # of ensembles. Also, in input.nml, the advance model
          # command must have -np N with N equal to the number
          # of processors this job is using.

          echo "calling model advance now:"
          ../shell_scripts/advance_model.csh 0 $noe filter_control00000 || exit 9

          echo "restarting filter."
          mpiexec -hostfile hf -n $nop ./wakeup_filter 

        else

          echo "main script: unexpected value received."
          break

        fi

    done

    echo "filter finished, removing pipes."
    \rm -f model_to_filter.lock filter_to_model.lock

    if [ -d $filterhome ] ; then
	rm -r $filterhome
    fi

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

