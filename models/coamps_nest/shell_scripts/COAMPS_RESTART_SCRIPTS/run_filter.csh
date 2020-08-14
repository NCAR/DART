#!/bin/tcsh 
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
######
#
# SCRIPT:	run_filter.csh
# AUTHOR:	DART Folks
#
# N.B. I have adjusted this to use the advance_wrapper script instead
#      of the advance_model script as was originally done.  Otherwise,
#      I haven't messed with anything here.
#
######
#
# Modified for PBS queue on ACESGrid by TRW
##=============================================================================
## This block of directives constitutes the preamble for the PBS queuing system
## 
## the normal way to submit to the queue is:    qsub runme_filter
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
##                     them both.  (ppn == Processors Per Node)
##=============================================================================
#PBS -N filter
#PBS -r n
#PBS -e filter.err
#PBS -o filter.out
#PBS -q long
#PBS -l nodes=4:ppn=2

# Import the job-specific resource commands
source ${PBS_O_WORKDIR}/job_setup.csh 

# Set how many processors we will use for filter
set filterprocs = 8

# if async=2, e.g. you are going to run './wrf.exe', single process
# (or possibly 'mpirun -np 1 ./wrf.exe'), so each processor advances 
# one ensemble independently of the others, leave this as false.
#
# if async=4, e.g.  all the processors advance each wrf.exe in turn with
# mpirun -np 64 wrf.exe (or whatever) for as many ensembles as you have,
# set this to "true"

# if async=4, also check that the call to advance_model.csh
# has the right number of ensemble members below; it must match
# the input.nml number.

set parallel_model = "true"
set num_ens = 16

cd $PBS_O_WORKDIR

# Set us up to do some logging
set LOGFILE = "run_filter.log"
if (-e ${LOGFILE}) then
    rm -f ${LOGFILE}
    date > ${LOGFILE}
endif

# each filter task advances the ensembles, each running on 1 proc.
if ( "$parallel_model" == "false" ) then
     mpirun ./filter
else
    # filter runs in parallel until time to do a model advance,
    # and then this script starts up the wrf.exe jobs, each one
    # running in parallel.  then it runs wakeup_filter to wake
    # up filter so it can continue.
   
    # Reset the pipes - need to keep these names the same as in the
	# filter program itself.
    echo "Making pipes" >> ${LOGFILE}
    rm -f model_to_filter.lock filter_to_model.lock
    mkfifo model_to_filter.lock filter_to_model.lock

    set filterhome = `mktemp -d -p .`
    echo "Using temporary directory ${filterhome}" >> ${LOGFILE}
     
    # this starts filter but also returns control back to
    # this script immediately.
      
    set MPI_CMD = "mpiexec -n $filterprocs"

    echo "Running filter" >> ${LOGFILE}
    (setenv HOME $filterhome; ${MPI_CMD} ./filter | tee filter.dump)  &
       
    while ( -e filter_to_model.lock )
	echo "Waiting for input" >> ${LOGFILE}
        set todo=`( echo $< ) < filter_to_model.lock`
	echo todo received, value = ${todo}
	echo todo received, value = ${todo} >> ${LOGFILE}

    if ( "${todo}" == "finished" ) then
         echo main script: filter done.
         wait
         break                                
    else if ( "${todo}" == "advance" ) then
         # the second number below must match the number
         # of ensembles.  and in input.nml, the advance model
         # command must have -np N with N equal to the number
         # of processors this job is using.
	     echo "Advancing the model..." >> ${LOGFILE}
         echo calling model advance now:                
         ./advance_wrapper.csh 0 ${num_ens} filter_control00000  true
         
         echo "Model is done: now wake up filter..." >> ${LOGFILE}
         echo "Need to wake up ${filterprocs} processes" >> ${LOGFILE}
         echo "Using PBS nodefile ${PBS_NODEFILE}" >> ${LOGFILE}
         cat ${PBS_NODEFILE} >> ${LOGFILE}

         echo restarting filter.     
         mpirun -machinefile ${PBS_NODEFILE} -np $filterprocs ./wakeup_filter
    else
         echo main script: unexpected value received.
         break
    endif
       
    end
     
    echo "filter is finished." >> ${LOGFILE}
    echo filter finished, removing pipes.
    rm -f model_to_filter.lock filter_to_model.lock
    rm -f ${LOGFILE}
    if ( -d $filterhome) rmdir $filterhome
endif

exit 0


