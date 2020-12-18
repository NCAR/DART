#!/bin/tcsh
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Script to handle integrating a group of ensemble members specified
# by the DART filter program.  This is flexible enough to handle either
# being called with a subset of the ensemble members (i.e. when async
# in the DART namelist is set to 2) or with the entire ensemble at once
# (i.e. async = 4).
#
# Instead of using command line arguments, this script uses a wrapper file
# (called 'groupwrapper') that contains the parent process, the number of
# states to integrate, and the control file name that contains information
# on which members to advance.  This is done so this can be run as a new
# job within PBS if necessary (since qsub does not like arguments)
#
# N.B. The sed script of the form "1p $filename" is used several times 
#      in this script and just pulls a specific line (line 1 in the 
#      example) out of $filename
#
######
# PBS instructions for MIT ACESGrid cluster
######
#PBS -l nodes=1:ppn=2
#PBS -q long
#PBS -N coamps_ensemble
######

echo "Entering advance_group.csh"

# Number of processors used for a model run
set NPROCS = 1

# Decide if we're running an interactive job or not
# Modify this to allow for the possibility that it just doesn't exist!
if ($?PBS_ENVIRONMENT) then
    switch ($PBS_ENVIRONMENT)
    case PBS_BATCH:
    echo "Batch mode detected."
    set WORKDIR = $PBS_O_WORKDIR
    breaksw
    case PBS_INTERACTIVE:
    echo "Interactive mode detected."
    set WORKDIR = `pwd`
    breaksw
    default:
        echo "$PBS_ENVIRONMENT found, but not set!"
        exit (1)
    breaksw
    endsw
else
    set WORKDIR = `pwd`
    echo "PBS_ENVIRONMENT is not defined, setting workdir to here"
endif
cd ${WORKDIR}

# Grab the command line arguments - this is the number of the calling
# process, the number of states to integrate, and the control file to
# use for the integration
set WRAPPERFILE = ${WORKDIR}/groupwrapper
echo "Using wrapper file $WRAPPERFILE"
set parent_process = `sed -ne "1p" $WRAPPERFILE`
set num_states     = `sed -ne "2p" $WRAPPERFILE`
set control_file   = `sed -ne "3p" $WRAPPERFILE`

# Output to confirm job characteristics
echo "Running on host "`hostname`
echo "Time is "`date`
echo "Directory is "`pwd`
echo "Working directory is:" $WORKDIR
echo "Parent process is:" $parent_process
echo "Number of states to integrate is:" $num_states
echo "Control file for this integration is:" $control_file

### Set up parallel ensemble handling

# Create a named pipe to handle the communication - randomize this so
# we don't end up with bad things happening between runs 
set PIPENAME = `mktemp -u $WORKDIR/ensemble_pipe.XXXXXXX`
echo "Using pipe $PIPENAME"
mkfifo $PIPENAME

# Parallel handling parameters are the maximum number of concurrent 
# runs we'll allow, the number of ensemble members currently running,
# and the number of ensemble members that have finished running.
set MAX_RUNNING  = 2
set num_running  = 0
set num_finished = 0

# Initialize the "member completed" array
set member_status = (`seq 1 $num_states`)
foreach ii (`seq 1 $num_states`)
    set member_status[$ii] = NOTRAN
end

# Initialize the "ensemble member" array - this is the same size as
# the number of states we are integrating, but contains the specific
# ensemble members.  I realize that this loop could be combined with
# the one just prior, but I think this separation is clearer.
set member_number = (`seq 1 $num_states`)
set curline = 1
foreach ii (`seq 1 $num_states`)
    set member_number[$ii] = `sed -ne "${curline}p" $control_file`
    @ curline = $curline + 3
end

# Flush out the progress bar
rm -f $WORKDIR/progress

# Start the parallel processing loop - do this as long as we still
# have things to run
while ($num_finished < $num_states)
    echo "Current states finished: $num_finished out of $num_states" \
    >> progress
    while ($num_running < $MAX_RUNNING)
    @ num_processing = $num_running + $num_finished

    # Only spawn a new process if doing so won't mean we have more
    # processes than members
    if ($num_processing >= $num_states) then
        break
    endif

    # If we can run another member, find the next one to run
    set ii = 1
    while ($member_status[$ii] == RUNNING || \
               $member_status[$ii] == RAN)
           @ ii++
        end
    echo "  Found member $ii to run next" >> progress
    set member_status[$ii] = RUNNING

    # Actually run the model
    set member = $member_number[$ii]
    set memberlog = `mktemp "member${member}.XXXXXX"`
    $WORKDIR/advance_model.csh $WORKDIR $member $NPROCS $PIPENAME \
        >& ${memberlog} &
    @ num_running++
    end

    # Now that we've spawned as many model runs as we can, listen for
    # them to finish via the named pipe
    foreach donemem (`cat < $PIPENAME`)
    set donemem_number = `echo "${donemem}" | gawk -F. '{print $1}'`
    set donemem_status = `echo "${donemem}" | gawk -F. '{print $2}'`

    set member_status[${donemem_number}] = ${donemem_status}
    if (${donemem_status} == RAN) then
        @ num_finished++
        echo "Member ${donemem} is done" >> progress
    else
        echo "Error for member ${donemem}!" >> progress
    endif

    @ num_running--
    end
end

echo "Done with running the group: exiting..." >> progress

# signal to async_filter.csh to continue
rm -f $control_file

# Clean up the pipe
rm -rf $PIPENAME

# Now, if the group wrapper trigger file exists, let the wrapper know
# that we're finished by deleting it
set GROUPLOCK = ens_running.lock
if (-e ${GROUPLOCK}) then
    rm -f ${GROUPLOCK}
    echo "Leaving advance_group.csh"
endif

exit 0


