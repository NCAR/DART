#!/bin/tcsh
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Wrapper interface between the DART programs perfect_model_obs and
# filter and the advance_group.csh script that handles advancing the
# COAMPS model.
#
# This script is the one that should be used in the adv_ens_command
# entry in the DART input.nml file, since DART insists on adding 
# command line arguments to that command.  If running on a system or 
# setup where you don't need to handle the ensemble integration in a 
# separately submitted job, you don't need to use this script.
#
# This script takes in command line arguments and writes them to disk
# for retrieval by the advance_group.csh script (which does the real
# work), creates a lock file, and checks that lock before terminating
# (since DART assumes that when this program is done, the ensemble has
# been integrated).
#
# The command line arguments specify the parent process number, the
# number of ensemble states to integrate, and the control file that
# gives instructions on which ensemble states to integrate.
######

set log_to_disk = yes

if ($log_to_disk == "yes") then
    set diskfile = "./wrapper_log"
    rm -f $diskfile
    date > $diskfile
endif

set parent       = $1
set num_states   = $2
set control_file = $3

if ($log_to_disk == "yes") then 
    echo "Wrapper to advance_group.csh..."     >> $diskfile 
    echo "    Parent is $parent"               >> $diskfile
    echo "    Number of states is $num_states" >> $diskfile
    echo "    Control file is $control_file"   >> $diskfile
endif

set WRAPPERFILE = groupwrapper

if ($log_to_disk == "yes") then
    echo "Writing arguments to wrapper file..." >> $diskfile
endif

rm -f ${WRAPPERFILE}
echo $parent       >  ${WRAPPERFILE}
echo $num_states   >> ${WRAPPERFILE}
echo $control_file >> ${WRAPPERFILE}

if ($log_to_disk == "yes") then
    echo "Creating ensemble lock file..." >> $diskfile
endif

set LOCKFILE = ens_running.lock
rm -f ${LOCKFILE}
touch ${LOCKFILE}

# Assume that advance_wrapper.csh and advance_group.csh are in the same
# directory - if this is not the case, change the call to qsub
if ($log_to_disk == "yes") then
    echo "Submitting job to PBS..." >> $diskfile
endif
qsub ./advance_group.csh

if ($log_to_disk == "yes") then
    qstat >> $diskfile
endif

# We need to put a block in here somehow so this doesn't quit until
# advance_group is done.  Need to do this with a "file exist" loop
# since named pipes don't work across machines.
echo "Waiting on lock file..."
while (-e ${LOCKFILE})
    sleep 5
end
echo "advance_wrapper done!"

exit 0


