#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# Shell script to work with *asynchronous* filter integration
# This script needs to be piped to the filter program with the
# filter namelist async variable set to .true.
#
# $Id$

# If this is first of recursive calls need to get rid of async_may_go
# Technically, this could lock, but seems incredibly unlikely
if ($?Notfirst) then
#  Do nothing if this is not the first time
# Actually, move data-model diagnostic file (extrapolations, data mismatches)
#    to a permanent storage
else
   setenv Notfirst yes
   rm -f async_may_go
# Clean up any assim_model_ic and ud files and temp directories
#   rm -f assim_model_ic* assim_model_ud*
    rm -f times batchflag .async_garb .batch_garb
# Call the model's initialization script to allow it to set up if needed
#   csh ./init_advance_model.csh
endif

while(1 == 1)
#   rm -f .async_garb
   ls async_may_go > .async_garb
   if ($status == 0) break
   echo 'waiting for async_may_go file'
   sleep 1
end

# create file to signal status of batch execution of ensemble
# cleanupgrade; machine specific. move elsewhere?  qsub ?
echo 'batch not done' > batchflag

# batch execution of ensemble
# qsub advance_ens.csh
./advance_ens.csh

# Wait for it to finish
while(1 == 1)
   ls batchflag >! .batch_garb
   if ($status != 0) break
   sleep 30
end
echo 'batch is done'

# Remove the semaphore file
rm -f async_may_go

# Cleaned up; let the filter know it can proceed
echo "All_done:Please_proceed"

# Doing recursive call
csh ./async_filter.csh
