#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
# 
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$
#
## || marks changes to make the script run CAM in parallel (and batch)
## 'sl' marks changes to make the script run on /scratch/local on each node
#
# Shell script to work with *syncronous* filter integration       ?
# This script needs to be piped to the filter program with the
# filter namelist async variable set to .true.

# set echo verbose

# If this is first of recursive calls need to get rid of async_may_go
# Technically, this could lock, but seems incredibly unlikely
if ($?First) then
#  Do nothing if this is not the first time
   echo 'doing nothing in this call to sync_submit.csh' >> dump
# Actually, move data-model diagnostic file (extrapolations, data mismatches)
#    to a permanent storage
else
   setenv First no
   echo 'setting First = ' $First ' and removing assim_model_ics' >> dump
###   rm -f async_may_go
# Clean up any assim_model_ic and ud files and temp directories
   rm -f assim_model_ic*
   rm -f assim_model_ud*
   rm -f times batchflag .*garb
#sl    rm -rf tempdir*
# Call the model's initialization script to allow it to set up if needed
   csh ./init_advance_model.csh
endif

echo 'before async_may_go section'
while(1 == 1)
   rm -f .async_garb
   ls async_may_go > .async_garb
   if ($status == 0) break
   echo 'waiting_for_async_may_go_file'
   sleep 15
end
echo 'found_async_may_go_file'

# create file to signal status of batch execution of ensemble
echo 'batch not done' > batchflag

# batch execution of ensemble
qsub advance_ens.csh

# Wait for it to finish
echo waiting_for_batch_advance_ens
while(1 == 1)
   ls batchflag >! .batch_garb
   if ($status != 0) break
   sleep 30
end
echo batch_is_done

# Remove the semaphore file
rm -f async_may_go

# Cleaned up; let the filter know it can proceed
echo All_done:Please_proceed

# Doing recursive call
# || csh ./async_filter.csh
csh ./sync_submit.csh
