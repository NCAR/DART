#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# Shell script to work with *synchronous* filter integration.
# This script is called in Aadvance_state (from assim_model_mod)
# when async = 2.
#
# $Id$

if ($?Notfirst) then
#  Do nothing if this is not the first time
# Actually, move data-model diagnostic file (extrapolations, data mismatches)
#    to a permanent storage
else
   setenv Notfirst yes
# Clean up any assim_model_ic and ud files and temp directories
#   rm -f assim_model_ic* assim_model_ud*
    rm -f times batchflag .batch_garb
# Call the model's initialization script to allow it to set up if needed
#   csh ./init_advance_model.csh
endif

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
