#!/bin/csh
# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

# This script assumes that an appropriate diag_table,
# RESTART directory, and 
# standard b-grid input.nml namelist file are
# available in the directory that is the parent of a temporary
# directory in which this script is executed.
# This script copies those files into the temporary directory
# and then executes the fortran program integrate_model
# which must have been compiled with the bgrid_model
# modules. 

# Might also want to edit in another directory as source for
# these.

cp -r ../RESTART .
cp ../diag_table ../input.nml .
../integrate_model 
