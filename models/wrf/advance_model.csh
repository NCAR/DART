#!/bin/csh
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

# This script copies the necessary files into the temporary directory
# and then executes the fortran program integrate_wrf.

# In construction.

../integrate_wrf
