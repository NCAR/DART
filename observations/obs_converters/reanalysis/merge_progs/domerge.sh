#!/bin/ksh

# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
#--------------------------------------------------------------
# DESCRIPTION:
# Merge NCEP BUFR obs (inc ACARS), GPS RO obs, and AIRS T/Q obs.
#
# this is the working script that expects to be called by a higher
# level script.  it should be given a working directory that already
# has an appropriate input.nml and linked executables in it,
# and the location of the output file directory.
#


# workdir
cd $1

# intended output dir
if [ ! -d $2 ]; then mkdir -p $2; fi

# do the merge
./obs_sequence_tool 
if [ $status == 0 ] then
   # remove the executables and what we expect to see.  
   rm -f advance_time obs_sequence_tool
   rm -f dart_log.out dart_log.nml olist input.nml
   # Otherwise, leave files that will cause the calling scripts to fail
   # because this directory still exists.
fi

# Try to remove the dir.  will fail if dir not empty.
cd ..
rmdir $1

exit 0

