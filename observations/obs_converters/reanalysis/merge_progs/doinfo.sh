#!/bin/ksh

# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
#--------------------------------------------------------------
# DESCRIPTION:
# run a months worth of obs_info to collect obs counts in the
# input obs_seq files for the reanalysis.
#
# this is the working script that expects to be called by a higher
# level script.  it should be given a working directory that already
# has an appropriate input.nml and linked executables in it,
# and the location of the output file directory.
#


# workdir
cd $1

# intended output file
if [ -f $2 ]; then rm -f $2; fi

# run the job
./obs_info

# move the output file up one dir for now
mv ocount.txt ../$2

# remove the executables and what we
# expect to see.  if nothing else is
# there, remove the dir.
rm -f obs_info flist 
rm -f dart_log.out dart_log.nml flist input.nml

# will fail if dir not empty.
cd ..
rmdir $1

exit 0

