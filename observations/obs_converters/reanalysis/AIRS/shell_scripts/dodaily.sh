#!/bin/ksh

# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
#--------------------------------------------------------------
# DESCRIPTION:
#  Convert a day's worth of AIRS T/Q obs into a DART obs_seq file.
#
# this is the working script that expects to be called by a higher
# level script.  it should be given a working directory that already
# has an appropriate input.nml and linked executables in it.
#


# workdir
cd $1

# do the conversion
./convert_airs_L2

# remove the executables and what we expect to see.  
# if no other files are there, remove the dir.
rm -f convert_airs_L2
rm -f dart_log.out dart_log.nml flist input.nml

# will fail if dir not empty.
cd ..
rmdir $1

exit 0

