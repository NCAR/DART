#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# check_tests.csh  run this after running test_dart.csh and it will
#                  look for errors in the build or run logs.
#                  some errors are expected (requires special files or
#                  libraries, or is intentionally provoking an error).


# if your system supports different options or needs to
# use a different location for these commands, set them here.
# they will be inherited by the other test scripts.
setenv REMOVE 'rm -f'
setenv RMDIR  'rmdir'
setenv COPY   'cp -p'
setenv MOVE   'mv -f'

# require we start running this from the developer_tests dir
if ( ! -d ../models ) then
   echo "../models directory does not exist. $0 must be run from"
   echo "the developer_tests directory."
   exit 2
endif

# number of lines of context to include
setenv nlines 4

# cd to the top level DART directory and 
# record where we are running this script
cd ..
set DARTHOME = `pwd`

if ( ! $?host) then
   setenv host `uname -n`
endif
echo "Running $0 on $host"
echo "The top-level DART directory is $DARTHOME"

# warning - the 'sampling_error_correction' file will match the simple
# string 'error' so you have to be more creative to find either 
# compile-time or run-time errors.

#----------------------------------------------------------------------
# setup complete
#----------------------------------------------------------------------

echo
echo
echo "=================================================================="
echo "Checking supported models"
echo "=================================================================="
echo
echo

#fgrep -C $nlines ERROR ${DARTHOME}/models/testing_logs/*
fgrep -C $nlines -i " error " ${DARTHOME}/models/testing_logs/*


#----------------------------------------------------------------------

echo
echo
echo "=================================================================="
echo "Checking observation converters"
echo "=================================================================="
echo
echo

echo "Not all observation converters are expected to build; you may"
echo "not have all the necessary supporting libraries.  So errors here"
echo "are not fatal."
echo ""

#fgrep -C $nlines ERROR ${DARTHOME}/observations/obs_converters/testing_logs/*
fgrep -C $nlines -i " error " ${DARTHOME}/observations/obs_converters/testing_logs/*


#----------------------------------------------------------------------

echo
echo
echo "=================================================================="
echo "Checking support programs"
echo "=================================================================="
echo
echo

#fgrep -C $nlines ERROR ${DARTHOME}/assimilation_code/programs/testing_logs/*
fgrep -C $nlines " error " ${DARTHOME}/assimilation_code/programs/testing_logs/*


#----------------------------------------------------------------------

echo
echo
echo "=================================================================="
echo "Checking developer tests" 
echo "=================================================================="
echo
echo

echo "A few developer tests are expected to fail."

#fgrep -C $nlines ERROR ${DARTHOME}/developer_tests/testing_logs/*
fgrep -C $nlines -i " error " ${DARTHOME}/developer_tests/testing_logs/*


#----------------------------------------------------------------------

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

