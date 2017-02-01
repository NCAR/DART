#!/bin/sh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# some versions of LSF do not reliably return an exit code from
# a serial section of a batch script, only the return from the 
# last call to mpirun.lsf.  so to make a batch script really exit
# in case of error, call this:
# 
#  setenv LSB_PJL_TASK_GEOMETRY "{(0)}"   # optional
#  setenv EXITCODE -1
#  mpirun.lsf shell_exit.sh
 
# if you call this inline with an arg:
if [[ $# -gt 0 ]]; then
   export EXITCODE=$1
fi

# if you set an env var (the normal way to pass in
# a value if using mpi):
if [[ "$EXITCODE" == "" ]]; then
   echo "incoming EXITCODE was not set, setting to 0"
   export EXITCODE=0
fi

echo exiting with status code $EXITCODE
exit $EXITCODE
   
# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

