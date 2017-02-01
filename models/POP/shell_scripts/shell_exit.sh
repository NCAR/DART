#!/bin/sh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# LSF does not reliably return an exit code from
# a serial section of a batch script, only the return
# from the last call to mpirun.lsf.  so to make a
# batch script exit, call this:
# 
#  setenv LSB_PJL_TASK_GEOMETRY "{(0)}"
#  setenv EXITCODE -1
#  mpirun.lsf shell_exit.sh
 
if [[ $# -gt 0 ]]; then
   EXITCODE=$1
fi

if [[ "$EXITCODE" == "" ]]; then
   EXITCODE=0
fi

echo exiting with status code $EXITCODE
exit $EXITCODE
   
# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

