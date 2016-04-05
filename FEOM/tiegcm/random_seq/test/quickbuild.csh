#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Script to manage the compilation of all tests in this directory.

#----------------------------------------------------------------------

\rm -f *.o *.mod

set TEST = "random_seq"

#----------------------------------------------------------------------
# Build all targets
#----------------------------------------------------------------------

@ n = 0
foreach TARGET ( mkmf_* )

   set PROG = `echo $TARGET | sed -e 's#mkmf_##'`

   @ n = $n + 1
   echo
   echo "---------------------------------------------------"
   echo "${TEST} build number ${n} is ${PROG}" 
   \rm -f ${PROG}
   csh $TARGET || exit $n
   make        || exit $n

end

\rm -f *.o *.mod input.nml*_default

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

