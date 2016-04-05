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

\rm -f preprocess *.o *.mod
\rm -f ../../../obs_def/obs_def_mod.f90
\rm -f ../../../obs_kind/obs_kind_mod.f90

set TEST = "utilities"

@ n = 1

echo
echo
echo "---------------------------------------------------------------"
echo "${TEST} build number ${n} is preprocess"

csh  mkmf_preprocess
make || exit $n

./preprocess || exit 99

#----------------------------------------------------------------------
# Build all targets
#----------------------------------------------------------------------

foreach TARGET ( mkmf_* )

   set PROG = `echo $TARGET | sed -e 's#mkmf_##'`

   switch ( $TARGET )
   case mkmf_preprocess:
      breaksw
   default:
      @ n = $n + 1
      echo
      echo "---------------------------------------------------"
      echo "${TEST} build number ${n} is ${PROG}"
      \rm -f ${PROG}
      csh $TARGET || exit $n
      make        || exit $n
      \rm -f *.o *.mod    # verify path_names - remove when working
      breaksw
   endsw
end

\rm -f *.o *.mod
\rm -f input.nml*_default


exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

