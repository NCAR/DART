#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2008, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$

# Script to manage the compilation of the executables in this directory.
#
# The 'preprocess' step constructs 2 modules which define which DART 
# observation types will be compiled into the code.  To add (or remove)
# obs types, edit the 'input.nml' namelist file, and find the &preprocess_nml
# section.  The 'input_types' item is an array of character strings listing
# the observation definition files which will be included when preprocess
# is compiled and run.  Add and remove filenames from this list to control
# the observation types.
#
#----------------------------------------------------------------------

\rm -f preprocess *.o *.mod
\rm -f ../../../../obs_def/obs_def_mod.f90
\rm -f ../../../../obs_kind/obs_kind_mod.f90

set MODEL = "create_real_obs"

@ n = 1

echo
echo
echo "---------------------------------------------------------------"
echo "${MODEL} build number ${n} is preprocess"

csh  mkmf_preprocess
make || exit $n

./preprocess || exit 99

#----------------------------------------------------------------------
# Build all the single-threaded targets
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
      echo "${MODEL} build number ${n} is ${PROG}" 
      \rm -f ${PROG}
      csh $TARGET || exit $n
      make        || exit $n
      breaksw
   endsw
end

echo "Success: All DART programs compiled."  
exit 0

