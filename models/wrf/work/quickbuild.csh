#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This script compiles all executables in this directory.

#----------------------------------------------------------------------
# 'preprocess' is a program that culls the appropriate sections of the
# observation module for the observations types in 'input.nml'; the 
# resulting source file is used by all the remaining programs, 
# so this MUST be run first.
#----------------------------------------------------------------------

\rm -f preprocess *.o *.mod
\rm -f ../../../obs_def/obs_def_mod.f90
\rm -f ../../../obs_kind/obs_kind_mod.f90

set MODEL = "wrf"

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

\rm -f *.o *.mod 
\rm -f input.nml*_default

if ( $#argv == 1 && "$1" == "-mpi" ) then
  echo "Success: All single task DART programs compiled."  
  echo "Script now compiling MPI parallel versions of the DART programs."
else if ( $#argv == 1 && "$1" == "-nompi" ) then
  echo "Success: All single task DART programs compiled."  
  echo "Script is exiting without building the MPI version of the DART programs."
  exit 0
else
  echo ""
  echo "Success: All single task DART programs compiled."  
  echo "Script now compiling MPI parallel versions of the DART programs."
  echo "Run the quickbuild.csh script with a -nompi argument or"
  echo "edit the quickbuild.csh script and add an exit line"
  echo "to bypass compiling with MPI to run in parallel on multiple cpus."
  echo ""
endif

#----------------------------------------------------------------------
# to disable an MPI parallel version of filter for this model, 
# call this script with the -nompi argument, or if you are never going to
# build with MPI, add an exit before the entire section above.
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Build the MPI-enabled target(s) 
#----------------------------------------------------------------------

\rm -f *.o *.mod 

set MPI_ENABLED = "filter wakeup_filter verify inflate perturb"

foreach TARGET ( $MPI_ENABLED )

   \rm -f $TARGET

   set MPROG = mkmf_${TARGET}

   @ n = $n + 1
   echo
   echo "---------------------------------------------------"
   echo "${MODEL} build number ${n} is ${TARGET}" 
   \rm -f ${TARGET}
   csh $MPROG -mpi || exit $n
   make            || exit $n
end

\rm -f *.o *.mod 
\rm -f input.nml*_default

echo "All MPI-enabled targets built"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

