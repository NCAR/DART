#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$

#----------------------------------------------------------------------
# 'preprocess' is a program that culls the appropriate sections of the
# observation module for the observations types in 'input.nml'; the 
# resulting source file is used by all the remaining programs, 
# so this MUST be run first.
#----------------------------------------------------------------------

\rm -f preprocess *.o *.mod
\rm -f ../../../obs_def/obs_def_mod.f90
\rm -f ../../../obs_kind/obs_kind_mod.f90

set MODEL = "template"

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

\rm -f *.o *.mod input.nml*_default

if ( $#argv == 1 && "$1" == "-mpi" ) then
  echo "Success: All single task DART programs compiled."  
  echo "Script now compiling MPI parallel versions of the DART programs."
else if ( $#argv == 1 && "$1" == "-nompi" ) then
  echo "Success: All single task DART programs compiled."  
  echo "Script is exiting without building the MPI version of the DART programs."
  exit 0
else
  echo ""
  echo "Success: All DART programs compiled."
  echo "Script is exiting before building the MPI version of the DART programs."
  echo "Run the quickbuild.csh script with a -mpi argument or"
  echo "edit the quickbuild.csh script and remove the exit line"
  echo "to compile with MPI to run in parallel on multiple cpus."
  echo ""
  exit 0
endif

#----------------------------------------------------------------------
# to enable an MPI parallel version of filter for this model, 
# call this script with the -mpi argument, or if you are going to build
# with MPI all the time, remove or comment out the entire section above.
#----------------------------------------------------------------------

\rm -f filter wakeup_filter

@ n = $n + 1
echo
echo "---------------------------------------------------"
echo "build number $n is mkmf_filter"
csh   mkmf_filter -mpi
make

if ($status != 0) then
   echo
   echo "If this died in mpi_utilities_mod, see code comment"
   echo "in mpi_utilities_mod.f90 starting with 'BUILD TIP' "
   echo
   exit $n
endif

@ n = $n + 1
echo
echo "---------------------------------------------------"
echo "build number $n is mkmf_wakeup_filter"
csh  mkmf_wakeup_filter -mpi
make || exit $n

\rm -f *.o *.mod input.nml*_default

echo
echo 'time to run filter here:'
echo ' for lsf run "bsub < runme_filter"'
echo ' for pbs run "qsub runme_filter"'
echo ' for lam-mpi run "lamboot" once, then "runme_filter"'
echo ' for mpich run "mpd" once, then "runme_filter"'

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

