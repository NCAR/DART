#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# Script to manage the compilation of all components for this model;

#----------------------------------------------------------------------
# Extra checks for missing files.
if ( ! -e ../get_name_info.f90 || ! -e ../coamps_intrinsic_mod.f90 ) then
 echo 'ERROR:'
 echo '../get_name_info.f90 and/or ../coamps_intrinsic_mod.f90 not found'
 echo ''
 echo 'To build the COAMPS/DART code, you are required to download the COAMPS'
 echo 'model separately.  See xxx for details.  This script can be run once'
 echo 'those files are in the $DART/models/coamps directory.'
 echo ''
 exit 1
endif
if ( ! -e input.nml ) then
 echo 'ERROR: input.nml not found'
 echo ''
 echo 'Copy ../templates/EXAMPLE.input.nml to ./input.nml (in the work dir)'
 echo 'and make any necessary changes.'
 echo ''
 exit 1
endif
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# 'preprocess' is a program that culls the appropriate sections of the
# observation module for the observations types in 'input.nml'; the
# resulting source file is used by all the remaining programs,
# so this MUST be run first.
#----------------------------------------------------------------------


\rm -f preprocess *.o *.mod Makefile .cppdefs

set MODEL = "coamps"

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

\rm -f *.o *.mod filter wakeup_filter Makefile .cppdefs

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

\rm -f *.o *.mod Makefile .cppdefs

exit 0


