#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#-------------------------------------------------------------------------------
# compile all programs in the current directory that have a mkmf_xxx file.
#
# usage: [ -mpi | -nompi ]
#
# environment variable options (which must be set before running this code):
#
#  "setenv CODE_DEBUG 1" (csh) or "export CODE_DEBUG=1" (bash)
#  to keep the .o and .mod files in the current directory instead of 
#  removing them at the end.  this usually improves runtime error reports 
#  and these files are required by most debuggers. This is a bit complicated
#  if you are trying to debug serial targets but also want to build MPI targets.
#  By virtue of building MPI targets, you must destroy the serial .o and .mod
#  files. So - to debug serial targets, you must use the -nompi option.
#
#  "setenv DART_TEST 1" (csh) or "export DART_TEST=1" (bash)
#  will remove all .o and .mod files between each build.
#  This is useful to ensure that the path_names_* files are sufficient.
#
#-------------------------------------------------------------------------------

set MODEL = "utilities test"

# programs which have the option of building with MPI:
set MPI_TARGETS = "error_handler_test"

# set default (override with -mpi or -nompi):
#  0 = build without MPI, 1 = build with MPI
set with_mpi = 1

#-------------------------------------------------------------------------------
# shouldn't have to modify this script below here.

if ( $#argv >= 1 ) then
   if ( "$1" == "-mpi" ) then
      set with_mpi = 1 
   else if ( "$1" == "-nompi" ) then
      set with_mpi = 0
   else
      echo usage: $0 '[ -mpi | -nompi ]'
      exit 0
   endif
endif

set preprocess_done = 0
set tdebug = 0
set cdebug = 0

if ( $?CODE_DEBUG && $?DART_TEST ) then
   echo ''
   echo "ERROR: Both environment variables CODE_DEBUG and DART_TEST are set."
   echo "ERROR: This is a contradictory situation."
   echo "ERROR: CODE_DEBUG indicates the .o and .mod files are to be kept."
   echo "ERROR: DART_TEST  indicates the .o and .mod files are to be removed."
   echo ''
   exit 1
endif

if ( $?CODE_DEBUG ) then
   set cdebug = $CODE_DEBUG
endif
if ( $?DART_TEST ) then
   set tdebug = $DART_TEST
endif

set nonomatch
\rm -f *.o *.mod Makefile .cppdefs preprocess

#-------------------------------------------------------------------------------
# Build all the single-threaded targets
# 'preprocess' is a program that culls the appropriate sections of the
# observation module for the observations types in 'input.nml'; the
# resulting source file is used by all the remaining programs,
# so this MUST be run first.
#-------------------------------------------------------------------------------
@ n = 0

foreach TARGET ( mkmf_preprocess mkmf_* )

   set PROG = `echo $TARGET | sed -e 's#mkmf_##'`

   if ( $PROG == "preprocess" && $preprocess_done ) goto skip

   if ( $with_mpi ) then
      foreach i ( $MPI_TARGETS )
         if ( $PROG == $i ) goto skip
      end
   endif

   @ n = $n + 1
   echo
   echo "---------------------------------------------------"
   echo "${MODEL} build number ${n} is ${PROG}"
   \rm -f ${PROG}
   csh $TARGET || exit $n
   make        || exit $n

   if ( $tdebug ) then
      echo 'removing all files between builds'
      \rm -f *.o *.mod Makefile .cppdefs
   endif

   # preprocess creates module files that are required by
   # the rest of the executables, so it must be run in addition
   # to being built.
   if ( $PROG == "preprocess" ) then
      ./preprocess || exit $n
      set preprocess_done = 1
   endif

skip:
end

if ( $cdebug ) then 
   echo 'preserving .o and .mod files for debugging unless compiling with MPI'
else
   \rm -f *.o *.mod Makefile .cppdefs
endif

\rm -f input.nml*_default preprocess

echo "Success: All single task DART programs compiled."  

if ( $with_mpi ) then
  echo "Script now compiling MPI parallel versions of the DART programs."
else
  echo "Script is exiting after building the serial versions of the DART programs."
  exit 0
endif

#-------------------------------------------------------------------------------
# Build the MPI-enabled target(s) 
#-------------------------------------------------------------------------------

\rm -f *.o *.mod Makefile .cppdefs

foreach PROG ( $MPI_TARGETS )

   set TARGET = `echo $PROG | sed -e 's/^/mkmf_/'`

   @ n = $n + 1
   echo
   echo "---------------------------------------------------"
   echo "${MODEL} (with MPI) build number ${n} is ${PROG}" 

   \rm -f ${PROG} using_mpi_for_${PROG}

   csh $TARGET -mpi || exit $n
   make             || exit $n

   if ( $tdebug ) then
      echo 'removing all files between builds'
      \rm -f *.o *.mod Makefile .cppdefs
   endif

end

if ( $cdebug ) then 
   echo 'preserving .o and .mod files for debugging'
else
   \rm -f *.o *.mod Makefile .cppdefs
endif
\rm -f input.nml*_default

echo "Success: All MPI parallel DART programs compiled."  

exit 0

