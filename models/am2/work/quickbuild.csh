#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# Script to manage the compilation of all components for this model;
# executes a known "perfect model" experiment using an existing
# observation sequence file (obs_seq.in) and initial conditions appropriate
# for both 'perfect_model_obs' (perfect_ics) and 'filter' (filter_ics).
# There are enough initial conditions for 80 ensemble members in filter.
# Use ens_size = 81 and it WILL bomb. Guaranteed.
# The 'input.nml' file controls all facets of this execution.
#
# 'create_obs_sequence' and 'create_fixed_network_sequence' were used to
# create the observation sequence file 'obs_seq.in' - this defines
# what/where/when we want observations. This script does not run these
# programs - intentionally.
#
# 'perfect_model_obs' results in a true_state.nc file that contains
# the true state, and obs_seq.out - a file that contains the "observations"
# that will be assimilated by 'filter'.
#
# 'filter' results in three files (depending on values in 'stages_to_write'):
# preassim.nc - all ensemble members prior to assimilation
#               i.e. the (potentially inflated) forecast,
# analysis.nc - all ensemble members after the assimilation (i.e. the analysis), and
# obs_seq.final - the ensemble members' estimates of the observations
#
# Once 'perfect_model_obs' has advanced the model and harvested the
# observations for the assimilation experiment, 'filter' may be run
# over and over by simply changing the namelist parameters in input.nml.
#
# The result of each assimilation can be explored in model-space with
# matlab scripts that directly read the netCDF output, or in observation-space.
# 'obs_diag' is a program that will create observation-space diagnostics
# for any result of 'filter' and results in a couple data files that can
# be explored with yet more matlab scripts.

#----------------------------------------------------------------------
# 'preprocess' is a program that culls the appropriate sections of the
# observation module for the observations types in 'input.nml'; the
# resulting source file is used by all the remaining programs,
# so this MUST be run first.
#----------------------------------------------------------------------

\rm -f preprocess *.o *.mod Makefile .cppdefs

set MODEL = "am2"

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

\rm -f *.o *.mod input.nml*_default Makefile .cppdefs

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

\rm -f *.o *.mod input.nml*_default Makefile .cppdefs

echo
echo 'time to run filter here:'
echo ' for lsf run "bsub < runme_filter"'
echo ' for pbs run "qsub runme_filter"'
echo ' for lam-mpi run "lamboot" once, then "runme_filter"'
echo ' for mpich run "mpd" once, then "runme_filter"'

exit 0


