#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# This script builds only perfect_model_obs and filter.  To build the rest
# of the executables, run './quickbuild.csh'.
#
# Executes a known "perfect model" experiment using an existing
# observation sequence file (obs_seq.in) and initial conditions appropriate
# for both 'perfect_model_obs' (perfect_input.nc) and 'filter' (filter_input.nc).
# There are enough initial conditions for a 80 member ensemble in filter.
# Use ens_size = 81 and it WILL bomb. Guaranteed.
# The 'input.nml' file controls all facets of this execution.
#
# 'create_obs_sequence' and 'create_fixed_network_sequence' were used to
# create the observation sequence file 'obs_seq.in' - this defines
# what/where/when we want observations. This script builds these
# programs in support of the tutorial exercises but does not RUN them.
#
# 'perfect_model_obs' results in a true_state.nc file that contains
# the true state, and obs_seq.out - a file that contains the 
# synthetic "observations" that will be assimilated by 'filter'.
#
# 'filter' results in three files (at least): preassim.nc - the state
# of all ensemble members prior to the assimilation (i.e. the forecast),
# analysis.nc - the state of all ensemble members after the
# assimilation, and obs_seq.final - the ensemble members'
# estimate of what the observations should have been.
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

# The input model states for both perfect_model_obs and filter come
# from netCDF files and must be built from the source .cdl files.

which ncgen > /dev/null
if ($status != 0) then
  echo "The required input netCDF files must be build using 'ncgen'"
  echo "'ncgen' is not currently available. It comes with every"
  echo "netCDF installation and is needed by DART. Stopping."
  exit 1
endif

if ( ! -e perfect_input.nc ) ncgen -o perfect_input.nc perfect_input.cdl
if ( ! -e  filter_input.nc ) ncgen -o  filter_input.nc  filter_input.cdl

#----------------------------------------------------------------------
# 'preprocess' is a program that culls the appropriate sections of the
# observation module for the observations types in 'input.nml'; the
# resulting source file is used by all the remaining programs,
# so it MUST be run first.
#----------------------------------------------------------------------

\rm -f preprocess *.o *.mod

set MODEL = "forced_lorenz_96"

echo 'building and running preprocess'

csh  mkmf_preprocess
make         || exit 1
./preprocess || exit 99

echo 'building create_obs_sequence'
csh mkmf_create_obs_sequence
make || exit 2

echo 'building create_fixed_network_seq'
csh mkmf_create_fixed_network_seq
make || exit 3

echo 'building perfect_model_obs'
csh mkmf_perfect_model_obs
make || exit 4

echo 'building filter'
csh mkmf_filter
make || exit 5

echo 'building obs_diag'
csh mkmf_obs_diag
make || exit 6

echo 'removing the compilation leftovers'
\rm -f *.o *.mod

#----------------------------------------------------------------------
# Tutorial section 20 states that perfect_model_obs has forcing
# fixed at 8.0, and can vary after that.
#----------------------------------------------------------------------

echo '/model_nml/'      >! ex_script
echo '/ forcing '       >> ex_script
echo 's;=.*;= 8.0,;'    >> ex_script
echo '/reset_forcing/'  >> ex_script
echo 's;=.*;= .true.,;' >> ex_script
echo 'wq'               >> ex_script

cat ex_script | ex input.nml || exit 40

echo 'running perfect_model_obs'
./perfect_model_obs || exit 41

#----------------------------------------------------------------------
# For forced L96, we want to allow filter to assimilate forcing.
# Use ex to change value of reset_forcing in namelist.
#----------------------------------------------------------------------

echo '/model_nml/'         >! ex_script
echo '/reset_forcing/'     >> ex_script
echo 's;=.*;= .false.,;'   >> ex_script
echo 'wq'                  >> ex_script

cat ex_script | ex input.nml || exit 50

echo 'running filter'
./filter || exit 51

\rm -f ex_script

exit 0


