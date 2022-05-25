#!/usr/bin/env bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
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

# Build the DART programs 
echo 'runing quickbuild.sh'
./quickbuild.sh nompi

echo 'running perfect_model_obs'
./perfect_model_obs || exit 41

echo 'running filter'
./filter            || exit 51

exit 0


