#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

main() {

export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildfunctions.sh

MODEL=wrf_unified
LOCATION=threed_sphere
EXCLUDE=experiments


programs=(
closest_member_tool
filter
model_mod_check
perfect_model_obs
)

serial_programs=(
#radiance_obs_seq_to_netcdf  # needs rttov
fill_inflation_restart
obs_diag
obs_sequence_tool
)


model_serial_programs=(
)

arguments "$@"

# clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# build and run preprocess before making any other DART executables
buildpreprocess

# build DART
buildit

# clean up
\rm -f -- *.o *.mod

}

main "$@"
