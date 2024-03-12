#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

main() {

export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildfunctions.sh

MODEL=aether_lat-lon
LOCATION=threed_sphere

programs=(
filter
model_mod_check
perfect_model_obs
)

serial_programs=(
create_fixed_network_seq
create_obs_sequence
obs_diag
obs_seq_to_netcdf
)

model_serial_programs=(
aether_to_dart
dart_to_aether)

arguments "$@"

# clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# build any NetCDF files from .cdl files
cdl_to_netcdf

# build and run preprocess before making any other DART executables
buildpreprocess

# build DART
buildit

# clean up
\rm -f -- *.o *.mod

}

main "$@"
