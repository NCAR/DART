#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

main() {

export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildfunctions.sh

MODEL=seir
LOCATION=oned


programs=(
closest_member_tool
filter
model_mod_check
perfect_model_obs
)

serial_programs=(
create_fixed_network_seq
create_obs_sequence
fill_inflation_restart
integrate_model
obs_common_subset
obs_diag
obs_sequence_tool
)

model_programs=(
)

model_serial_programs=(
)

# quickbuild arguments
arguments "$@"

# clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# build any NetCDF files from .cdl files
cdl_to_netcdf

# build and run preprocess before making any other DART executables
buildpreprocess

# build 
buildit

# clean up
\rm -f -- *.o *.mod

}

main "$@"
