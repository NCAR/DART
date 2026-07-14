#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

main() {


export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildfunctions.sh

MODEL="none"
EXTRA="$DART"/models/template/threed_model_mod.f90
LOCATION="threed_sphere"

serial_programs=(
obs_info
create_obs_grid
compare_states
gen_sampling_err_table
compute_error
obs_assim_count
obs_timejitter
obs_loop
obs_data_denial
obs_remove_dups
obs_sort
obs_total_error
obs_keep_a_few
)

# quickbuild arguments
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
