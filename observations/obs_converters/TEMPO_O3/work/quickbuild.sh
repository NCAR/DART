#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

main() {

export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildconvfunctions.sh

CONVERTER=ATMOS_CHEM/TEMPO_O3
LOCATION=threed_sphere


programs=(
tempo_o3_total_col_ascii_to_obs
tempo_o3_trop_col_ascii_to_obs
tempo_o3_profile_ascii_to_obs
tempo_o3_cpsr_ascii_to_obs
tempo_o3_profile_thinner
tempo_o3_cpsr_thinner
obs_sequence_tool
advance_time
)

# build arguments
arguments "$@"

# clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# build and run preprocess before making any other DART executables
buildpreprocess

# build 
buildconv


# clean up
\rm -f -- *.o *.mod

}

main "$@"
