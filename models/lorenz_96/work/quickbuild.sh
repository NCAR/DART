#!/bin/bash

main() {
set -e

[ -z "$DART" ] && echo "ERROR: Must set DART environment variable" && exit 9

MODEL=lorenz_96
LOCATION=oned
source $DART/build_templates/buildfunctions.sh

programs=( \
closest_member_tool \
filter \
model_mod_check \
perfect_model_obs
)

serial_programs=( \
create_fixed_network_seq \
create_obs_sequence \
fill_inflation_restart \
integrate_model \
obs_common_subset \
obs_diag \
obs_sequence_tool
)

# quickbuild arguments
arguments "$@"

# clean the directory
\rm -f *.o *.mod Makefile .cppdefs

# build and run preprocess before making any other DART executables
buildpreprocess

# build 
buildit

# clean up
\rm -f *.o *.mod

}

main "$@"
