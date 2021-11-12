#!/bin/bash

main() {
set -e

[ -z "$DART" ] && echo "ERROR: Must set DART environment variable" && exit 9

MODEL=cice
LOCATION=threed_sphere
source $DART/build_templates/buildfunctions.sh

# clean the directory
\rm -f *.o *.mod Makefile .cppdefs


programs=( \
advance_time \
closest_member_tool \
create_fixed_network_seq \
create_obs_sequence \
fill_inflation_restart \
filter \
model_mod_check \
obs_common_subset \
obs_diag \
obs_selection \
obs_seq_coverage \
obs_seq_to_netcdf \
obs_seq_verify \
obs_sequence_tool \
perfect_model_obs \
perturb_single_instance \
wakeup_filter \
)


model_serial_programs=(
cice_to_dart
dart_to_cice
)

arguments "$@"

# clean the directory
\rm -f *.o *.mod Makefile .cppdefs

# build and run preprocess before making any other DART executables
buildpreprocess

# build DART
buildit

# clean up
\rm -f *.o *.mod

}

main "$@"
