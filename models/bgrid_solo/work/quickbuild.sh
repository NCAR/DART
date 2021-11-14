#!/bin/bash

main() {
set -e

[ -z "$DART" ] && echo "ERROR: Must set DART environment variable" && exit 9

MODEL=bgrid_solo
LOCATION=threed_sphere
EXCLUDE=fms_src
EXTRA=extra_source
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

model_serial_programs=( \
column_rand \
id_set_def_stdin \
ps_id_stdin \
ps_rand_local
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
