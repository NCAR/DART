#!/usr/bin/env bash

main() {

export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildfunctions.sh

MODEL="none"
EXTRA="$DART"/models/template/threed_model_mod.f90
TEST=get_close
dev_test=1
LOCATION="threed_cartesian"

programs=(
)

serial_programs=(
test_get_close_init_zero_dist
test_get_close_init_zero_dists
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
