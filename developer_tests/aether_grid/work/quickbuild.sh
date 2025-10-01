#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

main() {

export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildfunctions.sh

MODEL="aether_cube_sphere"
LOCATION="threed_sphere"
dev_test=1
TEST="aether_grid"


programs=(
test_aether_grid
)

serial_programs=(
)

model_programs=(
)

model_serial_programs=(
aether_to_dart
dart_to_aether
)

# quickbuild arguments
arguments "$@"

# clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# build and run preprocess before making any other DART executables
buildpreprocess

# build
buildit

# clean up
\rm -f -- *.o *.mod

}

main "$@"
