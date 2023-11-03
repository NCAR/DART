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
dev_test=1
TEST="utilities"

serial_programs=(
PrecisionCheck
error_handler_test
file_utils_test
find_enclosing_indices_test
find_first_occurrence_test
nml_test
parse_args_test
sort_test
stacktest
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
