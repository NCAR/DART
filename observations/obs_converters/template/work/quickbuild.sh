#!/usr/bin/env bash
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

main() {

export DART=$(git rev-parse --show-toplevel)
source "$DART/build_templates/buildconvfunctions.sh"

CONVERTER=$(basename "$(dirname "$PWD")")
LOCATION=threed_sphere

# build arguments
arguments "$@"

# Clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# Build and run preprocess before making any other DART executables
buildpreprocess

buildconv
  
# Clean up
\rm -f -- *.o *.mod
}

main "$@"
