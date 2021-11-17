#!/bin/bash

main() {
set -e

[ -z "$DART" ] && echo "ERROR: Must set DART environment variable" && exit 9

CONVERTER=GSI2DART
LOCATION=threed_sphere
source $DART/build_templates/buildconvfunctions.sh

# overwrite mpi variables
mpisrc=mpi
m="-w" 

programs=( \
gsi_to_dart
)

# build arguments
arguments "$@"

# clean the directory
\rm -f *.o *.mod Makefile .cppdefs

# build and run preprocess before making any other DART executables
buildpreprocess

# build 
buildconv


# clean up
\rm -f *.o *.mod

}

main "$@"
