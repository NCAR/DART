#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

main() {

export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildconvfunctions.sh

CONVERTER=GSI2DART
LOCATION=threed_sphere


# overwrite mpi variables
mpisrc=mpi
m="-w" 

programs=(
gsi_to_dart
)

# build arguments
arguments "$@"

# clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# run the c preprocessor on  enkf/kinds.F90
\rm -f ../enkf/mykinds.f90
cpp -P -D_REAL8_ -traditional-cpp ../enkf/kinds.F90 > ../enkf/mykinds.f90

# build and run preprocess before making any other DART executables
buildpreprocess

# build 
buildconv


# clean up
\rm -f -- *.o *.mod

}

main "$@"
