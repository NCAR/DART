#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

main() {


export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildfunctions.sh

MODEL=lorenz_96
LOCATION=oned

# quickbuild arguments
arguments "$@"

# clean the directory
\rm -f -- libdart.a *.o *.mod Makefile .cppdefs

# build any NetCDF files from .cdl files
cdl_to_netcdf

# build and run preprocess before making any other DART executables
buildpreprocess

# build static library
buildlib  libdart.a

# clean up
\rm -f -- *.o 

}

main "$@"
