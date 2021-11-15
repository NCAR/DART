#!/bin/bash

main() {
set -e

[ -z "$DART" ] && echo "ERROR: Must set DART environment variable" && exit 9

CONVERTER=AIRS
LOCATION=threed_sphere
source $DART/build_templates/buildconvfunctions.sh

programs=( \
L1_AMSUA_to_netcdf
advance_time
convert_airs_L2
convert_amsu_L1
obs_sequence_tool
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
