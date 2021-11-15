#!/bin/bash

main() {
set -e

[ -z "$DART" ] && echo "ERROR: Must set DART environment variable" && exit 9

CONVERTER=MADIS
OBS_ERROR=ncep
LOCATION=threed_sphere
source $DART/build_templates/buildconvfunctions.sh

programs=( \
convert_madis_acars \
convert_madis_marine \
convert_madis_mesonet \
convert_madis_metar \
convert_madis_profiler \
convert_madis_rawin \
convert_madis_satwnd \
obs_seq_to_netcdf \
obs_sequence_tool \
advance_time \
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
