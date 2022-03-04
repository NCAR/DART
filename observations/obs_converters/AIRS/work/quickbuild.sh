#!/usr/bin/env bash

main() {

[ -z "$DART" ] && echo "ERROR: Must set DART environment variable" && exit 9
source "$DART"/build_templates/buildconvfunctions.sh

CONVERTER=AIRS
LOCATION=threed_sphere
EXTRA="$DART/observations/obs_converters/obs_error/ncep_obs_err_mod.f90"


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
\rm -f -- *.o *.mod Makefile .cppdefs

# build and run preprocess before making any other DART executables
buildpreprocess

# build 
buildconv


# clean up
\rm -f -- *.o *.mod

}

main "$@"
