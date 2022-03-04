#!/usr/bin/env bash

main() {

[ -z "$DART" ] && echo "ERROR: Must set DART environment variable" && exit 9
source "$DART"/build_templates/buildconvfunctions.sh

CONVERTER=var
LOCATION=threed_sphere
EXTRA="$DART/models/wrf/model_mod.f90 \
       $DART/models/wrf/module_map_utils.f90"


programs=( \
gts_to_dart \
littler_tf_dart \
rad_3dvar_to_dart \
obs_sequence_tool \
advance_time
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
