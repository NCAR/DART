#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

main() {

export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildconvfunctions.sh

CONVERTER=MEXICO_AQS
LOCATION=threed_sphere


programs=(
mexico_aqs_co_ascii_to_obs
mexico_aqs_no2_ascii_to_obs
mexico_aqs_o3_ascii_to_obs
mexico_aqs_pm10_ascii_to_obs
mexico_aqs_pm25_ascii_to_obs
mexico_aqs_so2_ascii_to_obs
obs_sequence_tool
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
