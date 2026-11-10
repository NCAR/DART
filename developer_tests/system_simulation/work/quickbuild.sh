#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

main() {


export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildfunctions.sh

MODEL="none"
EXTRA="$DART"/models/template/threed_model_mod.f90
dev_test=1
LOCATION="threed_sphere"
TEST="../assimilation_code/programs/system_simulation/"

serial_programs=(
sys_sim103
sys_sim3
sys_sim301
sys_sim2
sys_sim201
sys_sim302
sys_sim102
sys_sim104
obs_sampling_err
sys_sim501
sys_sim402
full_error
sys_sim202
sys_sim401
system_simulation
sys_sim502
sys_sim105
sys_sim4
sys_sim102b
sys_sim5
sys_sim104b
sampling_error
sys_sim101a
test_sampling_err_table
correl_error
sys_sim101
sys_sim203
)

# quickbuild arguments
arguments "$@"

# clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# build and run preprocess before making any other DART executables
buildpreprocess

# build DART
buildit

# clean up
\rm -f -- *.o *.mod

}

main "$@"
