#!/usr/bin/env bash

main() {


[ -z "$DART" ] && echo "ERROR: Must set DART environment variable" && exit 9
source "$DART"/build_templates/buildfunctions.sh

MODEL="template"
LOCATION="threed_sphere"
dev_test=1
TEST="io"

programs=(
test_cf_conventions
test_diag_structure
test_read_write_restarts
test_read_write_time
test_state_structure
)


# quickbuild arguments
arguments "$@"

# clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# build any NetCDF files from .cdl files
cdl_to_netcdf

# build and run preprocess before making any other DART executables
buildpreprocess

# build 
buildit

# clean up
\rm -f -- *.o *.mod

}

main "$@"
