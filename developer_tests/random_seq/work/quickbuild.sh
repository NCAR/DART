#!/usr/bin/env bash

main() {


[ -z "$DART" ] && echo "ERROR: Must set DART environment variable" && exit 9
source "$DART"/build_templates/buildfunctions.sh

MODEL="template"
LOCATION="threed_sphere"
dev_test=1
TEST="random_seq"

serial_programs=(
test_corr
test_diff
test_exp
test_gamma
test_gaussian
test_hist
test_inv_gamma
test_random
test_reseed
)


# quickbuild arguments
arguments "$@"

# clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# preprocess not needed for these tests

# build 
buildit

# clean up
\rm -f -- *.o *.mod

}

main "$@"
