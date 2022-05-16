#/usr/bin/env bash

main() {

export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildfunctions.sh

MODEL=template
LOCATION=threed_sphere
TEST=quad_interpolate
dev_test=1


programs=(
)

serial_programs=(
test_quad_irreg_interp
test_quad_reg_interp
)

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
