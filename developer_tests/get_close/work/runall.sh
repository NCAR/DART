#!/bin/bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# Usage:
#  ./runall.sh
#  ./runall.sh  | grep FAIL
#  ./runall.sh  | grep PASS

should_pass () {
if [[ $? -ne 0 ]]; then
  echo $1: "FAIL"
else
  echo $1: "PASS"
fi
}

function run_test () {
  ./$1 > out.$1 2>&1
  grep "$2" out.$1
  res=$?
  return $res
}

# for three_sphere
# Use sed to replace the line (compatible with both macOS and Linux)
if [[ "$OSTYPE" == "darwin"* ]]; then
  sed -i '' 's/^LOCATION=".*"/LOCATION="threed_sphere"/' quickbuild.sh
else
  sed -i 's/^LOCATION=".*"/LOCATION="threed_sphere"/' quickbuild.sh
fi

./quickbuild.sh

run_test test_get_close_init_zero_dist "bad maxdist value"; should_pass "3d sphere max dist <=0"
run_test test_get_close_init_zero_dists "bad maxdist_list value"; should_pass "3d sphere max dists <=0"


#  for threed_cartesian
if [[ "$OSTYPE" == "darwin"* ]]; then
  sed -i '' 's/^LOCATION=".*"/LOCATION="threed_cartesian"/' quickbuild.sh
else
  sed -i 's/^LOCATION=".*"/LOCATION="threed_cartesian"/' quickbuild.sh
fi

./quickbuild.sh

run_test test_get_close_init_zero_dist "bad maxdist value"; should_pass "3d cart max dist <=0"
run_test test_get_close_init_zero_dists "bad maxdist_list value"; should_pass "3d cart max dists <=0"

