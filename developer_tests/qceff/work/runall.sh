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

should_fail () {
if [[ $? -eq 0 ]]; then
   echo $1: "FAIL"
else
   echo $1: "PASS"
fi
}

./test_table_read ; should_pass "no table"

./test_table_read qcf_table.txt ; should_pass "correct v1 table"

./test_table_read qcf_table_v2.txt ; should_fail "detect wrong version"

./test_table_read qcf_table_extra_columns.txt ; should_pass "extra colums"

./test_table_read qcf_table_bad_qty.txt ; should_fail "bad qty"

./test_table_read qcf_table_broke.txt ; should_fail "bad value"

./test_table_read qcf_table_no_header.txt ; should_fail "no header"

./test_table_read qcf_table_lower_gt_upper.txt ; should_fail "upper bound less than lower"

./test_table_read qcf_table_lower_bound_only.txt ; should_pass "lower bound only"

./test_table_read qcf_table_no_bounds_with_values.txt ; should_pass "bounds false, values for bounds"
