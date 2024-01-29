#!/bin/bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# Usage:
#  ./runall.sh
#  ./runall.sh  | grep FAIL
#  ./runall.sh  | grep PASS

function set_table () {
echo "&algorithm_info_nml
qceff_table_filename = '$1'
/
$(cat input.nml)" > input.nml
}

function run_test () {
cp input.nml input.nml.orig
set_table $1
./test_table_read $1
res=$?
cp input.nml.orig input.nml
return $res
}

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

run_test ; should_pass "no table"

run_test qcf_table.txt ; should_pass "correct v1 table"

run_test qcf_table_v2.txt ; should_fail "detect wrong version"

run_test qcf_table_extra_columns.txt ; should_pass "extra colums"

run_test qcf_table_bad_qty.txt ; should_fail "bad qty"

run_test qcf_table_broke.txt ; should_fail "bad value"

run_test qcf_table_no_header.txt ; should_fail "no header"

run_test qcf_table_lower_gt_upper.txt ; should_fail "upper bound less than lower"

run_test qcf_table_lower_bound_only.txt ; should_pass "lower bound only"

run_test qcf_table_no_bounds_with_values.txt ; should_pass "bounds false, values for bounds"

run_test qcf_table_incorrect_filter_kind.txt ; should_fail "incorrect filter_kind"

run_test qcf_table_incorrect_distribution.txt ; should_fail "incorrect distribution"

run_test all_bnrhf_qceff_table.csv ; should_pass "lower case QTY"

run_test qcf_table_lower_case_dist.txt; should_pass "lower case dist_type"

