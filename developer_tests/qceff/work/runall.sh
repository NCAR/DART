#!/bin/bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# Emits TAP (one ok / not ok per case) and exits non-zero if any case fails.
# Usage:
#  ./runall.sh
#  ./runall.sh  | grep '^not ok'
#  ./runall.sh  | grep '^ok'
# run_tests.sh runs this script (see developer_tests/tests.list).

source "${DART:-$(git rev-parse --show-toplevel)}/developer_tests/tap.sh"

function set_table () {
echo "&algorithm_info_nml
qceff_table_filename = '$1'
/
$(cat input.nml)" > input.nml
}

function run_test () {
cp input.nml input.nml.orig
set_table $1
# keep program output out of the TAP stream; show it only for failing cases
last_out=$(./test_table_read $1 2>&1)
res=$?
mv input.nml.orig input.nml
return $res
}

should_pass () {
if [[ $? -ne 0 ]]; then
  tap_not_ok "$1" #; echo "$last_out" | tail -n 20 | sed 's/^/# /'
else
  tap_ok "$1"
fi
}

should_fail () {
if [[ $? -eq 0 ]]; then
  tap_not_ok "$1" #; echo "$last_out" | tail -n 20 | sed 's/^/# /'
else
  tap_ok "$1"
fi
}

cp input.nml input.nml.bak
cp input.nml.no_alg input.nml

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

mv input.nml.bak input.nml

tap_done
exit $?
