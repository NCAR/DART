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

./test_table_read qcf_table.txt /dev/null ; should_pass "correct v1 table"

./test_table_read qcf_table_v2.txt  /dev/null  ; should_fail "detect wrong version"

./test_table_read qcf_table_extra_columns.txt /dev/null  ; should_pass "extra colums"

./test_table_read qcf_table_bad_qty.txt  /dev/null  ; should_fail "bad qty"

./test_table_read qcf_table_broke.txt  /dev/null  ; should_fail "bad value"

