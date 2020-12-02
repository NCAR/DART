#!/bin/bash
domain_tag=florence_cutout
locality=cheyenne

# tests/experiment_test_files/florence_cutout
test_dir=../../
# this is relative to the test dir
test_domain_dir=experiment_test_files/$domain_tag
yaml_file=$test_domain_dir/${locality}.yaml
answer_file=$test_domain_dir/answers.yaml

cd $test_dir

pytest_cmd="pytest  -s  -vv  --pdb  --exp_yaml $yaml_file  --exp_answer_file $answer_file"

$pytest_cmd
# $pytest_cmd  --use_existing_build 

exit $?
