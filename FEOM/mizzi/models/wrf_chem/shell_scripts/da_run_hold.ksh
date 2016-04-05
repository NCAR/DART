#!/bin/ksh -x
#
# Script to hold script execution until all jobs 
# with the $1 job name have completed on bluefire
#   
bjobs > job_list
grep $1 job_list > test_list
while [[ -s test_list ]]; do
   sleep 30
   bjobs > job_list
   grep "$1" job_list > test_list
done
rm job_list test_list
#
