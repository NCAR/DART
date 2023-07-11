#!/bin/bash

# check if two arguments are provided
if [ $# -lt 2 ]; then
  echo "Error, usage: compare_obsseq <obs_seq file 1> <obs_seq file 2>"
  exit 1
fi

file1="$1"
file2="$2"

# check if both files exist
if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
  echo "Error: one or both files do not exist"
  exit 1
fi

# compare the files using diff
diff_output=$(diff "$file1" "$file2")

# check if the outputs are different
if [ -z "$diff_output" ]; then
  echo "Files are Identical"
else
  echo "Files are different"
  echo "$diff_output"
fi
