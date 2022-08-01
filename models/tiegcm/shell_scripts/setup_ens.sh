#!/bin/bash

# Set up directories for an ensemble of tiegcm runs

ens_size=5

for (( i=1; i<=$ens_size; i++))
do
  mem="mem"$(printf "%02d" $i)
  rsync -a --delete  mem.setup/ $mem
done
