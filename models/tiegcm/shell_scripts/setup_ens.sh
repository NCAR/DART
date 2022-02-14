#!/bin/bash


ens_size=5

for (( i==1; i<=$ens_size; i++))
do
  mem="mem"$(printf "%02d" $i)
  cp -r mem.setup $mem

done
