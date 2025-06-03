#!/bin/sh

dart_dir=`pwd`/..

mkdir -p $dart_dir/increments

prior_dir=$dart_dir/restart_files
poste_dir=$dart_dir/post_restart_files
blockfile=aether_restart
incs_dir=$dart_dir/increments

cd $incs_dir

# Get ensemble size
num_blocks=6  # can this change? 
blocks=$(($num_blocks - 1))

restarts=$(find $dart_dir/restart_files/ -type f | wc -l)
ens_size=$(($restarts / $num_blocks - 1))

# Compute the increments
for iens in $(seq -w 00 $ens_size); do
    echo Increment for member: $iens
    for k in $(seq 0 $blocks); do 
        ncdiff -O $prior_dir/${blockfile}_m00${iens}_g000$k.nc \
                  $poste_dir/${blockfile}_m00${iens}_g000$k.nc \
                  $incs_dir/m00${iens}_g000$k.nc
    done
done
