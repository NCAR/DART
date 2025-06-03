#!/bin/sh

dart_dir=`pwd`/..

# Make sure we're reading the correct restart files
sed -i '' "s#.*restart_directory.*#   restart_directory = '../restart_files/'#g" $dart_dir/work/input.nml

cd $dart_dir/work/

# Get ensemble size
num_blocks=6  # can this change? 

restarts=$(find $dart_dir/restart_files/ -type f | wc -l)
ens_size=$(($restarts / $num_blocks - 1)) 

# Run AETHER to DART
for iens in $(seq -w 00 $ens_size); do
    $dart_dir/work/aether_to_dart 00${iens} 
done

# Copy one of the members as a template
cp filter_input_0001.nc aether_restart.nc 
