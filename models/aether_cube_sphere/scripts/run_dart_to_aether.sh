#!/bin/sh

dart_dir=`pwd`/..

# Make sure we're modifying the correct restart files
sed -i '' "s#.*restart_directory.*#   restart_directory = '../post_restart_files/'#g" $dart_dir/work/input.nml

# Make the post DA directory if it's not there
# Usually, the restarts are overwritten with teh posteriors
# but I am saving the posteriors separately to be able 
# to look at the increments
mkdir -p $dart_dir/post_restart_files

# If the directory is empty (first time running filter)
# copy the priors and then overwrite them with the posteriors
if [ -z "$(find "$dart_dir/post_restart_files/" -mindepth 1 -print -quit)" ]; then
   cp $dart_dir/restart_files/* $dart_dir/post_restart_files/
fi

# Get ensemble size
num_blocks=6  # can this change? 

restarts=$(find $dart_dir/restart_files/ -type f | wc -l)
ens_size=$(($restarts / $num_blocks - 1))

cd $dart_dir/work/

# Run AETHER to DART
for iens in $(seq -w 00 $ens_size); do
    $dart_dir/work/dart_to_aether 00${iens} 
done
