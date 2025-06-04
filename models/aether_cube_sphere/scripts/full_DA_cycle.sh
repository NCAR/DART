#!/bin/sh

dart_dir=`pwd`/..

# 0. CLEANUP
rm $dart_dir/dart_log.out || true


# 1. AETHER_TO_DART
$dart_dir/scripts/run_aether_to_dart.sh


# 2. FILTER
# 2.1 Modify the input.nml (e.g., obs to assimilate, evaluate, ...)

# 2.2 Configure the obs_seq file (may need to use obs_sequence_tool)

# 2.3 Move the inflation files (may need to adjust inf configuration in the namelist)
cp $dart_dir/work/output_priorinf_mean.nc $dart_dir/work/input_priorinf_mean.nc
cp $dart_dir/work/output_priorinf_sd.nc $dart_dir/work/input_priorinf_sd.nc

# 2.4 run it 
$dart_dir/work/filter


# 3. DART_TO_AETHER
$dart_dir/scripts/run_dart_to_aether.sh


# 4. INCREMENTS
$dart_dir/scripts/diff_prior_post_blocks.sh


# 5. DIAGNOSTICS?! 
