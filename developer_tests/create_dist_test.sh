#!/usr/bin/env bash

# Check if an argument is passed to the script
if [ -z "$1" ]; then
  echo "Error: No distribution name provided."
  echo "Usage: $0 <distribution_name>"
  exit 1
fi

DIST="$1"
distdir="$DIST"_dist

# Step 1: Create the DIST directory
mkdir -p "$distdir"/work

# Step 2: Copy contents from gamma_dist to DIST
cp -r gamma_dist/work/quickbuild.sh $distdir/work/  
cp -r gamma_dist/test_gamma_dist.f90 $distdir/test_${DIST}_dist.f90
cp -r gamma_dist/work/input.nml $distdir/work/


# Step 3: Replace references of gamma with DIST in the copied files
find "$distdir" -type f -exec sed -i '' "s/gamma/$DIST/g" {} +
