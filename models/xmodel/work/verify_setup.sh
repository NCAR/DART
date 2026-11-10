#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#-------------------------
# Verify multi-model setup
# Tests preprocessing and configuration before full build
#-------------------------

set -e

export DART=$(git rev-parse --show-toplevel)

# Load configuration
XMODEL_DIR="$DART/models/xmodel/work"
source "$XMODEL_DIR/model_config.sh"

echo "========================================="
echo "Multi-Model Setup Verification"
echo "========================================="
echo ""

# Validate configuration
verify_config

# Test directory setup
echo "Testing directory structure..."
WORK_DIR="$DART/models/xmodel/work"
PREPROC_DIR="$WORK_DIR/preprocessed"
mkdir -p "$PREPROC_DIR"

# Test preprocessing for each model
echo ""
echo "Testing model_mod preprocessing..."
echo ""

all_success=1

for model in "${MODELS[@]}"; do
  short_name="${MODEL_SHORT_NAMES[$model]}"
  model_mod_src="$DART/models/$model/model_mod.f90"
  model_mod_dest="$PREPROC_DIR/${short_name}_model_mod.f90"
  
  echo "Processing: $model -> ${short_name}_model_mod.f90"
  
  # Attempt preprocessing
  if bash "$WORK_DIR/preprocess_model_mod.sh" "$short_name" "$model_mod_src" "$model_mod_dest" 2>&1; then
    echo "  ✓ Preprocessing successful"
    
    # Verify the file was created
    if [ -f "$model_mod_dest" ]; then
      echo "  ✓ Output file created"
      
      # Check for renamed module
      if grep -q "^module ${short_name}_model_mod" "$model_mod_dest"; then
        echo "  ✓ Module renamed correctly"
      else
        echo "  ✗ Module rename failed"
        all_success=0
      fi
      
      # Check for renamed namelist
      if grep -q "${short_name}_model_nml" "$model_mod_dest"; then
        echo "  ✓ Namelist renamed correctly"
      else
        echo "  ✗ Namelist rename failed"
        all_success=0
      fi
      
    else
      echo "  ✗ Output file not created"
      all_success=0
    fi
  else
    echo "  ✗ Preprocessing failed"
    all_success=0
  fi
  echo ""
done

# Summary
echo "========================================="
if [ $all_success -eq 1 ]; then
  echo "✓ All tests passed!"
  echo ""
  echo "Your multi-model setup appears to be configured correctly."
  echo "You can now run: ./quickbuild.sh"
else
  echo "✗ Some tests failed"
  echo ""
  echo "Please review the errors above and check your configuration."
fi
echo "========================================="
echo ""


echo ""
echo "Preprocessed files are in: $PREPROC_DIR"
echo "You can examine them to verify the renaming worked correctly."
echo ""

exit $((1 - all_success))
