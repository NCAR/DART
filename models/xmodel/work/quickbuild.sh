#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#-------------------------
# Multi-model quickbuild.sh
# Builds DART executables with multiple model_mod files
#-------------------------

main() {

export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildfunctions.sh

# Load multi-model configuration
XMODEL_DIR="$DART/models/xmodel/work"
source "$XMODEL_DIR/model_config.sh"

# verify configuration, setup lags and print for the user:
verify_config

# Build EXTRA from model configurations
EXTRA=""
for model in "${MODELS[@]}"; do
  if [ -n "${MODEL_EXTRAS[$model]}" ]; then
    EXTRA="$EXTRA ${MODEL_EXTRAS[$model]}"
  fi
done

# For xmodel, we use the special model designation
MODEL=xmodel

# Programs to build (standard DART programs)
programs=(
filter
model_mod_check
perfect_model_obs
)

serial_programs=(
)

# Model-specific programs can be included here if needed
model_serial_programs=(
)

arguments "$@"

# Set work/preprocessed directory for preprocessed model_mod files
WORK_DIR="$DART/models/xmodel/work"
PREPROC_DIR="$WORK_DIR/preprocessed"

# Clean the work directory, create preprocessed directory
\rm -f -- *.o *.mod Makefile .cppdefs
\rm -rf "$PREPROC_DIR"
mkdir -p "$PREPROC_DIR"

echo "================================================"
echo "Multi-Model DART Build Configuration"
echo "================================================"
echo "Models to include: ${MODELS[@]}"
echo "Location module: $LOCATION"
echo "Work directory: $WORK_DIR"
echo "Preprocessed files: $PREPROC_DIR"
echo "================================================"

# Preprocess each model's model_mod.f90
for model in "${MODELS[@]}"; do
  short_name="${MODEL_SHORT_NAMES[$model]}"
  if [ -z "$short_name" ]; then
    echo "ERROR: No short name mapping for model '$model'"
    echo "Please add it to MODEL_SHORT_NAMES array"
    exit 1
  fi
  
  model_mod_src="$DART/models/$model/model_mod.f90"
  model_mod_dest="$PREPROC_DIR/${short_name}_model_mod.f90"
  
  if [ ! -f "$model_mod_src" ]; then
    echo "ERROR: Cannot find model_mod.f90 for model '$model'"
    echo "Expected at: $model_mod_src"
    exit 1
  fi
  
  echo "Preprocessing $model model_mod.f90 -> ${short_name}_model_mod.f90"
  bash "$WORK_DIR/preprocess_model_mod.sh" "$short_name" "$model_mod_src" "$model_mod_dest"
  
  if [ $? -ne 0 ]; then
    echo "ERROR: Failed to preprocess model_mod.f90 for $model"
    exit 1
  fi
  
  # Add to EXTRA so it gets compiled
  EXTRA="$EXTRA $model_mod_dest"
done

# Generate assim_model_mod.f90 for the selected models
echo ""
echo "Generating assim_model_mod.f90..."
bash "$WORK_DIR/generate_assim_model_mod.sh" "$WORK_DIR/assim_model_mod.f90"

if [ $? -ne 0 ]; then
  echo "ERROR: Failed to generate assim_model_mod.f90"
  exit 1
fi
echo "Generated assim_model_mod.f90 for models: ${MODELS[@]}"
echo ""

# For each model, add only top-level .f90 files to avoid pulling in
# subdirectory programs (which have their own main() functions and cause
# duplicate symbol errors)
for model in "${MODELS[@]}"; do
  EXTRA="$EXTRA $DART/models/$model/*.f90"
done

echo "Building preprocess..."
# Build and run preprocess before making any other DART executables
buildpreprocess

echo "Building DART executables with multiple models..."

# Handle EXCLUDE directories - join all model excludes
EXCLUDE=""
for model in "${MODELS[@]}"; do
  if [ -n "${MODEL_EXCLUDES[$model]}" ]; then
    EXCLUDE="$EXCLUDE ${MODEL_EXCLUDES[$model]}"
  fi
done

# Build DART
buildit

# Clean up
\rm -f -- *.o *.mod

echo ""
echo "================================================"
echo "Multi-model build complete!"
echo "================================================"
echo "Models included: ${MODELS[@]}"
echo ""
echo "To change which models are included, edit "
echo "model_config.sh"
echo "================================================"

}

main "$@"
