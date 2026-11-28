#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#-------------------------
# Preprocess model_mod.f90 to add model-specific prefixes
# This allows multiple model_mods to be compiled into one executable
#-------------------------

set -e

if [ $# -lt 3 ]; then
  echo "Usage: $0 <model_name> <model_mod_path> <output_path>"
  echo "  model_name: short name for the model (e.g., camfv, wrf)"
  echo "  model_mod_path: path to the original model_mod.f90"
  echo "  output_path: path for the renamed model_mod.f90"
  exit 1
fi

MODEL_NAME=$1
MODEL_MOD_PATH=$2
OUTPUT_PATH=$3

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$OUTPUT_PATH")"

echo "Preprocessing $MODEL_MOD_PATH for model: $MODEL_NAME"
echo "Output will be written to: $OUTPUT_PATH"

# Start with the original file
cp "$MODEL_MOD_PATH" "$OUTPUT_PATH"

# Replace module name: module model_mod -> module {model}_model_mod
sed -i.bak "s/^module model_mod$/module ${MODEL_NAME}_model_mod/g" "$OUTPUT_PATH"
sed -i.bak "s/^end module model_mod$/end module ${MODEL_NAME}_model_mod/g" "$OUTPUT_PATH"

# Replace namelist name: model_nml -> {model}_model_nml
sed -i.bak "s/model_nml/${MODEL_NAME}_model_nml/g" "$OUTPUT_PATH"

# Clean up backup files
\rm -f "${OUTPUT_PATH}.bak"

echo "Preprocessing ${MODEL_NAME} complete!"
echo "Module renamed to: ${MODEL_NAME}_model_mod"
echo "Namelist renamed to: ${MODEL_NAME}_model_nml"
