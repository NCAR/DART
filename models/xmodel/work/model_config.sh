#!/usr/bin/env bash

# DART Multi-Model Configuration Template
# Copy this file and customize for your specific multi-model setup

# =============================================================================
# MODEL SELECTION
# =============================================================================
# List the models you want to include in this build
# Model names must match directory names in $DART/models/
# Examples:
#   MODELS=(cam-fv wrf)                    # Two atmospheric models
#   MODELS=(lorenz_96 lorenz_63)          # Two simple models
#   MODELS=(cam-fv wrf ROMS_rutgers)       # Three models
# =============================================================================

MODELS=(cam-fv wrf)
#MODELS=(lorenz_96 lorenz_63)

# =============================================================================
# MODEL SHORT NAMES
# =============================================================================
# Define short names for each model (used for variable prefixes)
# Format: MODEL_SHORT_NAMES["directory-name"]="shortname"
# 
# Short names should be:
#   - Lowercase
#   - No hyphens or special characters
#   - Unique across all models
#
# These will be used for:
#   - Preprocessor flags: -DUSE_SHORTNAME (uppercase)
#   - Module names: shortname_model_mod
#   - Function names: shortname_function_name
# =============================================================================

declare -A MODEL_SHORT_NAMES

# Atmospheric models
MODEL_SHORT_NAMES["cam-fv"]="camfv"
MODEL_SHORT_NAMES["cam-se"]="camse"
MODEL_SHORT_NAMES["wrf"]="wrf"
MODEL_SHORT_NAMES["mpas_atm"]="mpasatm"

# Ocean models
MODEL_SHORT_NAMES["MOM6"]="mom6"
MODEL_SHORT_NAMES["POP"]="pop"
MODEL_SHORT_NAMES["ROMS_rutgers"]="romsrutgers"
MODEL_SHORT_NAMES["mpas_ocn"]="mpasocn"

# Land models
MODEL_SHORT_NAMES["clm"]="clm"
MODEL_SHORT_NAMES["noah"]="noah"

# Simple/test models
MODEL_SHORT_NAMES["lorenz_63"]="l63"
MODEL_SHORT_NAMES["lorenz_96"]="l96"
MODEL_SHORT_NAMES["lorenz_04"]="l04"
MODEL_SHORT_NAMES["lorenz_96_2scale"]="l96_2s"

# Add your custom models here:
# MODEL_SHORT_NAMES["my-model"]="mymod"

# =============================================================================
# LOCATION MODULE
# =============================================================================
# All models in a multi-model build must use the same location module
# Common options:
#   - oned               : 1D locations
#   - threed_sphere      : 3D locations on a sphere (lat/lon/vertical)
#   - threed_cartesian   : 3D Cartesian coordinates
#   - channel            : Channel geometry
# =============================================================================

LOCATION=threed_sphere

# =============================================================================
# OBSERVATION QUANTITY ROUTING
# =============================================================================
# Configure which model should handle which observation quantities (QTYs)
# for the interpolate() routine.
#
# DEFAULT_INTERPOLATE_MODEL: Model to use when QTY is not explicitly mapped
# MODEL_QTYS: Space-separated list of QTY names each model should handle
#
# QTY names come from obs_kind_mod (e.g., QTY_TEMPERATURE, QTY_U_WIND_COMPONENT)
# Use the exact names as defined in assimilation_code/modules/observations/obs_kind_mod.f90
# =============================================================================

# Default model for unmapped quantities
DEFAULT_INTERPOLATE_MODEL="camfv"

declare -A MODEL_QTYS

# Example: CAM-FV handles atmospheric quantities
MODEL_QTYS["camfv"]="QTY_TEMPERATURE QTY_U_WIND_COMPONENT QTY_V_WIND_COMPONENT QTY_SURFACE_PRESSURE QTY_SPECIFIC_HUMIDITY QTY_PRESSURE"

# Example: WRF handles additional atmospheric quantities
MODEL_QTYS["wrf"]="QTY_VERTICAL_VELOCITY QTY_RAINWATER_MIXING_RATIO QTY_GRAUPEL_MIXING_RATIO"

# Add your model's QTYs here:
# MODEL_QTYS["my-model"]="QTY_SALINITY QTY_SEA_SURFACE_HEIGHT"

# =============================================================================
# MODEL-SPECIFIC DEPENDENCIES
# =============================================================================
# Some models require additional source files from other directories
# Add them here as needed
#
# Example: CAM-FV needs cam-common-code
# =============================================================================

declare -A MODEL_EXTRAS

# CAM models need common code
MODEL_EXTRAS["cam-fv"]="$DART/models/cam-common-code"
MODEL_EXTRAS["cam-se"]="$DART/models/cam-common-code"

# Add dependencies for your models here:
# MODEL_EXTRAS["my-model"]="$DART/models/my-model-support $DART/other/path"

# =============================================================================
# BUILD CONFIGURATION
# =============================================================================

# Directories to exclude from automatic source file collection
# (Space-separated list for each model if needed)
declare -A MODEL_EXCLUDES

# Example: WRF excludes some experimental directories
MODEL_EXCLUDES["wrf"]="experiments"

# Add exclusions for your models:
# MODEL_EXCLUDES["my-model"]="tests documentation"

# =============================================================================
# LAG CONFIGURATION
# =============================================================================
# To add lagged copies of a model to the state, set LAG_MODELS and NLAGS.
# For each model in LAG_MODELS, NLAGS lagged entries are added automatically:
#   - Symlinks are created: $DART/models/<model>_lag1 -> $DART/models/<model>
#   - MODELS array is extended with the lag directory names
#   - MODEL_SHORT_NAMES, MODEL_EXTRAS, and MODEL_EXCLUDES are propagated
#     from the base model to each lag instance
#
# Example (add 2 lagged WRF instances alongside the primary WRF):
#   MODELS=(cam-fv wrf)
#   LAG_MODELS=(wrf)
#   NLAGS=2
# This produces models: cam-fv  wrf  wrf_lag1  wrf_lag2
# with short names:     camfv   wrf  wrflag1   wrflag2
# =============================================================================

LAG_MODELS=()   # Models to add lagged copies of, e.g. LAG_MODELS=(wrf)
NLAGS=0         # Number of lags per model (0 = disabled)

# Internal: expand lag entries into MODELS, MODEL_SHORT_NAMES, EXTRAS, EXCLUDES
function setup_lags() {
  if [ ${#LAG_MODELS[@]} -eq 0 ] || [ "$NLAGS" -le 0 ]; then
    return
  fi

  for base_model in "${LAG_MODELS[@]}"; do
    local base_short="${MODEL_SHORT_NAMES[$base_model]}"
    if [ -z "$base_short" ]; then
      echo "ERROR: No short name defined for lag base model '$base_model'"
      exit 1
    fi

    for i in $(seq 1 "$NLAGS"); do
      local lag_dir="${base_model}_lag${i}"
      local lag_short="${base_short}lag${i}"

      # Create symlink in $DART/models/ pointing to the base model directory
      if [ ! -e "$DART/models/$lag_dir" ]; then
        echo "Creating symlink: $DART/models/$lag_dir -> $DART/models/$base_model"
        ln -s "$DART/models/$base_model" "$DART/models/$lag_dir"
      fi

      # Extend the MODELS array
      MODELS+=("$lag_dir")

      # Register short name
      MODEL_SHORT_NAMES["$lag_dir"]="$lag_short"

      # Propagate EXTRAS from base model to this lag
      if [ -n "${MODEL_EXTRAS[$base_model]}" ]; then
        MODEL_EXTRAS["$lag_dir"]="${MODEL_EXTRAS[$base_model]}"
      fi

      # Propagate EXCLUDES from base model to this lag
      if [ -n "${MODEL_EXCLUDES[$base_model]}" ]; then
        MODEL_EXCLUDES["$lag_dir"]="${MODEL_EXCLUDES[$base_model]}"
      fi
    done
  done
}

# =============================================================================
# VALIDATION
# =============================================================================

function validate_config() {
  local errors=0
  
  echo "Validating configuration..."
  
  # Check that we have at least one model
  if [ ${#MODELS[@]} -eq 0 ]; then
    echo "ERROR: No models specified in MODELS array"
    errors=$((errors + 1))
  fi
  
  # Check that each model has a short name
  for model in "${MODELS[@]}"; do
    if [ -z "${MODEL_SHORT_NAMES[$model]}" ]; then
      echo "ERROR: No short name defined for model '$model'"
      echo "       Add it to MODEL_SHORT_NAMES array"
      errors=$((errors + 1))
    fi
    
    # Check that model directory exists
    if [ ! -d "$DART/models/$model" ]; then
      echo "ERROR: Model directory not found: $DART/models/$model"
      errors=$((errors + 1))
    fi
    
    # Check that model_mod.f90 exists
    if [ ! -f "$DART/models/$model/model_mod.f90" ]; then
      echo "ERROR: model_mod.f90 not found for model '$model'"
      echo "       Expected at: $DART/models/$model/model_mod.f90"
      errors=$((errors + 1))
    fi
  done
  
  # Check for duplicate short names
  declare -A seen_short_names
  for model in "${MODELS[@]}"; do
    short="${MODEL_SHORT_NAMES[$model]}"
    if [ -n "${seen_short_names[$short]}" ]; then
      echo "ERROR: Duplicate short name '$short' used for models:"
      echo "       ${seen_short_names[$short]} and $model"
      errors=$((errors + 1))
    fi
    seen_short_names[$short]=$model
  done
  
  if [ $errors -eq 0 ]; then
    echo "Configuration is valid!"
    return 0
  else
    echo "Configuration has $errors error(s)"
    return 1
  fi
}

# =============================================================================
# CONFIGURATION SUMMARY
# =============================================================================

function print_config_summary() {
  echo ""
  echo "========================================="
  echo "Multi-Model Build Configuration Summary"
  echo "========================================="
  echo ""
  echo "Models to include (${#MODELS[@]}):"
  for model in "${MODELS[@]}"; do
    short="${MODEL_SHORT_NAMES[$model]}"
    echo "  - $model (short name: $short)"
    if [ -n "${MODEL_EXTRAS[$model]}" ]; then
      echo "    Extras: ${MODEL_EXTRAS[$model]}"
    fi
    if [ -n "${MODEL_EXCLUDES[$model]}" ]; then
      echo "    Excludes: ${MODEL_EXCLUDES[$model]}"
    fi
  done
  echo ""
  echo "Location module: $LOCATION"
  if [ -n "$CUSTOM_CPPFLAGS" ]; then
    echo "Custom preprocessor flags: $CUSTOM_CPPFLAGS"
  fi
  echo ""
  echo "Lag configuration:"
  if [ ${#LAG_MODELS[@]} -gt 0 ] && [ "$NLAGS" -gt 0 ]; then
    echo "  Lagged models (${#LAG_MODELS[@]}):"
    for model in "${LAG_MODELS[@]}"; do
      echo "    - $model: $NLAGS lag(s)"
    done
  else
    echo "  No lags configured (NLAGS=$NLAGS)"
  fi
  echo ""
  echo "Observation quantity routing:"
  echo "  Default interpolate model: $DEFAULT_INTERPOLATE_MODEL"
  echo ""
  echo "  Model-specific quantities:"
  for model in "${MODELS[@]}"; do
    short="${MODEL_SHORT_NAMES[$model]}"
    qtys="${MODEL_QTYS[$short]}"
    if [ -n "$qtys" ]; then
      echo "    $short:"
      # Print QTYs in a more readable format
      for qty in $qtys; do
        echo "      - $qty"
      done
    fi
  done
  echo ""
  echo "========================================="
  echo ""
}

# Export configuration (called by quickbuild.sh)
function verify_config() {
  setup_lags
  validate_config || exit 1
  print_config_summary
}
