# Multi-Model DART Support (xmodel)

## Overview

The `xmodel` directory provides infrastructure for compiling multiple DART model_mod files into a single executable. This allows for multi-model data assimilation where observations can be assimilated into multiple models simultaneously.

## Design

### Key Components

1. **preprocess_model_mod.sh** - Preprocesses model_mod.f90 files to add model-specific prefixes
2. **assim_model_mod.f90** - Multi-model version that wraps multiple model_mods
3. **quickbuild.sh** - Build script that handles multiple models

### How It Works

1. Each model's `model_mod.f90` is preprocessed to rename:
   - The module name: `model_mod` → `{model}_model_mod`
   - All public routines: `routine_name` → `{model}_routine_name`
   
2. The multi-model `assim_model_mod.f90`:
   - **Auto-generated** from models listed in `model_config.sh`
   - Imports all model_mods with renamed interfaces
   - Tracks state vector indices for each model
   - **Automatically tracks domains** added by each model during initialization
   - Routes calls to the appropriate model_mod based on:
     - State vector index (for state-related functions)
     - Domain ID (for file I/O operations)
     - Observation quantity (for forward operator)
   - Combines outputs from multiple models

3. The build system:
   - Preprocesses each model's model_mod.f90
   - Auto-generates assim_model_mod.f90 for selected models
   - Sets preprocessor flags for enabled models
   - Compiles all sources together

4. **Domain Tracking** (automatic):
   - Each model calls `add_domain()` during its `static_init_model()`
   - The wrapper automatically tracks which domains belong to which model
   - File I/O operations (`write_model_time`, `nc_write_model_atts`) route based on domain
   - No manual domain registration required

5. **Observation Routing** (configured):
   - Map observation quantities (QTYs) to models in `model_config.sh`
   - `interpolate()` automatically routes each observation to the correct model
   - Supports default fallback model for unmapped quantities

## Quick Start

### 1. Configure Models

Edit `work/quickbuild.sh` and modify the `MODELS` array:

```bash
# Include cam-fv and wrf
MODELS=(cam-fv wrf)

# Or include different models
MODELS=(lorenz_96 lorenz_63)
```

### 2. Add Model Mappings

If adding a new model, add its short name mapping:

```bash
declare -A MODEL_SHORT_NAMES
MODEL_SHORT_NAMES["cam-fv"]="camfv"
MODEL_SHORT_NAMES["wrf"]="wrf"
MODEL_SHORT_NAMES["ROMS_rutgers"]="romsrutgers"
MODEL_SHORT_NAMES["my-new-model"]="mynew"  # Add new mapping
```

### 3. Build

```bash
cd $DART/models/xmodel/work
./quickbuild.sh
```

## Configuration Details

### Selecting Models

The models to include are specified in the `MODELS` array in `quickbuild.sh`:

```bash
MODELS=(cam-fv wrf)  # Default: cam-fv and wrf
```

You can include any number of models (not limited to 2):

```bash
MODELS=(cam-fv wrf ROMS_rutgers)  # Three models
```

### Model Short Names

Each model needs a short name for variable prefixes. These are defined in the `MODEL_SHORT_NAMES` associative array:

```bash
MODEL_SHORT_NAMES["cam-fv"]="camfv"   # directory name -> short name
MODEL_SHORT_NAMES["wrf"]="wrf"
MODEL_SHORT_NAMES["ROMS_rutgers"]="romsrutgers"
```

The short name should be:
- Lowercase
- No special characters (will be used in variable names)
- Unique across all models

### Location Module

All models must be compatible with the same location module:

```bash
LOCATION=threed_sphere  # Set in model_config.sh
```

### Observation Quantity Routing

Configure which model should handle which observation quantities for the `interpolate()` function:

```bash
# Default model for unmapped quantities
DEFAULT_INTERPOLATE_MODEL="camfv"

# Assign specific quantities to each model
declare -A MODEL_QTYS
MODEL_QTYS["camfv"]="QTY_TEMPERATURE QTY_U_WIND_COMPONENT QTY_V_WIND_COMPONENT QTY_SURFACE_PRESSURE"
MODEL_QTYS["wrf"]="QTY_VERTICAL_VELOCITY QTY_RAINWATER_MIXING_RATIO QTY_GRAUPEL_MIXING_RATIO"
```

QTY names must match those defined in `obs_kind_mod.f90`. During initialization, the system:
- Creates a lookup table mapping each QTY to its responsible model
- Uses `DEFAULT_INTERPOLATE_MODEL` for any unmapped quantities
- Prints the complete QTY routing table for verification

### Extra Dependencies

If a model requires additional source files, add them to MODEL_EXTRAS
in `model_config.sh`, e.g.:

```bash
# CAM models need common code
MODEL_EXTRAS["cam-fv"]="$DART/models/cam-common-code"
MODEL_EXTRAS["cam-se"]="$DART/models/cam-common-code"
```

### Excludes

Similarly if a model has files that should be excluded from the build, you can specify those in MODEL_EXCLUDES.

```bash
# Example: WRF excludes some experimental directories
MODEL_EXCLUDES["wrf"]="experiments"
```

## File Structure

```
models/xmodel/
├── README.md                 # This file
├── QUICKSTART.md            # Quick start guide (5 steps)
├── SUMMARY.md               # Implementation summary
├── assim_model_mod.f90      # Template (overwritten during build)
└── work/
    ├── quickbuild.sh         # Main build script
    ├── model_config.sh       # Model selection configuration
    ├── generate_assim_model_mod.sh  # Auto-generates assim_model_mod.f90
    ├── preprocess_model_mod.sh  # Model_mod preprocessing script
    ├── verify_setup.sh       # Pre-build verification
    ├── show_config.sh        # Display current configuration
    ├── assim_model_mod.f90   # Auto-generated during build
    └── preprocessed/         # Generated during build
        ├── camfv_model_mod.f90  # Preprocessed cam-fv model_mod
        └── wrf_model_mod.f90    # Preprocessed wrf model_mod
```

## Utility Scripts

The `work/` directory contains several helper scripts:

- **quickbuild.sh** - Main build script (builds all executables)
- **model_config.sh** - Configuration file (edit this to select models)
- **generate_assim_model_mod.sh** - Auto-generates assim_model_mod.f90 (called by quickbuild.sh)
- **preprocess_model_mod.sh** - Preprocesses individual model_mod files (called by quickbuild.sh)
- **verify_setup.sh** - Tests configuration before building
- **show_config.sh** - Displays current configuration

## Adding a New Model

To add a new model to the multi-model system:

1. **Verify compatibility**: Ensure the model uses a compatible location module

2. **Edit model_config.sh**:
   ```bash
   # Add to MODELS array
   MODELS=(cam-fv wrf mynewmodel)
   
   # Add short name mapping
   MODEL_SHORT_NAMES["mynewmodel"]="mnm"
   
   # Configure which observation quantities it handles
   MODEL_QTYS["mnm"]="QTY_SALINITY QTY_SEA_SURFACE_HEIGHT"
   
   # Add any extra dependencies if needed
   MODEL_EXTRAS["mynewmodel"]="$DART/models/mynewmodel/support"
   ```

3. **Build**: That's it! The system will auto-generate assim_model_mod.f90
   ```bash
   cd $DART/models/xmodel/work
   ./verify_setup.sh
   ./quickbuild.sh
   ```

The system automatically:
- Preprocesses your model's model_mod.f90
- Generates the multi-model wrapper
- Tracks domains your model creates
- Routes observations based on your QTY configuration

The `assim_model_mod.f90` file is now **auto-generated** based on your configuration, so you don't need to manually edit it when adding new models

## Current Limitations

### Fully Implemented Features

The following features are **fully implemented** in the multi-model system:

1. **State vector management** - Combined state from all models with proper offset tracking
2. **Model initialization** (`static_init_model`) - Initializes all models and tracks their state
3. **Domain tracking** - Automatic registration of domains to models during initialization
4. **File I/O routing** - `write_model_time` and `nc_write_model_atts` route to correct model based on domain
5. **Forward operator** (`interpolate`) - Routes observations to models based on quantity (QTY) type
6. **Metadata access** (`get_state_meta_data`) - Routes based on state vector index

### QTY-Based Observation Routing

The `interpolate()` function now routes observations to models based on observation quantity:

```bash
# In model_config.sh
DEFAULT_INTERPOLATE_MODEL="camfv"  # Default for unmapped QTYs

MODEL_QTYS["camfv"]="QTY_TEMPERATURE QTY_U_WIND_COMPONENT QTY_V_WIND_COMPONENT"
MODEL_QTYS["wrf"]="QTY_VERTICAL_VELOCITY QTY_RAINWATER_MIXING_RATIO"
```

During initialization, the system builds a lookup table mapping each QTY to its responsible model.

### Time
The first model in the `MODELS` array is considered the "primary" model for time management. The `assim_model_mod` uses this model's time for read_model_time.

### Not Yet Implemented

The following features still need implementation:

1. **Model advance** (`adv_1step`) - Advancing multiple models simultaneously
2. **Perturbations** (`pert_model_copies`) - Perturbing multiple model states
3. **Get close** (`get_close_obs`, `get_close_state`) - Spatial localization across models
4. **Vertical conversion** - Converting between vertical coordinates

These are currently stubs that return errors. They need to be implemented based on specific multi-model requirements.

### Design Considerations

When implementing these features, consider:

1. **Which model handles which observations?**
   - Now handled via QTY-based routing in `interpolate()`
   - Based on observation type?
   - Based on user configuration?

2. **How do models interact?**
   - Are they independent?
   - Do they share state variables?
   - Do they operate on the same grid?

3. **Localization across models**
   - Should observations in one model affect another?
   - How to handle different model resolutions?

## Testing

To verify the multi-model build:

1. Check that preprocessing creates renamed modules:
   ```bash
   ls work/preprocessed/
   # Should show: camfv_model_mod.f90  wrf_model_mod.f90  romsrutgers_model_mod.f90
   ```

2. Run a test executable:
   ```bash
   ./model_mod_check
   # Should initialize both models and report their sizes
   ```

### Problem: Location module incompatibility

**Solution**: All models must use the same location module. Check each model's requirements.

