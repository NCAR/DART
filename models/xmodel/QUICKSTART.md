# Multi-Model DART Quick Start Guide

## What is This?

The xmodel directory provides infrastructure for compiling multiple DART model_mod files into a single executable. This enables multi-model data assimilation scenarios.

## Quick Start (5 steps)

### 1. Navigate to the xmodel work directory

```bash
cd $DART/models/xmodel/work
```

### 2. Review the configuration

Open `model_config.sh` and check:

**Which models to include:**
```bash
# Default configuration
MODELS=(cam-fv wrf)
```

**How quantities are routed:**
```bash
# Default model for unmapped observation quantities
DEFAULT_INTERPOLATE_MODEL="camfv"

# Which quantities each model handles
MODEL_QTYS["camfv"]="QTY_TEMPERATURE QTY_U_WIND_COMPONENT ..."
MODEL_QTYS["wrf"]="QTY_VERTICAL_VELOCITY ..."
```

You can view the full configuration with:
```bash
./show_config.sh
```

### 3. Verify the setup

Run the verification script to test your configuration:

```bash
./verify_setup.sh
```

This will:
- Validate your configuration
- Test preprocessing of model_mod files
- Show you what will be built
- Help catch configuration errors before building

### 4. Build the executables

If verification passed:

```bash
./quickbuild.sh
```

This will:
- Preprocess each model's model_mod.f90 with model-specific prefixes
- Generate preprocessor flags
- Compile all sources
- Build executables (filter, perfect_model_obs, model_mod_check)

### 5. Check the results

After building, you should have:

```bash
ls -1 filter model_mod_check perfect_model_obs
filter
model_mod_check  
perfect_model_obs
```

## Changing Which Models to Include

### Edit model_config.sh

Edit the `MODELS` array:

```bash
# Two models
MODELS=(cam-fv wrf)

# Three models
MODELS=(cam-fv wrf ROMS_rutgers)

# Different models
MODELS=(lorenz_96 lorenz_63)
```

## Configuring Observation Routing

The system routes observations to models based on observation quantity (QTY). Configure this in `model_config.sh`:

```bash
# Default model for any unmapped quantities
DEFAULT_INTERPOLATE_MODEL="camfv"

# Assign specific quantities to each model
declare -A MODEL_QTYS

# Atmospheric model handles atmospheric quantities
MODEL_QTYS["camfv"]="QTY_TEMPERATURE QTY_U_WIND_COMPONENT QTY_V_WIND_COMPONENT QTY_SURFACE_PRESSURE"

# WRF handles hydrometeor quantities
MODEL_QTYS["wrf"]="QTY_RAINWATER_MIXING_RATIO QTY_GRAUPEL_MIXING_RATIO"
```

**Key points:**
- QTY names must match those in `obs_kind_mod.f90`
- Quantities not explicitly mapped use the default model
- The system prints the routing table during initialization

## Common Use Cases

### Case 1: Atmospheric Models Only

```bash
# In model_config.sh
MODELS=(cam-fv wrf)
LOCATION=threed_sphere
```

### Case 2: 1D Toy Models

```bash
# In model_config.sh
MODELS=(lorenz_96 lorenz_63)
LOCATION=oned
```

### Case 3: Coupled System

```bash
# In model_config.sh
MODELS=(cam-fv MOM6)  # Atmosphere + Ocean
LOCATION=threed_sphere
```

## Troubleshooting

### Error: "No short name mapping for model"

**Problem**: You specified a model in `MODELS` that doesn't have a short name defined.

**Solution**: Add it to `MODEL_SHORT_NAMES` in model_config.sh:

```bash
MODEL_SHORT_NAMES["mymodel"]="mymod"
```

### Error: "model_mod.f90 not found"

**Problem**: The model directory doesn't exist or doesn't have a model_mod.f90.

**Solution**: 
1. Check the model name matches the directory: `ls $DART/models/`
2. Verify model_mod.f90 exists: `ls $DART/models/mymodel/model_mod.f90`

### Error: Location module incompatibility

**Problem**: Models use different location modules.

**Solution**: All models in a multi-model build must use the same location module. Check each model's requirements and choose compatible models.

### Build fails with compilation errors

**Problem**: Could be various issues.

**Solutions**:
1. Run `./verify_setup.sh` first to catch configuration issues
2. Check `preprocessed/` directory has renamed model_mod files
3. Review error messages for missing files or undefined symbols

## Understanding the Build Process

### Step 1: Preprocessing

For each model in `MODELS`, the build script:

1. Takes the original `model_mod.f90` from `$DART/models/MODEL/`
2. Renames the module: `model_mod` → `{short_name}_model_mod`
3. Renames all public functions: `function` → `{short_name}_function`
4. Saves the result to `work/preprocessed/{short_name}_model_mod.f90`

Example:
```fortran
! Original cam-fv/model_mod.f90
module model_mod
  subroutine static_init_model()
  
! Preprocessed camfv_model_mod.f90  
module camfv_model_mod
  subroutine camfv_static_init_model()
```

### Step 3: Compilation

The multi-model `assim_model_mod.f90`:
- Tracks state vector indices for each model
- Routes function calls to the appropriate model_mod

### Step 4: Linking

All sources are compiled and linked together into final executables.

## Next Steps

### Testing Your Build

1. **Run model_mod_check**:
   ```bash
   ./model_mod_check
   ```
   Should initialize all models and report their sizes.

2. **Check preprocess output**:
   ```bash
   ./preprocess
   ```
   Should list observation types from all models.

3. **Examine state vector**:
   Check that the combined state includes all models.

### Implementing Full Functionality

The current implementation has some stub functions that need to be completed based on your specific multi-model requirements:

- **Model advance** (`adv_1step`)
- **Forward operator** (`interpolate`)
- **Perturbations** (`pert_model_copies`)
- **Localization** (`get_close_obs`, `get_close_state`)
- **Vertical conversion** (`convert_vertical_obs`, `convert_vertical_state`)

See README.md section "Current Limitations" for details.

## Examples

### Example 1: Build with default models (cam-fv + wrf)

```bash
cd $DART/models/xmodel/work
./verify_setup.sh    # Check configuration
./quickbuild.sh      # Build everything
```

### Example 2: Build with Lorenz models

```bash
cd $DART/models/xmodel/work

# Edit model_config.sh:
#   MODELS=(lorenz_96 lorenz_63)
#   LOCATION=oned

./verify_setup.sh
./quickbuild.sh
```

### Example 3: Add a new model

```bash
# 1. Edit model_config.sh to add the model:
#    MODELS=(cam-fv wrf my-new-model)
#    MODEL_SHORT_NAMES["my-new-model"]="mynew"

# 2. Edit assim_model_mod.f90 to add support (see README.md)

# 3. Verify and build:
./verify_setup.sh
./quickbuild.sh
```

## Getting Help

- **Check README.md** for detailed documentation
- **Run verify_setup.sh** to diagnose configuration issues  
- **Examine preprocessed/** directory to see renamed files
- **Review build output** for specific error messages

## File Reference

```
models/xmodel/
├── README.md                    # Full documentation
├── QUICKSTART.md               # This file
├── assim_model_mod.f90         # Multi-model wrapper
└── work/
    ├── quickbuild.sh           # Main build script
    ├── model_config.sh         # Configuration (edit this!)
    ├── preprocess_model_mod.sh # Preprocessing script
    ├── verify_setup.sh         # Verification script
    └── preprocessed/           # Generated files (created during build)
```

## Tips

1. **Always run verify_setup.sh first** - catches most issues before building
2. **Start simple** - test with 1 or 2 models before adding more
3. **Check compatibility** - ensure models use the same location module
4. **Clean between builds** - `./quickbuild.sh clean` if you change configuration
5. **Examine preprocessed files** - helps understand what the scripts are doing

## Philosophy

The multi-model system is designed to be:
- **Flexible**: Not limited to 1 or 2 models - include as many as needed
- **Declarative**: Specify models in configuration, not hard-coded
- **Modular**: Each model_mod is preprocessed independently
- **Extensible**: Easy to add new models without changing existing code

