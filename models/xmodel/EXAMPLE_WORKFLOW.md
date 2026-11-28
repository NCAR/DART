# Example Workflow: Building Multi-Model DART Executables

This document walks through a complete example of setting up and building a multi-model DART executable with cam-fv and wrf.

## Initial Setup

### Step 1: Navigate to xmodel directory

```bash
cd /Users/hkershaw/DART/Projects/MultiModelComp/DART/models/xmodel/work
```

### Step 2: Check the current configuration

```bash
./show_config.sh
```

Expected output:
```
=========================================
Multi-Model Build Configuration Summary
=========================================

Models to include (2):
  - cam-fv (short name: camfv)
    Extras: /path/to/DART/models/cam-common-code
  - wrf (short name: wrf)
    Excludes: experiments

Location module: threed_sphere

Observation quantity routing:
  Default interpolate model: camfv

  Model-specific quantities:
    camfv:
      - QTY_TEMPERATURE
      - QTY_U_WIND_COMPONENT
      - QTY_V_WIND_COMPONENT
      - QTY_SURFACE_PRESSURE
      - QTY_SPECIFIC_HUMIDITY
      - QTY_PRESSURE
    wrf:
      - QTY_VERTICAL_VELOCITY
      - QTY_RAINWATER_MIXING_RATIO
      - QTY_GRAUPEL_MIXING_RATIO

=========================================
```

## Verification Phase

### Step 3: Verify the setup

```bash
./verify_setup.sh
```

This will:
1. Validate the configuration
2. Test preprocessing of each model_mod.f90
3. Check that renamed modules are created correctly

Expected output (abbreviated):
```
=========================================
Multi-Model Setup Verification
=========================================

Validating configuration...
Configuration is valid!

Testing model_mod preprocessing...

Processing: cam-fv -> camfv_model_mod.f90
  ✓ Preprocessing successful
  ✓ Output file created
  ✓ Module renamed correctly
  ✓ Functions appear to be renamed

Processing: wrf -> wrf_model_mod.f90
  ✓ Preprocessing successful
  ✓ Output file created
  ✓ Module renamed correctly
  ✓ Functions appear to be renamed

=========================================
✓ All tests passed!

Your multi-model setup appears to be configured correctly.
You can now run: ./quickbuild.sh
=========================================
```

### Step 4: Examine preprocessed files (optional)

```bash
ls -l preprocessed_test/
cat preprocessed_test/camfv_model_mod.f90 | head -50
```

This lets you see exactly how the preprocessing renamed the modules and functions.

## Build Phase

### Step 5: Run the build

```bash
./quickbuild.sh
```

The build process will:

1. **Display configuration**
```
================================================
Multi-Model DART Build Configuration
================================================
Models to include: cam-fv wrf
Location module: threed_sphere
Work directory: /path/to/work
Preprocessed files: /path/to/work/preprocessed
================================================
```

2. **Preprocess each model**
```
Preprocessing cam-fv model_mod.f90 -> camfv_model_mod.f90
Preprocessing complete!
Module renamed to: camfv_model_mod
All public routines prefixed with: camfv_

Preprocessing wrf model_mod.f90 -> wrf_model_mod.f90
Preprocessing complete!
Module renamed to: wrf_model_mod
All public routines prefixed with: wrf_
```

3. **Build preprocess**
```
Building preprocess...
```

4. **Compile and link**
```
Building DART executables with multiple models...
[compilation output...]
```

5. **Report completion**
```
================================================
Multi-model build complete!
================================================
Models included: cam-fv wrf

To change which models are included, edit the
MODELS array at the top of this script.
================================================
```

### Step 6: Verify the build succeeded

```bash
ls -lh filter model_mod_check perfect_model_obs
```

Expected output:
```
-rwxr-xr-x  1 user  group  2.1M Feb 25 10:30 filter
-rwxr-xr-x  1 user  group  1.8M Feb 25 10:30 model_mod_check
-rwxr-xr-x  1 user  group  1.9M Feb 25 10:30 perfect_model_obs
```

## Testing Phase

### Step 7: Test model initialization

input.nml will need to be created with appropriate namelist settings for both models. For example:

```fortran
&wrf_model_nml
...
/

&camfv_model_nml
...
/

&model_mod_check_nml
  input_state_files     = 'caminput.nc', 'wrfinput_d01'
  output_state_files    = 'mmc_output1.nc', 'mmc_output2.nc'
/

```

```bash
./model_mod_check
```



## Changing Configuration

### Example: Switch to different models

```bash
# Edit model_config.sh
vim model_config.sh

# Change MODELS array:
MODELS=(cam-fv wrf)

# Also update LOCATION if needed:
LOCATION=threed_sphere

# And the observation quantity routing if needed:
DEFAULT_INTERPOLATE_MODEL="camfv"
MODEL_QTYS["camfv"]=...
MODEL_QTYS["wrf"]=...
```

### Verify new configuration

```bash
./show_config.sh
./verify_setup.sh
```

### Clean and rebuild

```bash
./quickbuild.sh clean
./quickbuild.sh
```

## Troubleshooting Workflow

### If verify_setup.sh fails:

1. **Read the error messages carefully**
2. **Check model names match directories**:
   ```bash
   ls $DART/models/
   ```
3. **Verify short names are defined**:
   ```bash
   grep "MODEL_SHORT_NAMES" model_config.sh
   ```
4. **Check model_mod.f90 exists**:
   ```bash
   ls $DART/models/cam-fv/model_mod.f90
   ls $DART/models/wrf/model_mod.f90
   ```

### If quickbuild.sh fails:

1. **Run verify_setup.sh first**
2. **Clean and retry**:
   ```bash
   ./quickbuild.sh clean
   ./quickbuild.sh
   ```
3. **Check preprocessed files were created**:
   ```bash
   ls preprocessed/
   ```
4. **Review compilation errors** for specific issues

## File Locations After Build

```
work/
├── assim_model_mod.f90           # generated during build ├── preprocessed/                     # Preprocessed model_mod files
│   ├── camfv_model_mod.f90
│   └── wrf_model_mod.f90
├── filter
├── model_mod_check
├── perfect_model_obs
└── ... (and others if added to quickbuild.sh)
```

## Next Steps

After successfully building:

1. **Create namelists** for your multi-model experiment
2. **Prepare initial conditions** for each model
3. **Set up observations** compatible with both models
4. **Run filter** or other DART programs
5. **Analyze results** from multi-model assimilation

## Summary

The complete workflow is:

```bash
cd $DART/models/xmodel/work
./show_config.sh              # Check what's configured
./verify_setup.sh             # Verify before building
./quickbuild.sh               # Build everything
./model_mod_check             # Test the build

# To rebuild with different models:
# 1. Edit model_config.sh
# 2. Repeat the above steps
```

