## Multi-Model DART Implementation Summary

### Created Files and Directories

```
models/xmodel/
├── README.md                          # Comprehensive documentation
├── QUICKSTART.md                      # Quick start guide
├── assim_model_mod.f90                # Multi-model wrapper implementation, created during build
└── work/
    ├── quickbuild.sh                  # Main build script
    ├── model_config.sh                # Model selection configuration
    ├── preprocess_model_mod.sh        # Script to rename model_mod routines
    └── verify_setup.sh                # Pre-build verification script
```

### Key Features

1. **Declarative Configuration**
   - Models are selected via the `MODELS` array in `model_config.sh`
   - Not limited to 2 models - can include any number
   - Easy to add/remove models without code changes

2. **Automatic Preprocessing**
   - Each model's `model_mod.f90` is automatically renamed
   - Module name: `model_mod` → `{model}_model_mod`
   - All public routines get model-specific prefixes
   - Example: `static_init_model` → `camfv_static_init_model`

3. **Auto-Generated Multi-Model Wrapper**
   - `assim_model_mod.f90` is generated from `generate_assim_model_mod.sh`
   - No manual editing required when adding new models
   - Automatically includes all configured models

4. **State Vector Management**
   - Combined state vector includes all models
   - Tracks offsets and sizes for each model
   - Routes function calls to appropriate model based on state index

5. **Automatic Domain Tracking**
   - Monitors `add_domain()` calls during each model's initialization
   - Tracks which domains belong to which model
   - Routes file I/O operations based on domain ID
   - No manual domain registration needed

6. **Observation Quantity Routing**
   - Configurable QTY-to-model mapping for `interpolate()`
   - Supports default fallback model for unmapped quantities
   - Debug output shows complete routing table

7. **Multi-Domain File I/O**
   - `write_model_time` routes based on domain
   - `nc_write_model_atts` routes based on domain
   - Each model writes its own file structure

### How to Use

**Basic Usage:**
```bash
cd $DART/models/xmodel/work
./verify_setup.sh      # Check configuration
./quickbuild.sh        # Build executables
```

**To Change Models:**

Edit `model_config.sh`:
```bash
# Change this line:
MODELS=(cam-fv wrf)

# To include different models:
MODELS=(lorenz_96 lorenz_63)
# or
MODELS=(cam-fv wrf ROMS_rutgers)  # Three models
```

**To Configure Observation Routing:**

Edit `model_config.sh`:
```bash
# Set default model for unmapped quantities
DEFAULT_INTERPOLATE_MODEL="camfv"

# Assign quantities to models
MODEL_QTYS["camfv"]="QTY_TEMPERATURE QTY_U_WIND_COMPONENT"
MODEL_QTYS["wrf"]="QTY_RAINWATER_MIXING_RATIO"
```

Then rebuild:
```bash
./quickbuild.sh
```

The system automatically generates `assim_model_mod.f90` for whatever models you specify.

**To Add a New Model:**

1. Add to `MODELS` array
2. Add short name mapping to `MODEL_SHORT_NAMES`
3. Configure which observation quantities it handles in `MODEL_QTYS`
4. Add any extra dependencies to `MODEL_EXTRAS` if needed
5. Run `./verify_setup.sh` and `./quickbuild.sh`

Example in `model_config.sh`:
```bash
MODELS=(cam-fv wrf my-ocean-model)
MODEL_SHORT_NAMES["my-ocean-model"]="ocean"
MODEL_QTYS["ocean"]="QTY_SALINITY QTY_SEA_SURFACE_HEIGHT"
MODEL_EXTRAS["my-ocean-model"]="$DART/models/ocean-support"
```

### Design Pattern

The system uses a **preprocess-and-rename** pattern:

1. **Original model_mod.f90** (cam-fv example):
   ```fortran
   module model_mod
     public :: static_init_model
     subroutine static_init_model()
   ```

2. **Preprocessed camfv_model_mod.f90**:
   ```fortran
   module camfv_model_mod
     public :: camfv_static_init_model
     subroutine camfv_static_init_model()
   ```


### Auto-generation of assim_model_mod.f90

**Implemented:**
- ✅ Directory structure
- ✅ Preprocessing system
- ✅ Build infrastructure
- ✅ Configuration system
- ✅ Verification tools
- ✅ Documentation
- ✅ Basic assim_model_mod wrapper
- ✅ State vector management
- ✅ Model initialization (`static_init_model`)
- ✅ Model size tracking (`get_model_size`)
- ✅ State metadata access (`get_state_meta_data`)
- ✅ Automatic domain tracking
- ✅ Domain-based file I/O routing (`write_model_time`, `nc_write_model_atts`)
- ✅ QTY-based observation routing (`interpolate`)
- ✅ Initial conditions (`init_conditions`, `init_time`)
- ✅ Time management (`read_model_time`, `get_model_time_step`) **Note the first model's time is used as the dart time**

**Not Yet Implemented:**
- ❌ Model advance (`adv_1step`)
- ❌ Perturbations (`pert_model_copies`)
- ❌ Spatial localization (`get_close_obs`, `get_close_state`)
- ❌ Vertical coordinate conversion


### Testing

Run the verification script to test the setup:
```bash
cd $DART/models/xmodel/work
./verify_setup.sh
```

This will:
- Validate configuration
- Test preprocessing
- Show what will be built
- Identify configuration errors

### Documentation Files

- **QUICKSTART.md** - Quick start guide (5 steps to get running)
- **README.md** - Comprehensive documentation with examples
- **SUMMARY.md** - This file

