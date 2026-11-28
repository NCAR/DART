# Multi-Model DART - Documentation Index

Welcome to the DART multi-model system! This directory contains everything you need to compile multiple model_mod files into a single DART executable.

## Start Here

**New users**: Start with [QUICKSTART.md](QUICKSTART.md) - a 5-step guide to get running quickly.

**Detailed information**: See [README.md](README.md) - comprehensive documentation with examples and troubleshooting.

**Complete example**: See [EXAMPLE_WORKFLOW.md](EXAMPLE_WORKFLOW.md) - step-by-step walkthrough of building with cam-fv and wrf.

**Implementation details**: See [SUMMARY.md](SUMMARY.md) - technical summary of the design and implementation.

## Key Features

✅ **Auto-generated multi-model wrapper** - No manual editing required  
✅ **Automatic domain tracking** - Monitors which domains belong to which model  
✅ **QTY-based observation routing** - Configure which model handles which observations  
✅ **Multi-domain file I/O** - Each model writes its own files correctly  
✅ **Flexible configuration** - Easy to add/remove models  
✅ **Build verification** - Catch configuration errors before building

## Quick Links

### Getting Started
- [Quick Start (5 steps)](QUICKSTART.md#quick-start-5-steps)
- [Changing which models to include](QUICKSTART.md#changing-which-models-to-include)
- [Common use cases](QUICKSTART.md#common-use-cases)

### Configuration
- [Model selection](README.md#configuration-details)
- [Model short names](README.md#model-short-names)
- [Observation quantity routing](README.md#observation-quantity-routing)
- [Location module](README.md#location-module)
- [Extra dependencies](README.md#extra-dependencies)

### Building
- [Basic build](QUICKSTART.md#quick-start-5-steps)
- [Build with different models](EXAMPLE_WORKFLOW.md#changing-configuration)
- [Build options (MPI, single program)](EXAMPLE_WORKFLOW.md#common-scenarios)

### Adding Models
- [How to add a new model](README.md#adding-a-new-model)
- [Model configuration checklist](README.md#adding-a-new-model)
- [Updating assim_model_mod](README.md#adding-a-new-model)

### Troubleshooting
- [Common problems and solutions](QUICKSTART.md#troubleshooting)
- [Troubleshooting workflow](EXAMPLE_WORKFLOW.md#troubleshooting-workflow)
- [Error messages explained](README.md#troubleshooting)

### Implementation
- [Design pattern](SUMMARY.md#design-pattern)
- [How it works](README.md#how-it-works)
- [File structure](README.md#file-structure)
- [Current limitations](README.md#current-limitations)

## Essential Commands

```bash
# Show current configuration
./show_config.sh

# Verify setup before building
./verify_setup.sh

# Build all executables
./quickbuild.sh

# Clean build artifacts
./quickbuild.sh clean

# Build without MPI
./quickbuild.sh nompi

# Build single program
./quickbuild.sh filter
```

## Configuration File

- **model_config.sh** - Main configuration file (edit this to select models)

## Generated Files

During the build, the following are created in `work/`:
- **preprocessed/** - Directory containing renamed model_mod files
- **assim_model_mod.f90** - Copied from parent directory
- **executables** - Built DART programs (filter, perfect_model_obs, etc.)

## Documentation Files

| File | Purpose | When to Read |
|------|---------|--------------|
| [QUICKSTART.md](QUICKSTART.md) | Quick 5-step guide | First time using the system |
| [README.md](README.md) | Comprehensive documentation | Need detailed information |
| [EXAMPLE_WORKFLOW.md](EXAMPLE_WORKFLOW.md) | Complete walkthrough | Want to see full example |
| [SUMMARY.md](SUMMARY.md) | Technical summary | Understanding the design |
| [INDEX.md](INDEX.md) | This file | Finding documentation |

## Utility Scripts

| Script | Purpose | When to Use |
|--------|---------|-------------|
| quickbuild.sh | Build executables | Ready to build |
| verify_setup.sh | Test configuration | Before building |
| show_config.sh | Display config | Check current settings |
| model_config.sh | Configure models | Changing which models to include |
| generate_assim_model_mod.sh | Auto-generate wrapper | Called automatically by quickbuild |
| preprocess_model_mod.sh | Rename model_mod | Called automatically by quickbuild |

## Typical Workflow

1. **Configure**: Edit `model_config.sh` to select models
2. **Verify**: Run `./verify_setup.sh` to test configuration
3. **Build**: Run `./quickbuild.sh` to compile executables
4. **Test**: Run `./model_mod_check` to verify initialization
5. **Use**: Run DART programs (filter, perfect_model_obs, etc.)

## Examples by Use Case

### I want to build with the default models (cam-fv + wrf)
→ See [QUICKSTART.md](QUICKSTART.md#quick-start-5-steps)

### I want to use different models
→ See [QUICKSTART.md](QUICKSTART.md#changing-which-models-to-include)

### I want to add a new model
→ See [README.md](README.md#adding-a-new-model)

### I'm getting build errors
→ See [QUICKSTART.md](QUICKSTART.md#troubleshooting)

### I want to understand how it works
→ See [SUMMARY.md](SUMMARY.md#design-pattern)

### I want a complete step-by-step example
→ See [EXAMPLE_WORKFLOW.md](EXAMPLE_WORKFLOW.md)

## Getting Help

1. **Check the documentation** - Most questions are answered in one of the docs
2. **Run verify_setup.sh** - Catches most configuration issues
3. **Review error messages** - They often point to the specific problem
4. **Check the examples** - See working configurations

## Design Philosophy

The multi-model system is designed to be:
- **Flexible**: Support any number of models, not just 2
- **Declarative**: Specify models in configuration, not hard-coded
- **Modular**: Each model is preprocessed independently
- **Extensible**: Easy to add new models without changing existing code
- **Verifiable**: Test configuration before building

## Contributing

When adding new models or features:
1. Update `model_config.sh` with model mappings
2. Update `assim_model_mod.f90` with model-specific code
3. Test with `verify_setup.sh`
4. Document any special requirements in README.md

