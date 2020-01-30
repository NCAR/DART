# models/wrf_hydro: Contents and organization of this directory. 

Contacts: 
James McCreight(jamesmcc-at-ucar-dot-edu)
Tim Hoar (thoar-at-ucar-dot-edu).

```
(dart) jamesmcc@cheyenne6[1383]:~/wrf_hydro_dart/models/wrf_hydro> tree -L 1
.
├── ensemble_config_files/
│      # Files which configure ensembles in wrfhydropy.
├── experiment_config_files/
│      # File which configure hydro_dart_py experiments.
├── hydro_dart_py/
│      # Python package/library for configuring and executing experiments.
├── python/
│      # Python scripts for various purposes.
├── R/
│      # R scripts for various purposes.
├── shell_scripts/
│      # Shell scripts for various purposes.
├── templates/
│      # Obsolete?
├── work/
│      # Dart executables build directory and other testing.
├── model_mod.html
│      # The model_mod documentation.
├── model_mod.nml
│      # The model_mod namelist (subsumed by work/input.nml)
├── model_mod.f90
│      # The model_mod code.
├── noah_hydro_mod.f90
│      # Some model_mod interfaces more specific to Noah?
├── create_identity_streamflow_obs.f90
│      # For creating identity streamflow obs for the NHDPlus-based
│      # channel-network configuration of WRF-Hydro.
├── README.md
       # This file.
```

## Setup experiments
See python/experiment.
