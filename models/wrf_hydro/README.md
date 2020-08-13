
## models/wrf_hydro: Contents and organization of this directory. 

Contacts: 
James McCreight(jamesmcc-at-ucar-dot-edu)  
Moha Gharamti (mgharamti-at-ucar-dot-edu).

|                                      |                                 |
| ------------------------------------ | ------------------------------- |
| `README.md`                          | This file. |
| R/                                   | R scripts for various purposes. |
| ensemble_config_files/               | Files which configure ensembles in wrfhydropy. |
| experiment_config_files/             | File which configure hydro_dart_py experiments. |
| hydro_dart_py/                       | Python package/library for configuring and executing experiments. |
| python/                              | Python scripts for various purposes. |
| shell_scripts/                       | Shell scripts for various purposes. |
| templates/                           | Obsolete? |
| work/                                | Dart executables build directory and other testing. |
| `model_mod.html`                     | The model_mod documentation. |
| `model_mod.nml`                      | model_mod namelist defaults (subsumed by work/input.nml) |
| `model_mod.f90`                      | The model_mod code. |
| `noah_hydro_mod.f90`                 | Some model_mod interfaces more specific to Noah? |
| `create_identity_streamflow_obs.f90` | For creating identity streamflow obs for the NHDPlus-based channel-network configuration of WRF-Hydro. |

The modules used to test and run HydroDART on cheyenne are:

- ncarenv/1.3   
- ncarcompilers/0.5.0   
- intel/18.0.5
- openmpi/3.1.4
- netcdf/4.6.3   
- ncview/2.1.7   
- nco/4.7.9   

## Setup experiments
See python/experiment.
