# Scripts for DART - WRF-Hydro

This document is continually evolving and therefore may be out of date. If something doesnt work, please ask James (jamesmcc-at-ucar-dot-edu) or Tim (thoar-at-ucar-dot-edu).

## Dependencies
* Python 3.6+

## Experiment setup

TODO: add figure of experiment/run conceptualization.

To setup an experiment, change your working directory to your
`wrf_hydro_dart/models/wrf_hydro/shell_scripts` directory. The invocation:

```
   python setup_experiment.py experiment_config_files/your_config.yaml
```

if you want to use debugging, ipython offers the cleanest way:

```
   ipython --pdb setup_experiment.py experiment_config_files/your_config.yaml
```

Here:
* `python` (or `ipython`):
..Version 3.6+
* `setup_experiment.py`:
..The script run/exercised in this example.
* `your_config.yaml`:
..The YAML setup configuration file which contains \*all\* the moving parts of the experiment.
                           
Philosophy: The `your_config.yaml` file should encapsulate all aspects of the experiment setup. If you are finding this is not the case, please identify what parts of the setup are falling outside its scope to James.

An example case will be provide in a later section.

## Run filter experiment

At the end of the setup, instructions will be printed. 
