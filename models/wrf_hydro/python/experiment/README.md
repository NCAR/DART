# Setup experiments

Contacts:
James (jamesmcc-at-ucar-dot-edu)

## Dependencies

### Python 3.6+

#### Cheyenne
Please configure your python virtual environment using the instructions at https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/software/python. This is known to not work. The `pip uninstall -y seaborn` command is a check and the environment is broken if there's a permission error on the source of the virtual env. Please report this issue to cisl@ucar.edu or to us.


### Python package dependencies
#### pip
The dependencies can be installed via: 
```
pip install -r requirement.txt
```
Where the `requirements.txt` is under version control in this dir. It is possible that it is out of date if it does not provide all the requirements. 
#### Conda: 
.... todo

### wrf\_hydro\_dart
[https://github.com/NCAR/wrf\_hydro\_dart](https://github.com/NCAR/wrf\_hydro\_dart) is this repo.
The current DART code (RMA) is now used on the master branch. Please make a fork of this and use
the standard "pull from" and "pull-request to" workflow with NCAR/upstream. 

### hydro\_dart\_py
Inside wrf\_hydro\_dart lives the `hydro_dart_py` python package. This package is not available
outside the repo and must be installed by:
```
cd wrf_hydro_dart/hydro_dart_py
python setup.py develop
```
The above will pickup changes whenever the underlying `hydro_dart_py` repository is changed and does 
not need to be run again.

### wrf\_hydro\_nwm\_public
[https://github.com/NCAR/wrf\_hydro\_nwm\_public](https://github.com/NCAR/wrf\_hydro\_nwm\_public)
The WRF-Hydro model code. The aim is to keep wrf\_hydro\_dart synched with master on this repo.

### wrf_hydro_py
[https://github.com/NCAR/wrf\_hydro\_py](https://github.com/NCAR/wrf\_hydro\_py)
..Python API and tools for managing and running WRF-Hydro and its inputs and outputs (including
..ensembles. 
```
cd wrf_hydro_py
python setup.py develop
```
As above, this only is run once and updates are obtained by updating the repository.

###  WRF-Hydro Domain files in the proper convention 
Currently these can be obtained via
1. docker images: https://hub.docker.com/r/wrfhydro/domains/tags/
1. google drive tool (`wrf_hydro_nwm_public/tests/local/utils/gdrive_download.py`), though one has to know the identifier string for the domain currently.
1. On cheyenne
1. By personal contact with NCAR folks. 


## Platforms = linux
We generally run on 
* Cheyenne (PBS scheduler)
* Docker (interactively): wrfhydro/dev:conda, see https://hub.docker.com/r/wrfhydro/dev/tags/


## Experiment setup

To setup an experiment, change your working directory to your
`wrf_hydro_dart/models/wrfHydro/shell_scripts` directory. The invocation:

```
   python setup_filter_experiment.py ../../experiment_config_files/your_config.yaml
```

if you want to use debugging, ipython offers the cleanest way:

```
   ipython --pdb setup_filter_experiment.py ../../experiment_config_files/your_config.yaml
```
this lets you enter and navigate the stack when you get errors. Please see https://docs.python.org/3.6/library/pdb.html#debugger-commands for help on understanding using the interactive debugger.

In the issuing command:
* `python` (or `ipython`): Version 3.6+
* `setup_filter_experiment.py`: The script run/exercised in this example.
* `your_config.yaml`: The YAML setup configuration file which contains \*all\* the moving parts of the experiment.
                           
Philosophy: The `your_config.yaml` file should encapsulate all aspects of the experiment setup. If you are finding this is not the case, please identify what parts of the setup are falling outside its scope to James.

Note: Various parts of the experiment are placed in different areas. The most common errors when setting up an experiment are python saying "I will not overwrite existing files". There are options in the YAML in many of the sections to either "use_existing_build" or to "create". If those respective options are False and True, you may be asking the script to overwrite an existing file and it will complain. To understand what parts of the experiment are staged where, the following diagram is useful. 

<table class="image">
<tr><td><img src="https://drive.google.com/uc?export=view&id=1m3HrEArotXmx3H4dHVoMbstzYz_KfTPl"></td></tr>
<tr><caption align="bottom">Example directory and file structure resulting from setup_filter_experiment. </caption></tr>
</table>

This drawing is on google drive and can be edited with permission: https://docs.google.com/drawings/d/1FJXitueqJDQfD8OBhFXTR3bP8gZ0CEo9Ft9n9DFDs8w/edit

Note that at the very beginning of `setup_filter_experiment.py` it prints something like

```
expdir=/glade/work/jamesmcc/wrfhydro_dart/sixmile/experiments/test1
rundir=/glade/scratch/jamesmcc/wrfhydro_dart/sixmile/runs/test1
inidir=/glade/work/jamesmcc/wrfhydro_dart/sixmile/initial_ens/test1
allobsdir=/glade/work/jamesmcc/wrfhydro_dart/sixmile/obs_seqs/test1
```
to help you alias these directories in your shell. That convention is used in the figure above.

## Run filter experiment

At the end of the setup, instructions will be printed. 

## Examples of setup and run:
Currently the `../../experiment_config_files/cheyenne_sixmile.yaml` is the main working example. It requires the sixmile domain. The vaious `use_existing` and `create` options my not be in the correct state, but the overall experiment should run.

