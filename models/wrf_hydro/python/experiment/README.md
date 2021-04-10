# Setup experiments

Contacts:
James (jamesmcc-at-ucar-dot-edu)

## Dependencies

### Python 3.6+

#### Cheyenne
Please configure your python virtual environment using the following script (adapted from https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/software/python):

```
# Cheynne python virtual environments
# James McCreight

# This script installs and activates a new python virtual environment (venv)
# from the latest, stock cisl one.

# -----------------------------------------------------------------------------
# CONFIGURATION SECTION

# First determined the version to be installed
# Run:  module purge; module load ncarenv python/3.x.x; ncar_pylib -l
# Then copy the path to the desired version here.... (unfortunately, this has
# become unparseable and has changed a few times)
py_lib_path=/glade/u/apps/ch/opt/python/3.7.5/gnu/8.3.0/pkg-library/20200417
# The name or tag for the environment that you copy.
# The name will appear in () in your shell (PS1) prompt and the name
# will be used within the directory you use to keep all your different
# virtual environments. 
env_tag=hdp

# This is the dir where you will save different, named/tagged venvs
python_venv_dir=/glade/work/$USER/python_envs

# END CONFIGURATION SECTION
# -----------------------------------------------------------------------------


clone_dir=$python_venv_dir/$env_tag

# -------------------------------------------------------
# Make sure this script is being sourced
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "This script is to be sourced (not run in a sub shell)."
    exit 1
fi

# -------------------------------------------------------
# Install and config
deactivate 2&>1 > /dev/null
rm -rf $clone_dir
ncar_pylib -c $py_lib_path $clone_dir
# The above prints the requisite modules for that env

if [ ! -z $env_tag ]; then
    sed -i "s/(NPL)/(${env_tag})/" ${clone_dir}/bin/activate
fi
source ${clone_dir}/bin/activate || return 1

# -------------------------------------------------------
# Test
pip uninstall -y seaborn || return 1
pip install seaborn || return 1

# -------------------------------------------------------
# User actions
echo "You may want to note the modules listed when this virtual"
echo "env was installed and you may want to do something like:"
echo "module purge; module load <<list of modules printed above>>; module save $env_tag"
echo
echo "To activate this environment on login, put this at the"
echo "bottom of your bashrc:"
echo "source ${clone_dir}/bin/activate"
echo

return 0
```

### wrf\_hydro\_nwm\_public
[https://github.com/NCAR/wrf\_hydro\_nwm\_public](https://github.com/NCAR/wrf\_hydro\_nwm\_public)
The WRF-Hydro model code. The aim is to keep wrf\_hydro\_dart synched with master on this repo.

### wrf_hydro_py
This should be handled by the install of `hydrodartpy` below. Come here if not or if you need a development install of `wrfhydropy`.

[https://github.com/NCAR/wrf\_hydro\_py](https://github.com/NCAR/wrf\_hydro\_py) Python API and tools for managing and running WRF-Hydro and its inputs and outputs, including ensembles. 
```
cd wrf_hydro_py
python setup.py develop
```
As above, this only is run once and updates are obtained by updating the repository.


### Install and test hydrodartpy
```
cd /your/path/to/DART/models/wrf_hydro/hydro_dart_py
python setup.py install

# This is not ideal currently...
cd hydrodartpy/tests/data 
ln -s /your/path/to/wrf_hydro_nwm_public

cd ../local/cheyenne
./suite_1.sh
```
or use `develop` in lieu of `install` if you are editing the package source.


###  WRF-Hydro Domain files in the proper convention 
Currently these can be obtained via
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

