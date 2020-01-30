# ###################################
# HydroDart Experiment Setup Examples
# Overview:

# This is for a small domain, sixmile, using the WRF-Hydro NWM configuration
# and DART's RMA branch.

# Note there is a distinction between the "setup" and "run" phases of the experiment.
# This example is for the setup. 


# The HydroDart Experiment Setup is invoked from the wrf_hydro_dart/models/wrf_hydro/shell_scripts dir:
#
#   python setup_experiment.py experiment_config_files/hydro_dart_config.cheyenne.sixmile.yaml
#
# if you want to use debugging, ipython offers the cleanest way:
#
#   ipython --pdb setup_experiment.py experiment_config_files/hydro_dart_config.cheyenne.sixmile.yaml
#
# where:
#   python (or ipython)   : version 3.6.4+
#   setup_experiment.py   : the script run/exercised in this example.
#   hydro_dart_config.yaml: is the YAML setup configuration file
#                           which contains *all* the moving parts of the experiment.

# THEN CD to RUN_DIR and `chmod 755 run_filter.csh; ./run_filter.csh`

# Basically everything should be specified in the hydro_dart_config.yaml file.

# 1. The domain.
#    Domains with the *correct conventions* to work with wrfhydropy are kept here on
#    cheyenne:
#      /glade/p/work/jamesmcc/domains
#    If you are unsure that you should have access to the private subdir, please ask.
#
# 2. Setup the python 3.6.4 virutal env per cisl instructions.
#    Python 3.6.4+ is required.
#    Both sections: 
#       https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/software/python#modules
#       https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/software/python#library
#    including "Creating your own clone of the NCAR package library".
#
# 3. Install the wrfhydropy prequisites.
#    These are currently summarized here:
#       https://github.com/NCAR/wrf_hydro_docker/blob/160da2458e9be7313636910051fa8887776fb7be/dev/conda/Dockerfile#L40-L43
#    Currently (may be out of date) this is, for example:
#       pip install jupyter cartopy rasterio netcdf4 dask f90nml deepdiff xarray plotnine boltons pytest pytest-datadir-ng wrfhydropy
#    A development version of wrfhydropy is needed, you'll need to clone this repository to cheyenne
#    git clone wrf_hydro_py
#    cd /path/to/wrf_hydro_py/
#    git checkout james
#    pip uninstall -y wrfhydropy
#    python setup.py develop








# The hydro_dart_config.yaml in this example contains keys and values relevant to all  phases
# of the experiment:
#
# Build Phase: dart and wrf_hydro.
# Experiment setup phase: experiment.
# Ensemble construction phase: ensemble.
# Run Phase: run.
#
# The most complicated value is the file required to get to run phase. This is
# specified by the keys ensemble:setup_py. In this file, python commands operating
# on the WrfHydroEnsSetup object are specified to build and verify the desired
# ensemble. The object of this class is wrf_hydro_ens_setup and it has two useful
# methods:
#
#   setter: wrf_hydro_ens_setup.set_diffs_dict()
#   getter: wrf_hydro_ens_setup.diffs_dict()
# 
# Because specifying the commands in this file and verifying them  can take some work,
# when no file is specified, the script (setup_experiment.py) enters interactive mode
# for the ensemble construction phase. Here commands and be used to build the ensemble
# from the "vanilla" starting point provided by the script. When interactive mode is
# exited note that the script exits, it will never continue to the run phase. The commands
# used in the interactive portion of the ensemble construction phase must be put into
# a file which is then added to the YAML configuration file. The run phase can only be
# entered by specifying these in that file. 
#
# Determining the commands in ensemble:setup_py can take some time and you may wish to
# start over. When the build phase has already been completed, DART and WRF-Hydro are
# built and exist on disk and specifying the namelist option 'use_existing' for both
# dart and wrf_hydro means these will not be rebuilt, saving significant time.




