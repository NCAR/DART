# ###################################
# HydroDart Experiment Setup Examples
# Overview:

**** Currently this is stale for docker. ****

# This is for a small domain, croton_NY, using the WRF-Hydro NWM configuration
# and DART's RMA branch.

# Examples are shown below for both Docker and Cheyenne.

# Note there is a distinction between the "setup" and "run" phases of the experiment.
# This example is for the setup. 

# The HydroDart Experiment Setup is invoked:
#
#   python setup_experiment.py hydro_dart_config.yaml
#
# where:
#   python                : version 3.6.4+
#   setup_experiment.py   : the script exercised in this example.
#   hydro_dart_config.yaml: is the YAML setup configuration file
#                           which contains all the moving parts of the experiment.

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




# ###################################
# Docker setup + run
# (Recall that WRF-Hydro code mounted from case-INsensitive dirves
#  fail to compile and a "git reset --hard" is useful at that point.
# Pull the docker images
docker pull wrfhydro/domains:croton_NY
docker pull wrfhydro/dev:conda
# Estblis the docker data container
docker create --name croton wrfhydro/domains:croton_NY
# When done with the container: docker rm -v croton 
# Where are your repos? 
disk_dir='chimayoSpace'
if [ $HOSTNAME = yucatan.local ]; then disk_dir=jamesmcc; fi
host_repos=/Volumes/d1/${disk_dir}/git_repos
# Start the image
docker run -it \
    -e USER=docker \
    -v ${host_repos}/wrf_hydro_dart:/home/docker/wrf_hydro_dart \
    -v ${host_repos}/wrf_hydro_nwm_public:/home/docker/wrf_hydro_nwm_public \
    -v ${host_repos}/wrf_hydro_py:/home/docker/wrf_hydro_py \
    --volumes-from croton \
    wrfhydro/dev:conda
# To open an additional terminal to this docker image, this assumes that the
# above is the only hash returned by $(docker ps -q)
docker exec -it $(docker ps -q) bash
# Inside docker:
# This is necessary until the pip dist of wrfhydropy contains all the code.
cd ~/wrf_hydro_py/
pip uninstall -y wrfhydropy
python setup.py develop
# Run the script
cd ~/wrf_hydro_dart/models/wrf_hydro/shell_scripts
python setup_experiment.py hydro_dart_config.yaml


# ###################################
# Cheyenne
# 1. Get the domain.
#    Domains with the correct conventions to work with wrfhydropy are kept here on
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
#    If a development version of wrfhydropy is needed, you'll need to clone that repository to cheyenne, then do the following
#       cd /path/to/wrf_hydro_py/; pip uninstall -y wrfhydropy; python setup.py develop

#rm -rf /glade/scratch/jamesmcc/dart_runs/test_run 
#python setup_experiment.py hydro_dart_config_cheyenne.yaml 
