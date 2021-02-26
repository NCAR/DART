from boltons.iterutils import remap
import os
import pathlib
import re
import socket
import warnings
import yaml


this_file = pathlib.Path(__file__)
repo_dir = this_file.parent.parent.parent.parent.parent.parent


# ######################################################
# Remapping nested values
# http://sedimental.org/remap.html

def visit_expand_vars(path, key, value):
    if isinstance(value, str):
        return key, os.path.expanduser(os.path.expandvars(value))
    return key, value


config_abs_paths_list = [
    'experiment_dir',
    'run_dir',
    'path',
    'hydro_file_list.txt',
    'param_file_list.txt',
    'dart_src',
    'wrf_hydro_src',
    'domain_src',
    'setup_py',
    'noise_function',
    'output_dir',
    'input_dir'
]


def visit_abs_paths(path, key, value):
    if value is None:
        return True
    # Making a bet that the same keys at different hierarchical levels
    # will both be desired as pathlib.PosixPath objects.
    if key in config_abs_paths_list and type(value) is not dict:
        return key, pathlib.PosixPath(value)
    else:
        return True


def visit_abs_paths_to_str(path, key, value):
    if value is None:
        return True
    # Making a bet that the same keys at different hierarchical levels
    # will both be desired as pathlib.PosixPath objects.
    if key in config_abs_paths_list and type(value) is not dict:
        return key, str(value)
    else:
        return True


# ######################################################
# Generic spec establishment = YAML + remap_spec

def establish_config(spec_file):
    """1) Parse YAML,
       2) expand  ~ and $ using the system,
       3) turn selected strings into pathlib.Path objects
    """
    # 1) Parse YAML
    with open(spec_file) as ff:
        config0 = yaml.safe_load(ff)

    # 2) expand the ~ and  vars
    config1 = remap(config0, visit_expand_vars)

    # 3) Selected strings to pathlib.Path objects.
    # The following list and function turn strings in the config
    # into pathlib.PosixPath objs.
    config2 = remap(config1, visit_abs_paths)

    return(config2)


# #######################################################
# Get top-level domain dir from config.

def get_top_level_dir_from_config(config, sim):
    """Get top-level domain dir from config."""
    dom_config = config['wrf_hydro']['model_config']
    # I  like the former, commented appraoch b/c the later might have a
    # modified path. It seems like the solution is to stash a piece of metadata with
    # the sim since no in dividual namelist entry is guaranteed to be unchanged.
    #patches = sim.domain.hydro_namelist_patches
    #geo_file = patches['hydro_nlist']['geo_static_flnm']
    geo_file = sim.base_hydro_namelist['hydro_nlist']['geo_finegrid_flnm']
    top_level_dir = geo_file.split('/')[1]

    return top_level_dir


# #######################################################
# sed-like
def replace_in_file(file, str_find, str_replace):
    str_find = str(str_find)
    str_replace = str(str_replace)
    with open(file, "r") as sources:
        lines = sources.readlines()
    with open(file, "w") as sources:
        for line in lines:
            _ = sources.write(re.sub(str_find, str_replace, line))


# #######################################################
# Detect the host machine, generically.
def get_machine():
    hostname = socket.gethostname()
    if re.match('cheyenne', hostname):
        machine = 'cheyenne'
    else:
        # Need logic for geyser?
        if pathlib.PosixPath('/glade/').exists():
            machine = 'cheyenne'
        else:
            machine = 'docker'
            #if not in_docker():
            #    warnings.warn('This machine is not recognized, using docker defaults.')
    return machine

