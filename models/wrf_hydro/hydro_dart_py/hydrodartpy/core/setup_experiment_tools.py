from boltons.iterutils import remap
import os
import pathlib
import re
import socket
import warnings
import yaml


# ######################################################
# Remapping nested values
# http://sedimental.org/remap.html

def visit(path, key, value):
    if isinstance(value, str):
        return key, os.path.expanduser(os.path.expandvars(value))
    return key, value


def remap_spec(spec_file):
    return(remap(spec_file, visit))


# ######################################################
# Generic spec establishment = YAML + remap_spec

def establish_config(spec_file):
    """Parse YAML and expand  ~ and $
    """
    with open(spec_file) as ff:
        spec_dict = yaml.safe_load(ff)

    spec = remap_spec(spec_dict)

    return(spec)


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

