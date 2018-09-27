from boltons.iterutils import remap
import os
from pathlib import Path
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
