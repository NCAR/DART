import copy
import numpy as np
import os, sys
import shutil
import subprocess
import xarray as xr
import pathlib 
from pprint import pprint
import pywatershed as pws


#To DO:
# For now we are only perturn the manning n 

# This file creates
# 1. channel parameter distributons with physical constraints.
# 2. GW parameter distribution based on calibration stats.

def setup_pywatershed_ens(config):
    # Input
    n_mem = config['ensemble']['size']
    domain_dir = config['pywatershed']['domain_dir']
    param_path = domain_dir + '/' + config['pywatershed']['param_file']

    # Output
    out_dir = config['ensemble']['path']
    create_flag = config['ensemble']['create']

    if create_flag == True:
    # -------------------------------------------------------
        if out_dir.exists():
            print(f'The directory exists: {out_dir}')
            sys.exit(1)
        else:
            os.mkdir(str(out_dir))

        # PRMS native parameter file 
        params = pws.parameters.PrmsParameters.load(param_path)

        param_files = []  # get a list of written NetCDF files back at the end
        for ii in range(n_mem):
            param_ds = params.to_xr_ds()  # copies by default
            perturb_cfg = config['ensemble']['params_to_perturb']
            for name, func, min_val, max_val in zip(perturb_cfg['param_name'],perturb_cfg['perturb_function'], perturb_cfg['min'],perturb_cfg['max']):
                if func == 'uniform':
                    value = np.random.uniform(min_val, max_val)
                #elif func == ''
                else:
                    raise ValueError(f"Unsupported perturb_function: {func}")
            
                param_ds[name] *= value
                print(param_ds[name])
        
            param_file_name = (
                out_dir / f"perturbed_params_{str(ii).zfill(3)}.nc"
            )
            param_files += [param_file_name]
            param_ds.to_netcdf(param_file_name)
