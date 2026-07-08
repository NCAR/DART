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
# Perturbed forcings for pywatershed

def setup_pywatershed_forcings(config):
    # Input
    n_mem = config['ensemble']['size']
    out_dir = config['forcing_ensemble']['path']
    forcing_dir = pathlib.Path(config['forcing_ensemble']['perturb_forcing']['forcing_dir'])

    create_flag = config['forcing_ensemble']['perturb_forcing']['create']
    forcing_perturb_vars_config = config["forcing_ensemble"]["perturb_forcing"]["forcing_vars_to_perturb"]

    if create_flag == True:
    # -------------------------------------------------------
        if out_dir.exists():
            print(f'The directory exists: {out_dir}')
            sys.exit(1)
        else:
            os.mkdir(str(out_dir))


        # List of all forcing .nc files in input dir
        all_forcing_files = sorted(forcing_dir.glob("*.nc"))
        all_forcing_vars = [f.stem for f in all_forcing_files]

        forcing_vars_to_perturb = config['forcing_ensemble']['perturb_forcing']['forcing_vars_to_perturb']['param_name']
       
        # Vars to copy (not perturbed)
        vars_to_copy = [v for v in all_forcing_vars if v not in forcing_vars_to_perturb]

        for var in forcing_vars_to_perturb:
            forcing_file = str(forcing_dir) + "/" + var + ".nc"

            for ii in range(n_mem):
                forcing_ds = xr.open_dataset(forcing_file)

                for name, func, min_val, max_val in zip(forcing_perturb_vars_config['param_name'],forcing_perturb_vars_config['perturb_function'],forcing_perturb_vars_config['min'],forcing_perturb_vars_config['max']):
                    if func == 'uniform':
                        value = np.random.uniform(min_val, max_val)
                    #elif func == ''
                    else:
                        raise ValueError(f"Unsupported perturb_function: {func}")
            
                forcing_ds[var] = forcing_ds[var] * value
        
                forcing_member_dir = (out_dir/f"member_{str(ii).zfill(3)}")
                forcing_member_dir.mkdir(parents=True, exist_ok=True)
                forcing_file_name = (
                    out_dir / f"member_{str(ii).zfill(3)}/{var}.nc"
                )
                
                forcing_ds.to_netcdf(forcing_file_name)

        for var in vars_to_copy:
            src_file = forcing_dir / f"{var}.nc"
            for ii in range(n_mem):
                member_dir = out_dir / f"member_{str(ii).zfill(3)}"
                target = os.path.abspath(src_file)
                link_name = member_dir / f"{var}.nc"
                os.symlink(target, link_name)
