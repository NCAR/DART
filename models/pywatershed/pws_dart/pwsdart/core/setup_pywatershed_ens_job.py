import copy
import numpy as np
import os, sys
import shutil
import subprocess
import xarray as xr
from pprint import pprint
import pywatershed as pws
import re
from pathlib import Path

#To DO:
# at the setup time, we want to create a file .model_end_time and out the initial time in the file. 
# mke symlink to the domain directory 
# make symlinks to the forcing directory (if ensemble, then each member forcing directory should be to the member ensembl) FORCING_perturbed
# sy,ink to the parameter that we generated 
# copy the config control_file
# modify the control file to match the time (intial time), domain, forcing, parameter, restart file, duration of the simulation 

# For now we are only perturn the manning n 

# This file creates
# 1. channel parameter distributons with physical constraints.
# 2. GW parameter distribution based on calibration stats.

def setup_pywatershed_job(config):
    run_dir = Path(config["experiment"]["run_dir"])
    ensemble_size = config["ensemble"]["size"]
    domain_dir = Path(config["pywatershed"]["domain_dir"])
    control_file_name = config["pywatershed"]["control_file"]
    param_dir = Path(config["ensemble"]["path"])

    forcing_base = Path(config['forcing_ensemble']['path'])

    for i in range(ensemble_size):
        member_name = f"member_{i:03d}"
        member_dir = run_dir / member_name
        if not member_dir.exists():
            raise FileNotFoundError(f"Missing member directory: {member_dir}")

        # 1. Extract timestamp from restart and write .model_end_time
        restart_files = list(member_dir.glob("*.nc"))
        for file in restart_files:
            match = re.search(r"(\d{4}-\d{2}-\d{2})", file.name)
            if match:
                timestamp = f"{match.group(1)} 00:00:00"  #could the timestamp be changed?
                (member_dir / ".model_end_time").write_text(timestamp)
                print(f"{member_name}: wrote .model_end_time = {timestamp}")
                break
        else:
            raise ValueError(f"No date found in restart filenames in {member_dir}")
            # Do we need to udate model end time to something else if this is the case? Or this is never the case

        # 2. Symlink DOMAIN/
        domain_link = member_dir / "DOMAIN"
        if not domain_link.exists():
            domain_link.symlink_to(domain_dir.resolve(), target_is_directory=True)
            print(f"{member_name}: linked DOMAIN → {domain_dir}")

        # 3. Symlink FORCING/
        forcing_member_dir = forcing_base / member_name
        if not forcing_member_dir.exists():
            raise FileNotFoundError(f"{member_name}: Missing forcing dir at {forcing_member_dir}")
        forcing_link = member_dir / "FORCING"
        if not forcing_link.exists():
            forcing_link.symlink_to(forcing_member_dir.resolve(), target_is_directory=True)
            print(f"{member_name}: linked FORCING → {forcing_member_dir}")

        # 4. Symlink perturbed parameter file
        param_src = param_dir / f"perturbed_params_{i:03d}.nc"
        param_dst = member_dir / "perturbed_params.nc"
        if not param_src.exists():
            raise FileNotFoundError(f"{member_name}: missing perturbed param file {param_src}")
        if not param_dst.exists():
            param_dst.symlink_to(param_src.resolve())
            print(f"{member_name}: linked perturbed_params.nc → {param_src}")

    
        # 5. Copy and patch nhm.control
        control_src = domain_dir / control_file_name
        control_dst = member_dir / "nhm.control"
        
        shutil.copy(control_src, control_dst)
        

    print("✅ Setup complete for all ensemble members.")
