import argparse
from datetime import datetime, timedelta 
import multiprocessing
import os
import pathlib
import pickle
import shutil
import sys
import yaml
import pywatershed as pws
import numpy as np
from pwsdart.core.get_ensemble_time import get_ensemble_time

def run_member(member_dir, start_time, end_time):
    domain_dir = member_dir / "DOMAIN"
    param_file = member_dir / "perturbed_params.nc"
    params = pws.parameters.PrmsParameters.from_netcdf(param_file)

    control = pws.Control.load_prms(
              domain_dir / "nhm.control", warn_unused_options=False
               )

    control.options["input_dir"] = member_dir / "FORCING"
    control.options["restart_write"] = member_dir
    control.options["restart_write_freq"] = "d"
    control.options["restart_read"] = member_dir
    del control.options["netcdf_output_dir"]
    del control.options["netcdf_output_var_names"]
    control.edit_end_time(end_time)
    control.edit_init_start_times(start_time)

    control.options = control.options | {
     "budget_type": "warn",
     "calc_method": "numba",
     }

    nhm_processes = [
      pws.PRMSGroundwater,
      pws.PRMSChannel,
    ]

    model = pws.Model(
            nhm_processes,
            control=control,
            parameters=params,
            )

    model.run(finalize=True)

    # after the run is complete, change the date in the .model_end_time 
    dt = end_time.astype('datetime64[s]').tolist()
    timestamp = dt.strftime("%Y-%m-%d %H:%M:%S")
    (member_dir / ".model_end_time").write_text(timestamp)
    return str(member_dir)
 
def advance_ensemble(
    run_dir: pathlib.Path,
    advance_n_hours: int=None,
    ncores : int=None,
    ncores_available : int=None
):

    run_dir = pathlib.Path(run_dir)
    experiment_dir = run_dir / 'experiment_dir'
    member_dirs = sorted(pathlib.Path('.').glob("member_*"))

    # Run the current ensemble
    start_time = np.datetime64(get_ensemble_time(run_dir))
    end_time = np.datetime64(start_time + timedelta(hours=advance_n_hours))

    print(start_time)
    print(end_time)

    tasks = [(member_dir, start_time, end_time) for member_dir in member_dirs]

    with multiprocessing.Pool(processes=ncores_available) as pool:
        results = pool.starmap(run_member, tasks)

    print("Ensemble members advanced:", results)

if __name__ == "__main__":

    # advance_ensemble --help
    parser = argparse.ArgumentParser(
        description='Advance a model Ensemble'
    )

    parser.add_argument(
        '--run_dir',
        required=False,
        metavar='/abs/path/to/run/directory',
        help='Path to the experiment directory. (Default is director where this script or a ' +
             'symlink (not resolved) to it lives).',
        default= os.path.dirname(os.path.abspath(__file__))
    )

    parser.add_argument(
        '--model_start_time',
        required=False,
        metavar='YYYY-mm-dd HH:MM',
        type=str,
        help='The start time string as YYYY-mm-dd HH:MM defaults (set at run time to' +
             'last restart time available).',
        default=None
    )

    parser.add_argument(
        '--model_end_time',
        required=False,
        metavar='YYYY-mm-dd HH:MM',
        type=str,
        help="The end time string as YYYY-mm-dd HH:MM set at run time to" +
             "last model_start_time + config_file['run_experiment']['time']['window_size_hours']).",
        default=None
    )

    parser.add_argument(
        '--job_entry_cmd',
        required=False,
        metavar='entry-command',
        type=str,
        help='The command to run in run_dir before advancing.',
        default=None
    )

    parser.add_argument(
        '--job_exit_cmd',
        required=False,
        metavar='exit-command',
        type=str,
        help='The command to run in run_dir after advancing.',
        default=None
    )

    parser.add_argument(
        '--afterok',
        required=False,
        metavar='jobid',
        type=str,
        help='The job id of the script which must successfully complete before the advance.' +
             ' (Default = None).',
        default=None
    )

    parser.add_argument(
        '--afterany',
        required=False,
        metavar='jobid',
        type=str,
        help='The job id of the script which must complete before the advance. (Default = None).',
        default=None
    )

    parser.add_argument(
        '--hold',
        required=False,
        metavar='holdflag',
        type=bool,
        help='Flag to hold the job array. (Default = False).',
        default=False
    )

    args = parser.parse_args()

    config_file = sorted(experiment_dir.glob("original.*.yaml"))[0]
    # TODO this should use establish_config
    with open(config_file) as ff:
        config = yaml.safe_load(ff)

    advance_ensemble(
        run_dir=args.run_dir
    )

    # TODO: only one arg is being used... ?
    # Apparently the initial job has the stuff I thought I was passing here.
    
    sys.exit(0)
