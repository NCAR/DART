import argparse
import datetime
import f90nml
import os
import pathlib
import pickle
import pickle
import sys
import yaml
from wrfhydropy.core.ensemble_tools import get_ens_dotfile_end_datetime

parser = argparse.ArgumentParser(
    description='Get the time of a HydroDartRun Ensemble'
)

parser.add_argument(
    '--run_dir',
    required=False,
    metavar='/abs/path/to/run/directory',
    help='Path to the experiment directory. (Default is director where this script or a ' +
         'symlink (not resolved) to it lives).',
    default= os.path.dirname(os.path.abspath(__file__))
)

args = parser.parse_args()
run_dir = pathlib.PosixPath(args.run_dir)

valid_time = get_ens_dotfile_end_datetime(run_dir)

# Read the YAML config file.
experiment_dir = run_dir / 'experiment_dir'
config_file = sorted(experiment_dir.glob("*.original.*.yaml"))[0]
with open(config_file) as ff:
    config = yaml.safe_load(ff)

timedelta = datetime.timedelta
config_run_time = config['run_experiment']['time']
assim_window_hr = int(config_run_time['assim_window']['assim_window_size_hours'])
assim_half_window_sec = timedelta(seconds=(60*60/2)*assim_window_hr)

start_window = valid_time - assim_half_window_sec + timedelta(seconds=1)
end_window = valid_time + assim_half_window_sec

start_window_greg = start_window - datetime.datetime(1601, 1, 1)
end_window_greg = end_window - datetime.datetime(1601, 1, 1)

input_nml_file = run_dir / 'input.nml'
input_nml = f90nml.read(input_nml_file)

obs_seq_nml = input_nml['obs_sequence_tool_nml']
filter_nml = input_nml['filter_nml']

if (run_dir / 'obs_sequence_tool').exists():

    print('Using obs_sequence_tool.')
    obs_seq_nml['filename_seq'] = 'obs_seq.daily'
    obs_seq_nml['filename_out'] = 'obs_seq.window'
    obs_seq_nml['first_obs_days'] = start_window_greg.days
    obs_seq_nml['first_obs_seconds'] = start_window_greg.seconds
    obs_seq_nml['last_obs_days'] = end_window_greg.days
    obs_seq_nml['last_obs_seconds'] = end_window_greg.seconds

    filter_nml['obs_sequence_in_name'] = 'obs_seq.window'

else:

    print('Not using obs_sequence_tool.')
    filter_nml['obs_sequence_in_name'] = 'obs_seq.daily'
    filter_nml['first_obs_days'] = start_window_greg.days
    filter_nml['first_obs_seconds'] = start_window_greg.seconds
    filter_nml['last_obs_days'] = end_window_greg.days
    filter_nml['last_obs_seconds'] = end_window_greg.seconds


input_nml.write(input_nml_file, force=True)

sys.exit()
