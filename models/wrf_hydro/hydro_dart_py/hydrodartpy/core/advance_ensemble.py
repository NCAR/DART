import argparse
import os
import pickle
import sys
import wrfhydropy


# advance_ensemble --help
# Will print the help.
parser = argparse.ArgumentParser(
    description='Advance a HydroDartRun Ensemble'
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


args = parser.parse_args()

hydro_dart_run = pickle.load(open(args.run_dir + '/experiment_dir/HydroDartRun.pkl', 'rb'))

hydro_dart_run.advance_ensemble(
    model_start_time=args.model_start_time,
    model_end_time=args.model_end_time,
    job_entry_cmd=args.job_entry_cmd,
    job_exit_cmd=args.job_exit_cmd,
    afterok=args.afterok,
    afterany=args.afterany
)

sys.exit()
