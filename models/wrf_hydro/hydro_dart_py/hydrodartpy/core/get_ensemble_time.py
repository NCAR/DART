import argparse
import os
import pickle
import sys
from wrfhydropy.core.dartclasses import get_ens_last_restart_datetime

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
ens_run = pickle.load(open(args.run_dir + '/WrfHydroEnsembleRun.pkl', 'rb'))
ens_run.collect_ensemble_runs() # inlieu of a sweeper job # this is printing the annoying 0
print(get_ens_last_restart_datetime(ens_run).strftime('%Y-%m-%d_%H:%M'))

sys.exit()
