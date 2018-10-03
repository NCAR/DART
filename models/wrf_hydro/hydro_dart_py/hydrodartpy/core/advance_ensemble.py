import argparse
import multiprocessing
import os
import pathlib
import pickle
import shutil
import sys
from wrfhydropy.core.ensemble_tools import mute
import yaml

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

parser.add_argument(
    '--hold',
    required=False,
    metavar='holdflag',
    type=bool,
    help='Flag to hold the job array. (Default = False).',
    default=False
)


args = parser.parse_args()


def parallel_update_job(member_dir):

    orig_dir = os.getcwd()
    os.chdir(str(member_dir))

    shutil.rmtree('previous_job_ens_adv', ignore_errors=True)
    shutil.move('job_ens_adv', 'previous_job_ens_adv')
    os.mkdir('job_ens_adv')

    sim = pickle.load(open("WrfHydroSim.pkl", "rb"))
    job = sim.jobs[-1]

    new_start_time = job.model_end_time
    new_end_time = new_start_time + (job.model_end_time - job.model_start_time)

    print('new_start_time: ', new_start_time)
    print('new_end_time: ', new_end_time)

    job._model_start_time = new_start_time
    job._model_end_time = new_end_time
    job._set_hrldas_times()
    job._set_hydro_times()
    job._write_namelists()

    sim.jobs = [job]
    sim.pickle('WrfHydroSim.pkl')

    os.chdir(orig_dir)


run_dir = pathlib.Path(args.run_dir)
experiment_dir = run_dir / 'experiment_dir'
config_file = sorted(experiment_dir.glob("*.original.*.yaml"))[0]
with open(config_file) as ff:
    config = yaml.safe_load(ff)

# Run the current ensemble
n_concurrent = config['run_experiment']['dart']['scheduler']['ppn_max'] - 1
e = pickle.load(open('WrfHydroEnsembleSim.pkl', 'rb'))
e.run(n_concurrent=n_concurrent)

# Prep the next
pool = multiprocessing.Pool(n_concurrent, initializer=mute)
member_dirs = sorted(pathlib.Path('.').glob("member_*"))

if n_concurrent > 1:
    pool.map(parallel_update_job, (mm for mm in member_dirs))
else:
    for mm in member_dirs:
        parallel_update_job(mm)

sys.exit()
