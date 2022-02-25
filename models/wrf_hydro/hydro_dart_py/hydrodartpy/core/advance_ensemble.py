import argparse
import datetime
import multiprocessing
import os
import pathlib
import pickle
import shutil
import sys
from wrfhydropy.core.ensemble_tools import mute
import yaml


def parallel_update_job(arg_dict):
    # This function currently does 2 things
    # 1) if 'n_hours': is in the arg_dict: if extends an existing job to run that long.
    # 2) else: delete the existing job but use it's end time and khour to make a new job.
    member_dir = arg_dict['member_dir']

    orig_dir = os.getcwd()
    os.chdir(str(member_dir))

    sim = pickle.load(open("WrfHydroSim.pkl", "rb"))
    job = sim.jobs[-1]

    if 'n_hours' in arg_dict.keys():
        # Just adjust the existing job.
        n_hours = arg_dict['n_hours']
        new_end_time = job.model_start_time + datetime.timedelta(hours=n_hours)
        print('new_end_time: ', new_end_time)
        job._model_end_time = new_end_time
        job._write_namelists(mode='w')
        job._write_namelists(mode='w')

    else:
        # Create a fresh job.
        shutil.rmtree('previous_job_ens_adv', ignore_errors=True)
        shutil.move('job_ens_adv', 'previous_job_ens_adv')
        os.mkdir('job_ens_adv')

        new_start_time = job.model_end_time
        new_end_time = new_start_time + (job.model_end_time - job.model_start_time)

        print('new_start_time: ', new_start_time)
        print('new_end_time: ', new_end_time)

        job._model_start_time = new_start_time
        job._model_end_time = new_end_time
        job._write_namelists()

    sim.jobs = [job]
    sim.pickle('WrfHydroSim.pkl')

    os.chdir(orig_dir)


def advance_ensemble(
    run_dir: pathlib.Path,
    n_concurrent: int=None,
    teams_dict: dict=None,
    advance_n_hours: int=None
):

    run_dir = pathlib.Path(run_dir)
    experiment_dir = run_dir / 'experiment_dir'
    member_dirs = sorted(pathlib.Path('.').glob("member_*"))

    # Run the current ensemble
    e = pickle.load(open('WrfHydroEnsembleSim.pkl', 'rb'))

    if advance_n_hours is not None:

        if n_concurrent > 1:
            with multiprocessing.Pool(n_concurrent, initializer=mute) as pool:
                _ = pool.map(
                    parallel_update_job,
                    ({'member_dir': mm, 'n_hours': advance_n_hours} for mm in member_dirs)
                )
        else:
            _ = [
                parallel_update_job({'member_dir': mm, 'n_hours': advance_n_hours})
                for mm in member_dirs
            ]

    e.run(n_concurrent=n_concurrent,
          teams=(teams_dict is not None),
          teams_dict=teams_dict)

    # Prep the next
    if n_concurrent > 1:
        with multiprocessing.Pool(n_concurrent, initializer=mute) as pool:
            _ = pool.map(
                parallel_update_job,
                ({'member_dir': mm} for mm in member_dirs)
            )
    else:
        _ = [
            parallel_update_job({'member_dir': mm})
            for mm in member_dirs
        ]


if __name__ == "__main__":

    # advance_ensemble --help
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

    config_file = sorted(experiment_dir.glob("original.*.yaml"))[0]
    # TODO this should use establish_config
    with open(config_file) as ff:
        config = yaml.safe_load(ff)

    advance_ensemble(
        run_dir=args.run_dir,
        n_concurrent=config['run_experiment']['dart']['scheduler']['ppn_max'] - 1
    )

    # TODO: only one arg is being used... ?
    # Apparently the initial job has the stuff I thought I was passing here.
    
    sys.exit(0)
