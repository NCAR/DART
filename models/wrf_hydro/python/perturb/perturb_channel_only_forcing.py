import argparse
import datetime
import f90nml
import os
import pathlib
import pickle
import socket
import sys
import time
import xarray as xa

def perturb_channel_only_forcing(
    chrtout_file,
    qsfclat_perturb_func,
    qbucket_perturb_func,
    out_dir,
    seed: int=None
):

    # Input.
    # Open the single forcing file.
    drop_list = ['streamflow', 'nudge', 'q_lateral', 'velocity']
    chrtout = xa.open_dataset(chrtout_file, drop_variables=drop_list)

    # Apply the perturbations. This set maintains metadata.
    
    chrtout.qSfcLatRunoff.values = qsfclat_perturb_func(chrtout.qSfcLatRunoff.values, seed=seed)
    chrtout.qBucket.values = qbucket_perturb_func(chrtout.qBucket.values, seed=seed)

    # Output.
    # Annoyingly, xarray wants to force NaN on you for _FillValue.
    chrtout.qSfcLatRunoff.encoding.update({'_FillValue': None})
    chrtout.qBucket.encoding.update({'_FillValue': None})
    chrtout_basename = os.path.basename(chrtout_file)
    paths = out_dir / chrtout_basename

    chrtout.to_netcdf(paths, mode='w')
    return paths


if __name__ == "__main__":

    #ts = time.time()

    parser = argparse.ArgumentParser(
        description='Apply noise to channel-only forcing files implied by the current run.'
    )

    parser.add_argument(
        '--run_dir',
        metavar='run_dir',
        type=str,
        nargs=1,
        default=None
    )

    parser.add_argument(
        '--source_forcing_dir',
        metavar='source_forcing_dir',
        type=str,
        nargs=1,
        help='This dir should exist in the run_dir.'
    )

    parser.add_argument(
        '--perturb_forcing_dir',
        metavar='perturb_forcing_dir',
        type=str,
        nargs=1,
        help='This should just be the basename, the dir will be placed in the run_dir.'
    )

    parser.add_argument(
        '--qsfclat_perturb_func',
        metavar='qsfclat_perturb_func',
        type=str,
        nargs=1,
        help='A /path/to/file.py which contains a single function noise_model(). Which takes' + \
             ' a forcing file and applies noise to this specific variable.'
    )

    parser.add_argument(
        '--qbucket_perturb_func',
        metavar='qbucket_perturb_func',
        type=str,
        nargs=1,
        help='As for qsfclat_perturb_func.'
    )

    parser.add_argument(
        '--seed',
        metavar='seed_for_np.random',
        type=int,
        nargs=1,
        help='seed to control random variables.',
        default=None
    )

    args = parser.parse_args()

    run_dir = args.run_dir
    if run_dir is None:
        run_dir = os.getcwd() #os.path.dirname(__file__)
    else:
        run_dir = run_dir[0]

    source_forcing_dir = args.source_forcing_dir[0]
    perturb_forcing_dir = args.perturb_forcing_dir[0]
    qsfclat_perturb_func = args.qsfclat_perturb_func[0]
    qbucket_perturb_func = args.qbucket_perturb_func[0]

    seed = args.seed
    if seed is not None:
        seed = seed[0]

    # Note: Use the generic namelists in the run dir not from the object.
    # Now that wrfhydropy only stages the namelists on run, I resort to this hack
    # which relies on using the last job_* dir
    job_dir = sorted(pathlib.Path(run_dir).glob('job_*'))[-1]
    namelist_hrldas = f90nml.read(str(job_dir) + '/namelist.hrldas')

    # How many forcing files for the duration?
    # Kday or khour?
    ktime_key = set(namelist_hrldas['noahlsm_offline'].keys()) & set(['kday', 'khour'])
    ktime_key = list(ktime_key)[0]
    ktime_hour = int(namelist_hrldas['noahlsm_offline'][ktime_key])
    if ktime_key == 'kday':
        ktime_hour = ktime_hour * 24


    forcing_timestep = int(namelist_hrldas['noahlsm_offline']['forcing_timestep'])
    forcing_timestep_hour = forcing_timestep / (60*60)
    n_forcing_files = ktime_hour / forcing_timestep_hour
    if n_forcing_files != int(n_forcing_files):
        raise ValueError('ktime_hour / forcing_timestep_hour does not give integer')
    n_forcing_files = int(n_forcing_files)

    noahlsm_offline = namelist_hrldas['noahlsm_offline']
    start_year = int(noahlsm_offline['start_year'])
    start_month = int(noahlsm_offline['start_month'])
    start_day = int(noahlsm_offline['start_day'])
    start_hour = int(noahlsm_offline['start_hour'])
    start_min = int(noahlsm_offline['start_min'])

    model_time = datetime.datetime(
        start_year,
        start_month,
        start_day,
        start_hour,
        start_min
    )

    # YYYYmmddHHMM.CHRTOUT_DOMAIN1
    forcing_basenames = \
        [(model_time + datetime.timedelta(hours=(tt+1))).strftime('%Y%m%d%H%M') + '.CHRTOUT_DOMAIN1'
         for tt in range(n_forcing_files)]

    # Current forcing directory
    # Take this from the domain object NOT the namelist
    forcing_dir_existing = pathlib.PosixPath(run_dir) / source_forcing_dir

    # Create new forcing dir if it does not exist.
    forcing_dir_perturbed = pathlib.PosixPath(run_dir) / perturb_forcing_dir
    forcing_dir_perturbed.mkdir(parents=False, exist_ok=True)

    # Get the noise/perturbation functions
    exec(open(qsfclat_perturb_func).read())
    qsfclat_perturb_func = noise_model()

    exec(open(qbucket_perturb_func).read())
    qbucket_perturb_func = noise_model()

    # Loop over the forcing files
    for bb in forcing_basenames:
        perturb_channel_only_forcing(
            forcing_dir_existing / bb,
            qsfclat_perturb_func,
            qbucket_perturb_func,
            forcing_dir_perturbed,
            seed
        )

    # replace with the new/perturbed forcing dir
    # DO NOT change the run_obj.setup.domain.forcing_dir
    namelist_hrldas['noahlsm_offline']['indir'] = perturb_forcing_dir
    namelist_hrldas.write(str(job_dir) + '/namelist.hrldas', force=True)

    #te = time.time()
    #print('Timing: ' + str(round(te - ts, 2)) + ' seconds: ')
    #print(datetime.datetime.now())
    #print(socket.gethostname())

    sys.exit(0)
