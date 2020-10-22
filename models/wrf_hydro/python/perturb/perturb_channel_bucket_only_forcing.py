import argparse
import datetime
import f90nml
import importlib.util
import os
import pathlib
import pickle
import socket
import sys
import time
import xarray as xa


def import_noise_model(noise_file):
    spec = importlib.util.spec_from_file_location("noise_model", noise_file)
    noise_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(noise_module)
    return noise_module.noise_model


def perturb_channel_only_forcing(
    chrtout_file,
    qsfclat_perturb_func,
    qbtmvert_perturb_func,
    out_dir,
    seed: int=None
):

    # Input.
    # Open the single forcing file.
    drop_list = ['streamflow', 'nudge', 'q_lateral', 'velocity']
    chrtout = xa.open_dataset(chrtout_file, drop_variables=drop_list)

    # Apply the perturbations. This set maintains metadata.
    chrtout.qSfcLatRunoff.values = qsfclat_perturb_func(chrtout.qSfcLatRunoff.values, seed=seed)
    chrtout.qBtmVertRunoff.values = qbtmvert_perturb_func(chrtout.qBtmVertRunoff.values, seed=seed)

    # This is a check that the perturbations are not correlated. This needs run when
    # only the perturbations are returned (withouth the mean values)
    # import numpy as np
    # X = np.stack((chrtout.qSfcLatRunoff.values, chrtout.qBtmVertRunoff.values), axis=0)
    # np.corrcoef(X)
    
    # Output.
    # Annoyingly, xarray wants to force NaN on you for _FillValue.
    chrtout.qSfcLatRunoff.encoding.update({'_FillValue': None})
    chrtout.qBtmVertRunoff.encoding.update({'_FillValue': None})
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
        '--qbtmvert_perturb_func',
        metavar='qbtm_vert_perturb_func',
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
    qbtmvert_perturb_func = args.qbtmvert_perturb_func[0]

    seed = args.seed

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
    qsfclat_perturb_func = import_noise_model(qsfclat_perturb_func)
    qbtmvert_perturb_func = import_noise_model(qbtmvert_perturb_func)

    # Loop over the forcing files
    for bb in forcing_basenames:

        if seed is not None:
            if not isinstance(seed, int):
                seed = seed[0]
        else:
            member_number = int(pathlib.Path(run_dir).name.split('_')[1])
            the_time = datetime.datetime.strptime(bb.split('.')[0], '%Y%m%d%H%M')
            greg_time = the_time - datetime.datetime(1601,1,1)

            # Idea 1
            # cyclic every 27.397 years (10,000 days)
            # minutes is finest time resolution
            # assuming < 100 ens members.
            # greg_min = greg_time.seconds / 60
            # seed = (greg_time.days % 10000) * 1000000 + \
            #                        greg_min *     100 + \
            #                   member_number

            # Idea 2
            # cyclic every 19.29 years (10,000,00 hours)
            # minutes is finest time resolution
            # assuming < 100 ens members.
            greg_time_min = greg_time.days * 24*60 + greg_time.seconds / 60
            seed = int(((greg_time_min %10000000) * 100) + member_number)
            
        perturb_channel_only_forcing(
            forcing_dir_existing / bb,
            qsfclat_perturb_func,
            qbtmvert_perturb_func,
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
