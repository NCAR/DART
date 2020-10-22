#!/usr/bin/env python

import argparse
import dask
from multiprocessing.pool import Pool
import pathlib
import sys
import time
from typing import Union
import wrfhydropy
import xarray as xr

def collect_run_output(
    run_dir: Union[pathlib.Path,str],
    n_cores: int,
    spatial_indices: list=None,
    drop_variables: list=None
):

    run_dir = pathlib.Path(run_dir)

    t_start = time.time()

    # The DART_cleanup_add_time.csh step can be done in the preprocess steps below.
    # JLM: I dont understand why the member files dont have "out" equivalents with time...

    # -------------------------------------------------------
    # 1. These files only concatenate along the time dimension.
    #    We can use xr.open_mfdataset out of the box!
    stage_list = ['input', 'preassim', 'analysis', 'output']

    type_base_list = ['mean', 'sd', 'priorinf_mean', 'priorinf_sd', 'postinf_mean', 'postinf_sd']
    domain_list = ['d01', 'd02']
    type_list = type_base_list
    for domain in domain_list:
        type_list = type_list + [(typ + '_' + domain) for typ in type_base_list]

    for stage in stage_list:
        for typ in type_list:

            # This is the correct glob
            #in_files = sorted((run_dir / 'output').glob('*/' + stage + '_' + typ + '.*.nc'))
            # This one is more restrictive b/c of the DART_cleanup *out* files.
            in_files = sorted((run_dir / 'output').glob('*/' + stage + '_' + typ + '.*[0-9].nc'))
            if len(in_files) == 0:
                continue

            out_file = run_dir / ('all_' + stage + '_' + typ + '.nc')

            # Do have to add the time dim to each variable to get the correct result.
            def preproc_time(ds):
                for key in ds.variables.keys():
                    if 'time' not in ds[key].dims:
                        ds[key] = ds[key].expand_dims('time')
                return ds

            the_pool=Pool(n_cores)
            with dask.config.set(scheduler='processes', pool=the_pool):
                ds = xr.open_mfdataset(in_files, parallel=True, preprocess=preproc_time)
            the_pool.close()
            # Feel like this should go in the above with/context. But it errors.
            # Xarray sets nan as the fill value when there is none. Dont allow that...
            for key, val in ds.variables.items():
                if '_FillValue' not in ds[key].encoding:
                    ds[key].encoding.update({'_FillValue': None})
            ds.to_netcdf(out_file)
            del ds

    #-------------------------------------------------------
    # 2. Collect members. This replaces DART_cleanup_pack_members.csh and DART_cleanup.csh
    #    The explicit handling of individual members happens in wrfhydropy.open_dart_dataset
    stage_list = ['input', 'preassim', 'analysis', 'output']
    domain_list = ['', '_d0', '_d01', '_d02']

    for stage in stage_list:
        for domain in domain_list:

            # This is the correct glob
            #in_files = sorted((run_dir / 'output').glob('*/' + stage + '_' + typ + '.*.nc'))
            # This one is more restrictive b/c of the DART_cleanup *out* files.
            in_files = sorted(
                (run_dir / 'output').glob(
                    '*/' + stage + '_member_*' + domain + '.*[0-9].nc'
                )
            )
            if len(in_files) == 0:
                continue

            the_pool=Pool(n_cores)
            with dask.config.set(scheduler='processes', pool=the_pool):
                ds = wrfhydropy.ioutils.open_dart_dataset(in_files)
            the_pool.close()

            out_file = run_dir / ('all_' + stage + '_ensemble' + domain + '.nc')
            ds.to_netcdf(out_file)

    #-------------------------------------------------------
    # Wrape it up.
    t_end = time.time()
    print("Wrote collected output to : ", out_file)
    print('DART data collection took: %2.4f sec' % (t_end-t_start))

    return(True)


if __name__=='__main__':

    parser = argparse.ArgumentParser(
        description='Script interface to wrfhydropy.ioutils.open_nwm_dataset.'
    )
    requiredNamed = parser.add_argument_group('Required, named arguments')
    requiredNamed.add_argument(
        '--run_dir',
        required=True,
        action='store',
        help='The path to the dart run_dir.'
    )
    requiredNamed.add_argument(
        '--n_cores',
        required=True,
        action='store',
        help='The number of cores to use in the processing.'
    )
    # Not sure these two will ever be used, but they are easy to leave until later.
    parser.add_argument(
        '--spatial_indices',
        action='store',
        default=None,
        help='A comma separated list of spatial/feature_id indices (quoted in the shell) or "None".'
    )
    parser.add_argument(
        '--drop_variables',
        action='store',
        default=None,
        help='A comma separated list of spatial/feature_id indices (quoted in the shell) or "None".'
    )
    args = parser.parse_args()

    run_dir = args.run_dir
    n_cores = int(args.n_cores)
    spatial_indices = args.spatial_indices
    drop_variables = args.drop_variables

    if spatial_indices == 'None':
        spatial_indices = None
    else:
        spatial_indices = [int(ind) for ind in spatial_indices.replace(" ", "").split(',')]

    if drop_variables == 'None':
        drop_variables = None
    else:
        drop_variables = [var for var in drop_variables.replace(" ", "").split(',')]

    success = collect_run_output(
        run_dir,
        n_cores,
        spatial_indices=spatial_indices,
        drop_variables=drop_variables
    )

    # Shell translation.
    sys.exit(int(not success))
