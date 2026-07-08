#!/usr/bin/env python

import argparse
import dask
import dask.bag as db
from multiprocessing.pool import Pool
import pathlib
import sys
import time
from typing import Union
import xarray as xr
import itertools
import collections


def is_not_none(x):
    return x is not None

def group_member(ds: xr.Dataset) -> int:
    return ds.member.item(0)

def merge_time(ds_list: list) -> xr.Dataset:
    return xr.concat(ds_list, dim='time', coords='minimal')

def merge_member(ds_list: list) -> xr.Dataset:
    return xr.concat(ds_list, dim='member', coords='minimal')


def preprocess_dart_data(
    path,
    chunks: dict = None,
    spatial_indices: list = None,
    drop_variables: list = None
) -> xr.Dataset:

    # This non-optional is different from preprocess_nwm_data
    # I kinda dont think this should be optional for dart experiment/run collection.
    # try:
    ds = xr.open_dataset(path)
    # except OSError:
    #    print("Skipping file, unable to open: ", path)
    #    return None

    # May need to add time... do this before changing any dimensions.
    for key in ds.variables.keys():
        if 'time' not in ds[key].dims:
            ds[key] = ds[key].expand_dims('time')

    if drop_variables is not None:
        to_drop = set(ds.variables).intersection(set(drop_variables))
        if to_drop != set():
            ds = ds.drop_vars(to_drop)

    # This member definition is different from preprocess_nwm_data
    member = int(ds.attrs['DART_file_information'].split()[-1])
    ds.coords['member'] = member

    # Spatial subsetting
    if spatial_indices is not None:
        ds = ds.isel(feature_id=spatial_indices)

    # Chunk here?

    return ds


def open_dart_dataset(
    paths: list,
    domain: list = None,
    chunks: dict = None,
    spatial_indices: list = None,
    drop_variables: list = None,
    npartitions: int = None,
    attrs_keep: list = None
) -> xr.Dataset:
    """Open a multi-file ensemble pywatershed output dataset
    Args:
paths: List ,iterable, or generator of file paths to pywatershed netcdf output files
        chunks: chunks argument passed on to xarray DataFrame.chunk() method
        preprocess_member: A function that identifies the member from the file or filename.
        attrs_keep: A list of the global attributes to be retained.
    Returns:
        An xarray dataset of dask arrays chunked by chunk_size along the feature_id
        dimension concatenated along the time and member dimensions.
    """

    # Explanation:
    # Xarray currently first requires concatenation along existing dimensions (e.g. time)
    # over the individual member groups, then it allows concatenation along the member
    # dimensions.

    # Set partitions
    # This is arbitrary
    if npartitions is None:
        npartitions = dask.config.get('pool')._processes * 4

    #print(paths)
    paths_bag = db.from_sequence(paths)
    print(paths_bag)

    ds_list = paths_bag.map(
        preprocess_dart_data,
        chunks=chunks,
        spatial_indices=spatial_indices,
        drop_variables=drop_variables
    ).filter(is_not_none).compute()

    the_sort = sorted(ds_list, key=group_member)
    ds_groups = [list(it) for k, it in itertools.groupby(the_sort, group_member)]
    group_bag = db.from_sequence(ds_groups)  # , npartitions=npartitions)
    #print(ds_groups)
    #print(group_bag)
    ds_list = group_bag.map(merge_time).compute()
    del group_bag, ds_groups, the_sort
    dart_dataset = merge_member(ds_list)
    del ds_list

    # Xarray sets nan as the fill value when there is none. Dont allow that...
    for key, val in dart_dataset.variables.items():
        if '_FillValue' not in dart_dataset[key].encoding:
            dart_dataset[key].encoding.update({'_FillValue': None})

    # Clean up attributes
    new_attrs = collections.OrderedDict()
    if attrs_keep is not None:
        for key, value in dart_dataset.attrs.items():
            if key in attrs_keep:
                new_attrs[key] = dart_dataset.attrs[key]

    dart_dataset.attrs = new_attrs

    # The existing DART convention.
    if domain ==  '_d02':
       dart_dataset = dart_dataset.transpose('time', 'member', 'nhru')
    else:
       dart_dataset = dart_dataset.transpose('time', 'member', 'nsegment')
    
    

    # Break into chunked dask array
    if chunks is not None:
        dart_dataset = dart_dataset.chunk(chunks=chunks)

    return dart_dataset

def collect_run_output(
    run_dir: Union[pathlib.Path,str],
    n_cores: int,
    spatial_indices: list=None,
    drop_variables: list=None
):

    run_dir = pathlib.Path(run_dir)

    t_start = time.time()

    # The DART_cleanup_add_time.csh step can be done in the preprocess steps below.

    # -------------------------------------------------------
    # 1. These files only concatenate along the time dimension.
    #    We can use xr.open_mfdataset out of the box!
    stage_list = ['input', 'preassim', 'analysis', 'output']

    type_base_list = ['mean', 'sd', 'priorinf_mean', 'priorinf_sd', 'postinf_mean', 'postinf_sd', 'hybridweight_mean', 'hybridweight_sd']
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
    #    The explicit handling of individual members happens in open_dart_dataset
    stage_list = ['input', 'preassim', 'analysis', 'output']
    domain_list = ['_d01', '_d02']

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
                #print(in_files)
                ds = open_dart_dataset(in_files, domain)
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
        description='Script to collect all the pywatershed outputs.'
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
