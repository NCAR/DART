import os
import pathlib
import pytest
import hydrodartpy
import hydrodartpy.core.setup_experiment_tools as hdp_tools
import xarray as xr
import yaml

def calc_file_stats(file):
    ds = xr.open_dataset(file)
    data_vars = set(ds.variables).difference(set(ds.coords))
    stats_dict = {}
    for vv in data_vars:
        # dtype check?
        # if ds[vv].dtype
        stats_dict[vv] = {
            'mean': ds[vv].mean().values.tolist(),
            'std': ds[vv].std().values.tolist()}
    return stats_dict


def test_setup_experiment(config_file):
    setup_return_code = hydrodartpy.setup_experiment(config_file=config_file)
    assert setup_return_code == 0


def test_run_experiment(config_dict):
    run_dir = pathlib.Path(config_dict['experiment']['run_dir'])
    os.chdir(run_dir)
    run_return_code = hydrodartpy.run_filter_experiment(run_dir)
    assert run_return_code == 0


def test_experiment_results(config_dict, answer_file):
    run_dir = config_dict['experiment']['run_dir']

    # generate answer outputs
    # Note that this depends very much on the dependencies
    # of wrfhydropy among other factors.
    # check_files = (
    #     run_dir / ('output/2018081503')).glob(
    #     'analysis_member_00*.2018081503.nc')
    # answer = {
    #     str(ff.relative_to(run_dir)): calc_file_stats(ff) for
    #     ff in check_files }
    # with open(answer_file, 'w') as out_file:
    #     _ = yaml.dump(answer, out_file)
    # answer_orig = answer
    # end generate answer outputs

    with open(run_dir / answer_file) as in_file:
        answer = yaml.load(in_file, Loader=yaml.FullLoader)

    check_files = list(answer.keys())
    candidates = {
        str(ff): calc_file_stats(run_dir / ff) for
        ff in check_files }

    assert candidates == answer
