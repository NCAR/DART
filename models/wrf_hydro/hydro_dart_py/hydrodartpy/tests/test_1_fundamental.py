import os
import pathlib
import pytest
import hydrodartpy
import hydrodartpy.core.setup_experiment_tools as hdp_tools


def test_setup_experiment(config_file):
    setup_return_code = hydrodartpy.setup_experiment(config_file=config_file)
    assert setup_return_code == 0


def test_run_experiment(config_dict):
    run_dir = pathlib.Path(config_dict['experiment']['run_dir'])
    os.chdir(run_dir)
    run_return_code = hydrodartpy.run_filter_experiment(run_dir)
    assert run_return_code == 0

#def test_experiment_results(config_dict, answer_file):
# check output/

