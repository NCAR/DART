from boltons.iterutils import remap
import hydrodartpy
import os
import pathlib
import pytest
import shlex
import shutil
import subprocess
import yaml

import hydrodartpy.core.setup_experiment_tools as hdp_tools


repo_dir = hdp_tools.repo_dir


def if_e_rm(path: pathlib.Path):
    if path.is_file():
        _ = path.unlink()
    elif path.is_dir():
        _ = shutil.rmtree(str(path))
    else:
        return 1
    return 0


def pytest_addoption(parser):
    parser.addoption(
        "--exp_yaml",
        required=True,
        type=str,
        nargs=1,
        action='store',
        help=("The path to the YAML file for the experiment."))

    parser.addoption(
        "--exp_answer_file",
        type=str,
        nargs=1,
        required=True,
        action='store',
        help=("The directory where expected experiment outputs live."))

    parser.addoption(
        "--scratch_dir",
        type=str,
        nargs=1,
        required=False,
        default=['None'],
        action='store',
        help=("The local scratch where experiment output will go."))

    parser.addoption(
        "--use_existing_build",
        required=False,
        default=False,
        action='store_true',
        help=("Use existing dart and wrfhydro builds?"))
    # Could also use existing obs and/or initial ensemble


@pytest.fixture(scope="session")
def test_dir():
    return pathlib.Path(os.getcwd()).resolve()


@pytest.fixture(scope="session")
def answer_file(request, test_dir):
    exp_answer_file = test_dir / (
        request.config.getoption("--exp_answer_file")[0])
    return exp_answer_file


@pytest.fixture(scope="session")
def config_file(request, test_dir):
    exp_yaml = test_dir / (
        request.config.getoption("--exp_yaml")[0])
    scratch_dir = request.config.getoption("--scratch_dir")[0]
    use_existing_build = request.config.getoption("--use_existing_build")

    # scratch guessing
    if scratch_dir == 'None':
        tmp_dir = os.getenv('TMPDIR')
        if tmp_dir == '':
            raise ValueError('TMPDIR environment variable not found, '
                             'you must specify --scratch_dir')
        tmp_dir = pathlib.Path(tmp_dir)
        if not tmp_dir.exists():
            raise ValueError('The tmp_dir does not exist',
                             str(tmp_dir))
        scratch_dir = tmp_dir

    # read the yaml
    exp_config = hdp_tools.establish_config(exp_yaml)

    ## check the yaml compiler matches the current modules
    module_list = subprocess.run(
        'module list', shell=True, stderr=subprocess.PIPE)
    the_compiler =  exp_config['wrf_hydro']['compiler']
    wrong_compiler_msg = (
        'The current compiler does not appear to match the compiler in the '
        'YAML file: ' + the_compiler+ '. Please quit python and change '
        'your modules or change the YAML:' + exp_config['wrf_hydro']['compiler'])
    if the_compiler == 'ifort':
        assert 'intel' in module_list.stderr.decode('utf-8'), wrong_compiler_msg
    elif the_compiler == 'gfort':
        assert 'intel' in module_list.stderr.decode('utf-8'), wrong_compiler_msg
    else:
        raise ValueError("What compiler are you using???")

    out_dir = scratch_dir / ('hydrodartpy_tests/' + exp_yaml.parent.name)

    # alter the config values
    exp_dir = out_dir / 'experiment'
    exp_config['experiment']['experiment_dir'] = str(exp_dir)
    run_dir = out_dir / 'run'
    exp_config['experiment']['run_dir'] = str(run_dir)

    ini_dir = out_dir / 'initial_ensemble'
    exp_config['initial_ens']['path'] = str(ini_dir)

    all_obs_dir = out_dir / 'obs_seqs'
    exp_config['observation_preparation']['all_obs_dir'] = str(all_obs_dir)
    ## Some semi-hard-coded patches
    if 'USGS_daily' in exp_config['observation_preparation'].keys():
        exp_config['observation_preparation']['USGS_daily']['output_dir'] = str(
            all_obs_dir / 'USGS_daily')

    if ('noise_function_files' in
        exp_config['run_experiment']['perturb_forcing'].keys()):
        pf = exp_config['run_experiment']['perturb_forcing']
        nffs = pf['noise_function_files']
        nffs = [str(repo_dir.joinpath(nff)) for nff in nffs]
        pf['noise_function_files'] = nffs

    test_dir = pathlib.Path(os.getcwd())
    exp_config['dart']['dart_src'] = str(
        test_dir.parent.parent.parent.parent.parent)
    ## could clone/update a model repo here
    wrf_hydro_src = test_dir / 'data/wrf_hydro_nwm_public'
    exp_config['wrf_hydro']['wrf_hydro_src'] = str(wrf_hydro_src)

    print('\n\n'
          '**********************************\n'
          'WRF-Hydro repository information:\n')
    if use_existing_build:
        if not wrf_hydro_src.exists():
            raise ValueError('Repository missing, unfortunately cant '
                             'use existing build for this test.')
    elif not wrf_hydro_src.exists():
            clone_cmd = 'git clone https://github.com/NCAR/wrf_hydro_nwm_public.git'
            clone_result = subprocess.run(
                shlex.split(clone_cmd), cwd=wrf_hydro_src.parent)
        # else:
            # I have Decided against continuous integration.
            #update_cmd = 'git pull origin master'
            #update_result = subprocess.run(
            #    shlex.split(update_cmd), cwd=wrf_hydro_src)

    checkout_cmd = 'git checkout -f nwm-v2.1-beta3'
    sub_result = subprocess.run(shlex.split(checkout_cmd), cwd=wrf_hydro_src)

    commit_cmd = 'git --no-pager log -n 1 --pretty=oneline'
    sub_result = subprocess.run(shlex.split(commit_cmd), cwd=wrf_hydro_src)

    print('\n\n'
          '**********************************\n')

    exp_config['ensemble']['constructor'] = str(exp_yaml.parent / 'constructor.py')

    test_yaml = out_dir / 'test.yaml'

    # Set up test directories
    exp_config['dart']['use_existing_build'] = use_existing_build
    exp_config['wrf_hydro']['use_existing_build'] = use_existing_build
    if not use_existing_build:
        if out_dir.exists():
            _ = shutil.rmtree(str(out_dir))
        _ = out_dir.mkdir(parents=True)
    else:
        _ = if_e_rm(test_yaml)
        _ = if_e_rm(run_dir)
        _ = if_e_rm(ini_dir)
        _ = if_e_rm(all_obs_dir)
        exp_dir_files = exp_dir.glob('*')
        builds_keep = [
            exp_config['dart']['build_dir'],
            exp_config['wrf_hydro']['build_dir']]
        _ = [if_e_rm(path) for path in exp_dir_files
             if path.name not in builds_keep]

    strify = hdp_tools.visit_abs_paths_to_str
    exp_config_out = remap(exp_config, strify)
    with open(test_yaml, 'w') as out_file:
        _ = yaml.dump(exp_config_out, out_file)

    return test_yaml


@pytest.fixture(scope="session")
def config_dict(config_file):
    return hdp_tools.establish_config(config_file)
