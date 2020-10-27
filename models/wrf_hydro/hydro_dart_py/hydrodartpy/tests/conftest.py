import pathlib
import shutil

import pytest
import hydrodartpy


def pytest_addoption(parser):

    # Required args:
    parser.addoption(
        '--domain_dir',
        required=True,
        action='store',
        help='domain directory'
    )

    parser.addoption(
        '--output_dir',
        required=True,
        action='store',
        help='test output directory'
    )

    parser.addoption(
        '--candidate_dir',
        required=True,
        action='store',
        help='candidate model directory'
    )

    parser.addoption(
        '--reference_dir',
        required=True,
        action='store',
        help='reference model directory'
    )

    parser.addoption(
        "--config",
        required=True,
        action='store',
        help=("The configuration to test, "
              "must be one listed in trunk/NDHMS/hydro_namelist.json keys.")
    )

    # Optional args:
    parser.addoption(
        '--compiler',
        default='gfort',
        required=False,
        action='store',
        help='compiler, options are ifort or gfort'
    )

    parser.addoption(
        "--option_suite",
        required=False,
        action='store',
        help=("An option suite to test on top of the specified configuration,"
              "must be one listed in hydro_option_suites.json")
    )

    parser.addoption(
        '--ncores',
        default='2',
        required=False,
        action='store',
        help='Number of cores to use for testing'
    )

    parser.addoption(
        '--scheduler',
        action='store_true',
        help='Use PBS scheduler on cheyenne'
    )

    parser.addoption(
        '--nnodes',
        default='2',
        required=False,
        help='Number of nodes to use for testing if running on scheduler'
    )

    parser.addoption(
        '--account',
        default='NRAL0017',
        required=False,
        action='store',
        help='Account number to use if using a scheduler.'
    )

    parser.addoption(
        '--walltime',
        default='02:00:00',
        required=False,
        action='store',
        help='Wall clock time for each test run in hh:mm:ss format'
    )

    parser.addoption(
        '--queue',
        default='regular',
        required=False,
        action='store',
        help='Queue to use if running on NCAR Cheyenne, options are regular, '
        'premium, or shared'
    )


@pytest.fixture(scope="session")
def config(request):

    configuration = request.config.getoption("--config")
    config = hydrodartpy.setup_experiment_tools.establish_config(spec_file=configuration)
    return config


@pytest.fixture(scope="session")
def candidate_channel_only_sim(request):

    domain_dir = request.config.getoption("--domain_dir")
    compiler = request.config.getoption("--compiler")
    candidate_dir = request.config.getoption("--candidate_dir")
    configuration = request.config.getoption("--config")
    option_suite = request.config.getoption("--option_suite")
    ncores = request.config.getoption("--ncores")
    nnodes = request.config.getoption("--nnodes")
    scheduler = request.config.getoption("--scheduler")
    account = request.config.getoption("--account")
    walltime = request.config.getoption("--walltime")
    queue = request.config.getoption("--queue")

    candidate_channel_only_sim = _make_sim(
        domain_dir=domain_dir,
        compiler=compiler,
        source_dir=candidate_dir,
        configuration=configuration,
        option_suite=option_suite,
        ncores=ncores,
        nnodes=nnodes,
        scheduler=scheduler,
        account=account,
        walltime=walltime,
        queue=queue
    )

    # Channel and bucket mode is forc_typ = 10.
    candidate_channel_only_sim.base_hrldas_namelist['wrf_hydro_offline']['forc_typ'] = 10
    return candidate_channel_only_sim


@pytest.fixture(scope="session")
def reference_sim(request):

    domain_dir = request.config.getoption("--domain_dir")
    compiler = request.config.getoption("--compiler")
    reference_dir = request.config.getoption("--reference_dir")
    configuration = request.config.getoption("--config")
    option_suite = request.config.getoption("--option_suite")
    ncores = request.config.getoption("--ncores")
    nnodes = request.config.getoption("--nnodes")
    scheduler = request.config.getoption("--scheduler")
    account = request.config.getoption("--account")
    walltime = request.config.getoption("--walltime")
    queue = request.config.getoption("--queue")

    reference_sim = _make_sim(
        domain_dir=domain_dir,
        compiler=compiler,
        source_dir=reference_dir,
        configuration=configuration,
        option_suite=option_suite,
        ncores=ncores,
        nnodes=nnodes,
        scheduler=scheduler,
        account=account,
        walltime=walltime,
        queue=queue
    )

    return reference_sim


@pytest.fixture(scope="session")
def output_dir(request):
    configuration = request.config.getoption("--config")
    output_dir = request.config.getoption("--output_dir")

    output_dir = pathlib.Path(output_dir)
    output_dir = output_dir / configuration

    if output_dir.is_dir() is True:
        shutil.rmtree(str(output_dir))

    output_dir.mkdir(parents=True)
    return output_dir


@pytest.fixture(scope="session")
def ncores(request):
    ncores = request.config.getoption("--ncores")

    return ncores
