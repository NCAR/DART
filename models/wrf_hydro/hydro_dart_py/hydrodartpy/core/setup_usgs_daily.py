import copy
import datetime
import f90nml
import pathlib
import pickle
import shlex
import shutil
import subprocess
import time
import yaml

from .gregorian import gregorian
from .setup_experiment_tools import replace_in_file, get_top_level_dir_from_config
from .create_usgs_daily_obs_seq import create_usgs_daily_obs_seq

def setup_usgs_daily(
    config,
    config_file: None
):

    print('Preparing USGS daily observations.')

    # Setup the namelist and establish the "outputdir". The create_usgs_daily_obs_seq
    # actually creates the obs sequences.

    usgs_daily_config = config['observation_preparation']['USGS_daily']
    input_dir = usgs_daily_config['input_dir']
    output_dir = usgs_daily_config['output_dir']
    # Output directory: make if DNE
    #output_dir.mkdir(exist_ok=False, parents=True)
    output_dir.mkdir(exist_ok=True, parents=True)

    # converter: identity or regular obs converter?
    # Check that the desired obs converter is in the dart build
    exp_dir = config['experiment']['experiment_dir']
    dart_build_dir = config['dart']['build_dir']
    dart_compile = pickle.load(open(exp_dir / dart_build_dir / 'DartCompile.pkl', 'rb'))

    if usgs_daily_config['identity_obs']:
        ocp = dart_compile.models__wrf_hydro__work.exes['create_identity_streamflow_obs']
    else:
        ocp = dart_compile.observations__obs_converters__USGS__work.exes['convert_streamflow']

    obs_conv_prog = ocp
    _ = shutil.copy(obs_conv_prog, output_dir / obs_conv_prog.name)

    # input.nml: patch.
    converter_nml = str(obs_conv_prog.name) + '_nml'
    obs_conv_patches = usgs_daily_config['input_nml_patches'][converter_nml]
    input_nml = f90nml.read(obs_conv_prog.parent / 'input.nml')
    internal_patches = ['input_files', 'location_file']
    special_patches = ['gages_file_list']
    for kk in usgs_daily_config['input_nml_patches'].keys():
        if kk in internal_patches + special_patches:
            if kk in internal_patches:
                warnings.warn("USGS observation converter namelist patch is applied internally: " + kk)
            pass
        input_nml[kk].update(usgs_daily_config['input_nml_patches'][kk])

    # input.nml gage_file_list: Allow a file or a list: link or construct file, set in input.nml. 
    wanted_gages = usgs_daily_config['wanted_gages']
    if type(wanted_gages) in [str, pathlib.PosixPath]:
        wanted_gages = pathlib.PosixPath(wanted_gages)
        (output_dir / wanted_gages.name).symlink_to(wanted_gages)
        input_nml[convert_nml]['gages_list_file'] = wanted_gages.name
    elif type(wanted_gages) is list:
        default_filename = 'wanted_gages_list.txt'
        input_nml[converter_nml]['gages_list_file'] = default_filename
        with open(output_dir / default_filename, 'w') as opened_file:
            for gg in wanted_gages:
                _ = opened_file.write(str(gg) + '\n')
    else:
        raise ValueError("wanted_gages must be either string or list type.")

    # input.nml location_file/Routelink: link file, edit input.nml
    run_dir = config['experiment']['run_dir']
    m0 = pickle.load(open(run_dir / "member_000/WrfHydroSim.pkl", 'rb'))
    route_link_f = run_dir / 'member_000' / m0.base_hydro_namelist['hydro_nlist']['route_link_f']
    if not route_link_f.is_file():
        (output_dir / route_link_f.name).symlink_to(route_link_f)
    input_nml[converter_nml]['location_file'] = route_link_f.name

    #input.nml input_files: create a list of files in the start and end range.
    in_start_time = datetime.datetime.strptime(str(usgs_daily_config['start_date']), '%Y-%m-%d')
    in_end_time = datetime.datetime.strptime(str(usgs_daily_config['end_date']), '%Y-%m-%d')
    all_input_files = sorted(input_dir.glob("*.usgsTimeSlice.ncdf"))
    input_files_requested = []
    for ff in all_input_files:
        file_time = datetime.datetime.strptime(ff.name.split('.')[0], '%Y-%m-%d_%H:%M:%S')
        # For end_time, add in a day and use a stricly less than... 
        if file_time >= in_start_time and file_time < (in_end_time + datetime.timedelta(days=1)):
            input_files_requested.append(ff)

    default_filename = 'list_of_obs_files.txt'
    input_nml[converter_nml]['input_files'] = default_filename
    with open(output_dir / default_filename, 'w') as opened_file:
        for ff in input_files_requested:
            _ = opened_file.write(str(ff) + '\n')

    # For identity obs, the model_mod and the model_nml namelist are used....
    # This requires hydro_rst file, hydro.namelist, attendant file.... 
    # Parameter and LSM files should not be needed, so just set the hydro restart file.
    if usgs_daily_config['identity_obs']:

        hydro_rst_file = run_dir / 'member_000' / m0.base_hydro_namelist['hydro_nlist']['restart_file']
        if not hydro_rst_file.is_file():
            (output_dir / hydro_rst_file.name).symlink_to(hydro_rst_file)
        input_nml['model_nml']['domain_order'] = 'hydro'
        input_nml['model_nml']['domain_shapefiles'] = str(hydro_rst_file.name)

        f90nml.Namelist(m0.base_hydro_namelist).write(output_dir / 'hydro.namelist', force=True)
        top_level_dir = get_top_level_dir_from_config(config, m0)
        nwm_dir = config['wrf_hydro']['domain_src'] / top_level_dir
        if not nwm_dir.is_dir():
            (output_dir / top_level_dir).symlink_to(config['wrf_hydro']['domain_src'] / top_level_dir)

    # Now we are done editing it, write the input.nml back out.
    out_input = output_dir / 'input.nml'
    if not out_input.is_file():
        input_nml.write(output_dir / 'input.nml')

    # Symlink the config file into the output_dir so the default yaml file name
    # can be used by create_usgs_daily_obs_seq.
    if config_file is None:
        config_file = sorted(exp_dir.glob('original.*.yaml'))[0]
    if not config_file.is_file():
        (output_dir / 'config_file.yaml').symlink_to(config_file)

    # Stage the file that does the batch processing.
    this_file = pathlib.Path(__file__)
    batcher_base = 'create_usgs_daily_obs_seq.py'
    pyscript = this_file.parent / batcher_base
    if not pyscript.is_file():
        (output_dir / batcher_base).symlink_to(this_file.parent / batcher_base)

    # Setup the scheduled script.
    orig_submit_script = this_file.parent / 'submission_scripts/submit_usgs_daily_obs_converter.sh'
    this_submit_script = output_dir / 'submit_usgs_daily_obs_converter.sh'
    shutil.copy(orig_submit_script, this_submit_script)

    # Set the PBS directives (cheyenne)

    # PBS options from config
    # Short-hand
    usgs_sched = config['observation_preparation']['USGS_daily']['scheduler']

    if usgs_sched is not None and usgs_sched != 'None':
    
        # The easy ones.
        replace_in_file(this_submit_script, 'JOB_NAME_TEMPLATE', usgs_sched['job_name'])
        replace_in_file(this_submit_script, 'ACCOUNT_TEMPLATE', usgs_sched['account'])
        replace_in_file(this_submit_script, 'EMAIL_WHO_TEMPLATE', usgs_sched['email_who'])
        replace_in_file(this_submit_script, 'EMAIL_WHEN_TEMPLATE', usgs_sched['email_when'])
        replace_in_file(this_submit_script, 'QUEUE_TEMPLATE', usgs_sched['queue'])

        # Wall time
        usgs_walltime = usgs_sched['walltime']
        if len(usgs_walltime.split(':')) == 2:
            usgs_walltime = usgs_walltime + ':00'
        usgs_walltime = 'walltime=' + usgs_walltime
        replace_in_file(this_submit_script, 'WALLTIME_TEMPLATE', usgs_walltime)

        # Select statement
        # Right now, only single node processing
        select_stmt = 'select=1:ncpus={ncpus}:mpiprocs={mpiprocs}:mem={reqmem}GB'.format(
            **{
                'ncpus': usgs_sched['ncpus'],
                'mpiprocs': usgs_sched['mpiprocs'],
                'reqmem': usgs_sched['reqmem']
            }
        )
        replace_in_file(this_submit_script, 'PBS_SELECT_TEMPLATE', select_stmt)

        wait_file = output_dir / '.this_submit_script_not_complete'
        replace_in_file(this_submit_script, 'WAIT_FILE_TEMPLATE', str(wait_file))

        proc = subprocess.Popen(
            shlex.split('touch ' + wait_file.name),
            cwd=output_dir
        )
        proc.wait()

        proc = subprocess.Popen(
            shlex.split('qsub ' + this_submit_script.name),
            cwd=output_dir
        )
        proc.wait()

        print('Job submitted. \nWait file: ' + str(wait_file) + ' ...')
        while wait_file.exists():
            msg = 'Last check for wait file {twirl}: ' + \
                  datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            for twirl in ['|', '/', '-', '\\', '|', '/', '-', '\\']:
                print(msg.format(**{'twirl':twirl}), end='\r')
                time.sleep(10/8)

    else:

        result = create_usgs_daily_obs_seq(config)

    # Link the obs_seq files to the "all_obs_dir" for the experiment.    
    all_obs_dir = pathlib.PosixPath(config['observation_preparation']['all_obs_dir'])
    all_obs_seq = output_dir.glob('obs_seq.*')
    for oo in all_obs_seq:
        if not oo.is_file():
            (all_obs_dir / oo.name).symlink_to(oo)

    return 0
