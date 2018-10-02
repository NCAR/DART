import datetime
import f90nml
import pathlib
import pickle
import shlex
import shutil
import subprocess
import time
import wrfhydropy

from .gregorian import gregorian
from .setup_experiment_tools import replace_in_file, get_top_level_dir_from_config

def setup_usgs_daily(config):

    print('Preparing USGS daily observations.')

    usgs_daily_config = config['observation_preparation']['USGS_daily']

    input_dir = usgs_daily_config['input_dir']
    output_dir = usgs_daily_config['output_dir']
    # Output directory: make if DNE
    output_dir.mkdir(exist_ok=False, parents=True)

    # stage the scripts
    make_daily = output_dir / 'makedaily.csh'
    parallel_daily_batch = output_dir / 'parallel_daily.batch.csh'
    
    dart_repo_usgs_scripts_path = \
        config['dart']['dart_src'] / 'observations/obs_converters/USGS/scripts'
    
    _ = shutil.copy(dart_repo_usgs_scripts_path / 'makedaily.TEMPLATE.csh', make_daily)
    _ = shutil.copy(
        dart_repo_usgs_scripts_path / 'parallel_daily.batch.TEMPLATE.csh',
        parallel_daily_batch
    )

    # identity or regular obs converter?
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
    replace_in_file(make_daily, 'CONVPROG_TEMPLATE', str(obs_conv_prog.name))

    replace_in_file(make_daily, 'DART_DIR', str(dart_build_dir))

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

    replace_in_file(
        make_daily,
        'GAGES_LIST_FILE_TEMPLATE',
        input_nml[converter_nml]['gages_list_file']
    )

    # input.nml location_file/Routelink: link file, edit input.nml
    run_dir = config['experiment']['run_dir']
    m0 = pickle.load(open(run_dir / "member_000/WrfHydroSim.pkl", 'rb'))
    route_link_f = run_dir / 'member_000' / m0.base_hydro_namelist['hydro_nlist']['route_link_f']
    (output_dir / route_link_f.name).symlink_to(route_link_f)
    input_nml[converter_nml]['location_file'] = route_link_f.name
    replace_in_file(make_daily, 'ROUTELINK_TEMPLATE', route_link_f.name)

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
        (output_dir / hydro_rst_file.name).symlink_to(hydro_rst_file)
        input_nml['model_nml']['domain_order'] = 'hydro'
        input_nml['model_nml']['domain_shapefiles'] = str(hydro_rst_file.name)
        replace_in_file(
            make_daily,
            'HYDRO_RST_TEMPLATE',
            'ln -s $origindir/' + hydro_rst_file.name + ' .'
        )

        f90nml.Namelist(m0.base_hydro_namelist).write(output_dir / 'hydro.namelist', force=True)
        top_level_dir = get_top_level_dir_from_config(config, m0)
        (output_dir / top_level_dir).symlink_to(config['wrf_hydro']['domain_src'] / top_level_dir)
        replace_in_file(
            make_daily,
            'TOP_LEVEL_CONFIG_TEMPLATE',
            'ln -s $origindir/' + top_level_dir + ' .'
        )

    else:

        replace_in_file(make_daily, 'HYDRO_RST_TEMPLATE', '')
        replace_in_file(make_daily, 'TOP_LEVEL_CONFIG_TEMPLATE', '')

    # Now we are done editing it, write the input.nml back out.
    input_nml.write(output_dir / 'input.nml')

    # Create the command files here
    # Create all the commands then write to files.
    current_time = in_start_time
    all_cmds=[]
    while current_time <= in_end_time:
        cmd_str = "csh ./makedaily.csh {0} " + str(input_dir)
        cmd_str = cmd_str.format(current_time.strftime('%Y%m%d'))
        all_cmds.append(cmd_str)
        current_time = current_time + datetime.timedelta(days=1)
        
    # The length of all_cmds must be zero mod ppn, if it is not, pad it with dummy commands.
    ppn = usgs_daily_config['scheduler']['ppn']
    remainder = len(all_cmds) % ppn
    n_pad = ppn - remainder
    for pp in range(n_pad):
        all_cmds.append("echo 'This command is just padding this last cmdfile.'")

    n_cmd_files = len(all_cmds) // ppn

    for nn in range(n_cmd_files): 
        cmd_filename = 'day_cmdfile.' + str(nn+1)
        with open(output_dir / cmd_filename, 'w') as opened_file:
            for cc in range(ppn*nn, ppn*(nn+1)):
                _ = opened_file.write(all_cmds[cc] + '\n')

    # Set the PBS directives (cheyenne)

    # PBS options from config
    # Short-hand
    usgs_sched = usgs_daily_config['scheduler']

    # The easy ones.
    replace_in_file(parallel_daily_batch, 'JOB_NAME_TEMPLATE', usgs_sched['job_name'])
    replace_in_file(parallel_daily_batch, 'ACCOUNT_TEMPLATE', usgs_sched['account'])
    replace_in_file(parallel_daily_batch, 'EMAIL_WHO_TEMPLATE', usgs_sched['email_who'])
    replace_in_file(parallel_daily_batch, 'EMAIL_WHEN_TEMPLATE', usgs_sched['email_when'])
    replace_in_file(parallel_daily_batch, 'QUEUE_TEMPLATE', usgs_sched['queue'])

    # Wall time
    usgs_walltime = usgs_sched['walltime']
    if len(usgs_walltime.split(':')) == 2:
        usgs_walltime = usgs_walltime + ':00'
    usgs_walltime = 'walltime=' + usgs_walltime
    replace_in_file(parallel_daily_batch, 'WALLTIME_TEMPLATE', usgs_walltime)

    wait_file = output_dir / '.parallel_daily_batch_not_complete'
    proc = subprocess.Popen(
        shlex.split('touch ' + wait_file.name),
        cwd=output_dir
    )
    proc.wait()
   
    proc = subprocess.Popen(
        shlex.split('qsub ' + parallel_daily_batch.name),
        cwd=output_dir
    )

    while wait_file.exists():
        print('setup_usgs_daily: \nWaiting for removal of wait file: ' + str(wait_file))
        time.sleep(20)

    # Link the obs_seq files to the "all_obs_dir" for the experiment.    
    all_obs_dir = pathlib.PosixPath(config['observation_preparation']['all_obs_dir'])
    all_obs_seq = output_dir.glob('obs_seq.*')
    for oo in all_obs_seq:
        (all_obs_dir / oo.name).symlink_to(oo)
        
        
    return True
