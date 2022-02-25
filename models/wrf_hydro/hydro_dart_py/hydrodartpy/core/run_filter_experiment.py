import datetime
import f90nml
import itertools
import hydrodartpy
import math
import multiprocessing
import operator
import os
import pathlib
import shutil
import shlex
import subprocess
import sys
import time
import wrfhydropy

hostname = subprocess.run('hostname', stdout=subprocess.PIPE).stdout.rstrip().decode('utf-8')

major_delim = '\n\n===========================================================================\n'
minor_delim = '\n---------------------------------------------------------------------------\n'

# Date formats used.
datestr_fmts = {
    'separated'  : '%Y-%m-%d_%H:%M',
    'compact_min': "%Y%m%d%H%M",
    'compact_hr' : "%Y%m%d%H",
    'compact_day': "%Y%m%d"
}

# TODO: put in a utility module
def domain_restart_file_pattern(domain_name: str, date: datetime.datetime):
    """Generate a restart file basename from a domain tag and a date."""
    if domain_name == 'hydro':
        return "HYDRO_RST.{0}_DOMAIN1".format(date.strftime(datestr_fmts['separated']))
    elif domain_name == 'lsm':
        return "RESTART.{0}_DOMAIN1".format(date.strftime(datestr_fmts['compact_hr']))
    elif domain_name == 'param':
        return "parameter_restart.{0}.nc".format(date.strftime(datestr_fmts['separated']))
    else:
        raise ValueError("Invalid domain_name, should be one of ['hydro', 'lsm', 'param'].")


# These parallel copy and move routines could be move to a utiltiy module
def parallel_cp_file(arg_dict):
    if 'if_exists' in arg_dict.keys():
        if arg_dict['if_exists'] and not arg_dict['src_file'].exists():
            return
    shutil.copy2(str(arg_dict['src_file']), str(arg_dict['dest_file']))


def parallel_mv_file(arg_dict):
    if 'if_exists' in arg_dict.keys():
        if arg_dict['if_exists'] and not arg_dict['src_file'].exists():
            return
    arg_dict['src_file'].rename(arg_dict['dest_file'])


# TODO: move the stopwatch to a utility module.
# A decorator/closure to check timings.
def stopwatch(the_func):
    def the_closure(*args, **kw):
        ts = time.time()
        result = the_func(*args, **kw)
        te = time.time()
        print('Timing: ' + str(round(te - ts, 2)) + ' seconds: ' + the_func.__name__)
        return result
    return the_closure


@stopwatch
def prep_obs_seq_for_filter(
    run_dir: pathlib.Path,
    obs_dir: pathlib.Path,
    curr_ens_time: datetime.datetime,
    assim_window_hr: int,
    skip_missing_obs_hours: bool,
    skip_missing_obs_days: bool
):

    # This script edits the obs_sequence_tool_nml and establish
    # the files specified in that nml.
    #
    # Major (known) assumptions:
    # 1) the input obs_seq files come in daily chunks and 
    # 2) the obs_seq input files are named as obs_seq.YYYYmmdd

    # The behavior depends on if obs_sequence_tool is in the run_dir.
    # If obs_sequence_tool is present, then daily obs are substed to "window" obs. Otherwise,
    # filter_nml is passed the start and end times for the obs in the daily obs_seq file.
    # In both cases the name of the correct obs_seq is set in filter_nml.

    print(minor_delim)
    print("Prepare the obs_seq for filter.\n")

    # Solve the assimilation time information. 
    timedelta = datetime.timedelta
    assim_half_window_sec = timedelta(seconds=(60*60/2)*assim_window_hr)

    start_window = curr_ens_time - assim_half_window_sec + timedelta(seconds=1)
    end_window = curr_ens_time + assim_half_window_sec

    start_window_greg = start_window - datetime.datetime(1601, 1, 1)
    end_window_greg = end_window - datetime.datetime(1601, 1, 1)

    # The input.nml
    # Read/edit/write the input.nml for both obs_sequence_tool_nml and filter_nml.
    input_nml_file = run_dir / 'input.nml'
    input_nml = f90nml.read(input_nml_file)
    # Shorthands
    obs_seq_nml = input_nml['obs_sequence_tool_nml']
    filter_nml = input_nml['filter_nml']

    # Clean up previous obs_seq files in the run dir
    filename_seq = obs_seq_nml['filename_seq']
    if not isinstance(filename_seq, list): filename_seq = [filename_seq]
    hard_coded_mess = ['obs_seq.previous', 'obs_seq.daily', 'obs_seq.window']
    cleanup_files = filename_seq + [obs_seq_nml['filename_out']] + hard_coded_mess
    for file in cleanup_files:
        rm_file = run_dir / file
        if rm_file.exists():
            rm_file.unlink()

    # Does the assim window straddle two days?
    assim_window_two_days = start_window.day != end_window.day

    # Daily error message.
    missing_obs_day_err_msg = 'Required daily input obs sequence file(s) are missing: '

    # Identify the input obs_seq files for the current day.
    obs_seq_daily_src = obs_dir / ('obs_seq.'+curr_ens_time.strftime(datestr_fmts['compact_day']))
    obs_seq_daily = run_dir.joinpath('obs_seq.daily')
    if obs_seq_daily_src.exists():
        obs_seq_daily.symlink_to(obs_seq_daily_src)
    else:
        if not skip_missing_obs_days:
            raise FileNotFoundError(missing_obs_day_err_msg + str(obs_seq_daily_src))
        else:
            return 24

    obs_seq_nml['filename_seq'] = 'obs_seq.daily'
    
    if assim_window_two_days:
        prev_time = curr_ens_time - datetime.timedelta(days=1)
        obs_seq_previous_src = obs_dir / ('obs_seq.'+prev_time.strftime(datestr_fmts['compact_day']))
        obs_seq_previous = run_dir.joinpath('obs_seq.previous')
        if obs_seq_previous_src.exists():
            obs_seq_previous.symlink_to(obs_seq_previous_src)
            obs_seq_nml['filename_seq'] = ['obs_seq.previous', 'obs_seq.daily']
            
        #else:
            # This else is a bit draconian, this could be neglected. This is only likely
            # to error when the model is initialized at the beginning of an analysis. If
            # the previous hour ran then this wont happen...
            
            # if not skip_missing_obs_days:
            #     raise FileNotFoundError(missing_obs_day_err_msg + str(obs_seq_previous_src))

    # Set the output file name
    # This file will exist if the obs_sequence tool is used.
    obs_seq_window = run_dir / 'obs_seq.window'
    obs_seq_nml['filename_out'] = obs_seq_window.name

    if (run_dir / 'obs_sequence_tool').exists():
        print('Using obs_sequence_tool.')
        obs_seq_nml['filename_out'] = 'obs_seq.window'
        obs_seq_nml['first_obs_days'] = start_window_greg.days
        obs_seq_nml['first_obs_seconds'] = start_window_greg.seconds
        obs_seq_nml['last_obs_days'] = end_window_greg.days
        obs_seq_nml['last_obs_seconds'] = end_window_greg.seconds

        filter_nml['obs_sequence_in_name'] = 'obs_seq.window'

    else:
        print('Not using obs_sequence_tool.')
        filter_nml['obs_sequence_in_name'] = 'obs_seq.daily'
        filter_nml['first_obs_days'] = start_window_greg.days
        filter_nml['first_obs_seconds'] = start_window_greg.seconds
        filter_nml['last_obs_days'] = end_window_greg.days
        filter_nml['last_obs_seconds'] = end_window_greg.seconds

    # write the input.nml
    input_nml.write(input_nml_file, force=True)

    obs_seq_tool = run_dir.joinpath("obs_sequence_tool")
    if obs_seq_tool.exists():
        # run obs_sequence_tool
        proc = subprocess.run('./obs_sequence_tool', cwd=run_dir)
        
        if not obs_seq_window.exists():
            if not skip_missing_obs_hours:
                missing_obs_hours_err_msg = \
                    'Required hourly data in obs sequence file(s) are missing: '
                raise FileNotFoundError(missing_obs_hours_err_msg + str(obs_seq_daily_src))

            else:
                # This means advance mode by advance_model_hours
                return 1

    return 0

@stopwatch
def parameter_time_advance(
    run_dir: pathlib.Path,
    prev_ens_time: datetime.datetime,
    curr_ens_time: datetime.datetime,
    ens_size: int,
    ncores: int=1
):

    # Parameter "advance". Take the posterior parameters and make them priors.
    # This could involve a model one day. For now, "cp" is the model.
    def mk_param_file_names(the_member: int, the_date: datetime.datetime):
        var_str = 'member_{0}/parameter_restart.{1}.nc'
        the_date_fmt = the_date.strftime(datestr_fmts['separated'])
        return run_dir.joinpath(var_str.format("%03d" % the_member, the_date_fmt))

    param_rst_prev_list = [mk_param_file_names(ee, prev_ens_time) for ee in range(ens_size)]
    param_rst_curr_list = [mk_param_file_names(ee, curr_ens_time) for ee in range(ens_size)]

    if ncores > 1:
        with multiprocessing.Pool(ncores) as pool:
            _ = pool.map(
                parallel_cp_file,
                ({'src_file': src, 'dest_file': dest, 'if_exists': True}
                 for src,dest in zip(param_rst_prev_list, param_rst_curr_list))
            )

    else:
            _ = [parallel_cp_file({'src_file': src, 'dest_file': dest, 'if_exists': True}) 
                 for src,dest in zip(param_rst_prev_list, param_rst_curr_list)]


@stopwatch
def prep_run_dir_for_filter(
    run_dir: pathlib.Path,
    curr_ens_time: datetime.datetime,
    ens_size: int,
    use_hydro_rst: bool,
    use_lsm_rst: bool,
    use_param_rst: bool,
    ncores: int=1
):

    # These are small files, should be no need to do this in parallel.

    def rm_rundir_file(run_dir: pathlib.Path, basename_str: str):
        the_file = run_dir.joinpath(basename_str)
        if the_file.exists():
            the_file.unlink()

    rm_file_list = ['dart_log.out', 'dart_log.nml', 'hydro_file_list.txt',
                    'lsm_file_list.txt', 'param_file_list.txt']
    for ff in rm_file_list:
        rm_rundir_file(run_dir, ff)

    # -------------------------------------------------------
    # Build the individual list files for assimilation.
    # The name of the list file (key) could be taken from the input_nml instead of hard coded.
    domain_list_file_data = {
        'hydro_file_list.txt': {'in_use': use_hydro_rst, 'domain_name': 'hydro'},
        'lsm_file_list.txt': {'in_use': use_lsm_rst, 'domain_name': 'lsm'},
        'param_file_list.txt': {'in_use': use_param_rst, 'domain_name': 'param'}
    }

    def add_restart_to_domain_list(
        run_dir: pathlib.Path,
        curr_ens_time: datetime.datetime,
        ens_size: int,
        list_file_name: str,
        domain_list_file_data: dict
    ):
        if not domain_list_file_data['in_use']:
            return

        file_basename = domain_restart_file_pattern(
            domain_name=domain_list_file_data['domain_name'],
            date=curr_ens_time
        )
        rst_list = sorted(run_dir.glob('member_*/' + file_basename))
        if len(rst_list) != ens_size:
            raise FileExistsError("The number of existing " + domain_list_file_data['domain_name'] +
                                  " restart files does not match the ensemble size.")

        with open(run_dir.joinpath(list_file_name), 'w') as opened_file:
            for mm in rst_list:
                _ = opened_file.write(str(mm) + '\n')

    for key in domain_list_file_data.keys():
        add_restart_to_domain_list(
            run_dir=run_dir,
            curr_ens_time=curr_ens_time,
            ens_size=ens_size,
            list_file_name=key,
            domain_list_file_data=domain_list_file_data[key]
        )


@stopwatch
def run_filter(
    run_dir: pathlib.Path,
    nproc: int,
    cmd: str,
):

    # TODO (JLM): right now hostname is just the master node.
    cmd_to_run = cmd.format(
        **{'hostname':hostname,
           'nproc': nproc,
           'cmd': './filter'
           }
    )
    proc = subprocess.run(shlex.split(cmd_to_run), cwd=str(run_dir))
    if proc.returncode != 0:
        raise ValueError('Filter did not return 0')


@stopwatch
def manage_filter_output(
    run_dir: pathlib.Path,
    curr_ens_time: datetime.datetime,
    n_domains: int,
    ncores: int=1
):

    datestr = curr_ens_time.strftime(datestr_fmts['compact_hr'])

    output_dir_date = run_dir.joinpath('output/' + datestr)
    output_dir_date.mkdir(parents=True, exist_ok=False)

    # Manage inflation files.
    # This is currently not parallel, but could be if the copy and move are separated.
    if n_domains == 1:
        domain_tags = ['']
    else:
        domain_tags = ['_d01', '_d02', '_d03']

    for domain_tag in domain_tags:
        prior_mean_file = run_dir.joinpath('output_priorinf_mean{0}.nc'.format(domain_tag))
        prior_sd_file = run_dir.joinpath('output_priorinf_sd{0}.nc'.format(domain_tag))
        if prior_mean_file.exists():
            shutil.copy2(
                str(prior_mean_file),
                str(output_dir_date / ('output_priorinf_mean{0}.{1}.nc'.format(domain_tag, datestr)))
            )
            shutil.copy2(
                str(prior_sd_file),
                str(output_dir_date / ('output_priorinf_sd{0}.{1}.nc'.format(domain_tag, datestr)))
            )
            prior_mean_file.rename('input_priorinf_mean{0}.nc'.format(domain_tag))
            prior_sd_file.rename('input_priorinf_sd{0}.nc'.format(domain_tag))

	# MEG: Posterior inflation section 
        post_mean_file = run_dir.joinpath('output_postinf_mean{0}.nc'.format(domain_tag))
        post_sd_file = run_dir.joinpath('output_postinf_sd{0}.nc'.format(domain_tag))
        if post_mean_file.exists():
            shutil.copy2(
                str(post_mean_file),
                str(output_dir_date / ('output_postinf_mean{0}.{1}.nc'.format(domain_tag, datestr)))
            )
            shutil.copy2(
                str(post_sd_file),
                str(output_dir_date / ('output_postinf_sd{0}.{1}.nc'.format(domain_tag, datestr)))
            )
            post_mean_file.rename('input_postinf_mean{0}.nc'.format(domain_tag))
            post_sd_file.rename('input_postinf_sd{0}.nc'.format(domain_tag))

    # Other output files.
    potential_output_file_globs = [
        'forecast*mean*nc',   'forecast*sd*nc',  'forecast_member_*.nc',
        'preassim*mean*nc',   'preassim*sd*nc',  'preassim_member_*.nc',
        'postassim*mean*nc',  'postassim*sd*nc', 'postassim_member_*.nc',
        'analysis*mean*nc',   'analysis*sd*nc',  'analysis_member_*.nc',
        'output*mean*nc',     'output*sd*nc'
    ]
    # Combine all of these into a single generator
    full_generator = run_dir.glob(potential_output_file_globs[0])
    for glob in potential_output_file_globs[1:]:
        full_generator = itertools.chain(full_generator, run_dir.glob(glob))

    # "Instantiate" the generator and edit the nameist
    src_output_files = sorted(full_generator)

    # TODO: Move this to utilities.
    def insert_before_extension(the_file: pathlib.Path, the_insert: str):
        ext_split = os.path.splitext(str(the_file))
        return pathlib.Path(ext_split[0] + the_insert + ext_split[1])
    date_insert = '.' + datestr
    dest_output_files = [
        output_dir_date.joinpath(insert_before_extension(file.name, date_insert))
        for file in src_output_files
    ]

    # Farm out the mv to multiple processors.
    if ncores > 1:
        with multiprocessing.Pool(ncores) as pool:
            _ = pool.map(
                parallel_mv_file,
                ({'src_file': src, 'dest_file': dest, 'if_exists': True}
                 for src,dest in zip(src_output_files, dest_output_files))
            )

    else:
            _ = [parallel_mv_file({'src_file': src, 'dest_file': dest, 'if_exists': True}) 
                 for src,dest in zip(src_output_files, dest_output_files)]

    run_dir.joinpath('obs_seq.final').rename(
        output_dir_date.joinpath('obs_seq.final.' + datestr)
    )


@stopwatch
def advance_ensemble(
    run_dir: pathlib.Path,
    n_concurrent: int=None,
    teams_dict: dict=None,
    advance_n_hours: int=None
):

    hydrodartpy.advance_ensemble(
        run_dir=run_dir,
        n_concurrent=n_concurrent,
        teams_dict=teams_dict,
        advance_n_hours=advance_n_hours
    )


# ############################################################################
def run_filter_experiment(run_dir: pathlib.Path):

    print(major_delim)
    print("Starting run_filter_experiment(.py)")

    # Establish the configuration
    config_file = sorted(
        run_dir.joinpath('experiment_dir').glob('original.*.yaml')
    )[0]
    config = hydrodartpy.establish_config(config_file)

    # make some config variebale local for ease of use
    run_dir = pathlib.Path(config['experiment']['run_dir'])
    obs_dir = pathlib.Path(config['observation_preparation']['all_obs_dir'])
    ini_dir = pathlib.Path(config['initial_ens']['path'])
    cycle_end_time = datetime.datetime.strptime(
        config['run_experiment']['time']['end_time'],
        datestr_fmts['separated']
    )
    assim_window_hr = int(
        config['run_experiment']['time']['assim_window']['assim_window_size_hours']
    )
    model_adv_hr = int(config['run_experiment']['time']['advance_model_hours'])
    skip_missing_obs_hours = config['run_experiment']['time']['skip_missing_obs_hours']
    skip_missing_obs_days = config['run_experiment']['time']['skip_missing_obs_days']
    ens_size = int(config['ensemble']['size'])

    # Report some basics before getting started
    print("Run dir: ", run_dir)
    print("Obs dir: ", obs_dir)
    print("ens_size: ", ens_size)

    # Figure out how the job execution is laid out.
    nnodes = int(config['run_experiment']['job_execution']['scheduler']['nnodes'])
    ppn_use = int(config['run_experiment']['job_execution']['scheduler']['ppn_use'])
    ncores_wrf_hydro = int(config['run_experiment']['job_execution']['wrf_hydro']['nproc'])
    ncores_dart = config['run_experiment']['job_execution']['dart']['nproc']
    if ncores_dart is None:
        ncores_dart = nnodes * ppn_use
    else:
        ncores_dart = int(ncores_dart)
    dart_exe_cmd = config['run_experiment']['job_execution']['dart']['exe_cmd']
    
    # Create teams based on PBS_NODEFILE?
    pbs_nodefile = os.environ.get('PBS_NODEFILE')
    if pbs_nodefile is not None:
        n_concurrent=ppn_use
        pbs_nodes = []
        with open(pbs_nodefile, 'r') as infile:
            for line in infile:
                pbs_nodes.append(line.rstrip())

        n_total_processors = len(pbs_nodes) # less may be used.
        n_teams = min(math.floor(len(pbs_nodes) / ncores_wrf_hydro), ens_size)
        teams_dict = {}

        # Map the members on to the teams (this seems overly complicated, should prob
        # consider using pandas:
        member_dirs = ["member_%03d" % the_member for the_member in range(ens_size)]
        member_teams = [the_member % n_teams for the_member in range(ens_size)]
        member_team_seq = [ [dir, team] for dir,team in zip(member_dirs, member_teams)]
        member_team_seq.sort(key = operator.itemgetter(1))
        team_groups = itertools.groupby(member_team_seq, operator.itemgetter(1))
        team_members = [[item[0] for item in data] for (key, data) in team_groups]

        # Map the nodes on to the teams
        # Homogonization step here to avoid communication across nodes...
        unique_nodes = [node.split('.')[0] for node in list(set(pbs_nodes))]
        print("Running on nodes: " + ', '.join(unique_nodes))
        del pbs_nodes
        pbs_nodes = []
        for i_team in range(n_teams):
            pbs_nodes = pbs_nodes + ([unique_nodes[i_team % len(unique_nodes)]] * ncores_wrf_hydro)
        node_teams = [the_node // ncores_wrf_hydro for the_node in range(len(pbs_nodes))]
        node_team_seq = [ [node, team] for node,team in zip(pbs_nodes, node_teams)]
        
        node_team_seq.sort(key = operator.itemgetter(1))
        team_groups = itertools.groupby(node_team_seq, operator.itemgetter(1))
        team_nodes = [[item[0] for item in data] for (key, data) in team_groups]

        entry_cmd = config['run_experiment']['perturb_forcing']['noise_cmd']
        exit_cmd = None
        exe_cmd = config['run_experiment']['job_execution']['wrf_hydro']['exe_cmd']
        env = config['run_experiment']['job_execution']['wrf_hydro']['env']

        for team in range(n_teams):
            teams_dict.update({
                team: {
                    'members': team_members[team],
                    'nodes': team_nodes[team],
                    'entry_cmd': entry_cmd,
                    'exit_cmd': exit_cmd,
                    'exe_cmd': exe_cmd,
                    'env': env
                }
            })

        print('\nPBS_NODEFILE present: ')
        print('    ' + str(len(unique_nodes)) + ' nodes with')
        print('    ' + str(n_total_processors) + ' processors requested.')

        print('\nWRF-Hydro advance parallelization over:')
        print('    ' + str(n_teams) + ' concurrent teams each using')
        print('    ' + str(ncores_wrf_hydro) + ' processors.')

    else: 
        n_concurrent = ppn_use // ncores_wrf_hydro
        teams_dict = None
        print('\nWRF-Hydro advance parallelization over:')
        print('    ' + str(n_concurrent) + ' concurrent jobs each using')
        print('    ' + str(ncores_wrf_hydro) + ' processors.')

    print('\nDART Filter parallelization over:')
    print('    ' + str(ncores_wrf_hydro) + ' processors.')

    curr_ens_time = hydrodartpy.get_ensemble_time(run_dir)
    prev_ens_time = curr_ens_time - datetime.timedelta(hours=model_adv_hr)

    # Which restarts are being used? As implied by the inital ensemble
    use_hydro_rst = len(list(ini_dir.glob('member_000/HYDRO_RST*'))) > 0
    use_lsm_rst = len(list(ini_dir.glob('member_000/RESTART*'))) > 0
    use_param_rst = len(list(ini_dir.glob('member_000/parameter_restart*'))) > 0

    # -------------------------------------------------------
    # The analysis cycling loop
    while curr_ens_time < cycle_end_time:
        print(major_delim)
        print('Assimilation cycle at time: ', curr_ens_time.strftime(datestr_fmts['separated']))

        no_obs_adv_n_hrs = prep_obs_seq_for_filter(
            run_dir=run_dir,
            obs_dir=obs_dir,
            curr_ens_time=curr_ens_time,
            assim_window_hr=assim_window_hr,
            skip_missing_obs_hours=skip_missing_obs_hours,
            skip_missing_obs_days=skip_missing_obs_days
        )

        parameter_time_advance(
            run_dir=run_dir,
            prev_ens_time=prev_ens_time,
            curr_ens_time=curr_ens_time,
            ens_size=ens_size,
            ncores=ppn_use
        )

        if no_obs_adv_n_hrs == 0:
            skip_filter = False
            no_obs_adv_n_hrs = model_adv_hr
            
        else:
            skip_filter = True
            if no_obs_adv_n_hrs == 1:
                no_obs_adv_n_hrs = model_adv_hr

        if skip_filter:
            print(
                '\nObservations not present, advancing the model ' + \
                str(no_obs_adv_n_hrs) + ' hrs.\n'
            )
            
        else:
            prep_run_dir_for_filter(
                run_dir=run_dir,
                curr_ens_time=curr_ens_time,
                ens_size=ens_size,
                use_hydro_rst=use_hydro_rst,
                use_lsm_rst=use_lsm_rst,
                use_param_rst=use_param_rst,
                ncores=ppn_use
            )

            run_filter(
                run_dir=run_dir,
                nproc=ncores_dart,
                cmd=dart_exe_cmd
            )

            manage_filter_output(
                run_dir=run_dir,
                curr_ens_time=curr_ens_time,
                n_domains=use_hydro_rst + use_lsm_rst + use_param_rst
            )
            
        advance_ensemble(
            run_dir=run_dir,
            n_concurrent=n_concurrent,
            teams_dict=teams_dict,
            advance_n_hours=no_obs_adv_n_hrs
        )

        prev_ens_time, curr_ens_time = hydrodartpy.get_ensemble_time(run_dir, model_adv_hr)

    return 0


if __name__ == "__main__":

    run_dir = pathlib.Path(os.getcwd())
    run_filter_experiment(run_dir)
    sys.exit(0)
