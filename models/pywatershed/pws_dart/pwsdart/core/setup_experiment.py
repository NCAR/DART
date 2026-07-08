import datetime
import f90nml
import os
import pathlib
from pprint import pprint
import shutil
import subprocess
import warnings

from .setup_dart import setup_dart
from .setup_experiment_tools import establish_config, replace_in_file, get_machine
from .setup_initial_ens import setup_initial_ens
from .setup_obs_prep import setup_obs_prep
from .setup_pywatershed_ens import setup_pywatershed_ens
from .setup_pywatershed_forcings import setup_pywatershed_forcings
from .setup_pywatershed_ens_job import setup_pywatershed_job


def setup_experiment(config_file):

    config = establish_config(config_file)

    # Report the important dirs in the experiment structure.
    dirs_string = """
expdir={0}
rundir={1}
inidir={2}
allobsdir={3}
    """.format(
        config['experiment']['experiment_dir'],
        config['experiment']['run_dir'],
        config['initial_ens']['path'],
        config['observation_preparation']['all_obs_dir']
        )
    print(dirs_string)

    # -----------------------------------
    # Setup experiment_dir
    print('Create experiment directory.')
    config['experiment']['experiment_dir'].mkdir(
        parents=True,
        exist_ok=config['experiment']['experiment_dir_mode'] == 'w'
    )
    # Setting up and experiment should be a totally automated process,
    # so if anything needs to change (in the yaml), then just resetup and rerun.

    # Copy the config file to the experiment dir
    # copy it here and rename as a dot file: ".original.<file_name>" then stash that
    # path in the internal config file representation. Maybe also stash the path to the original file?
    file_dir = os.path.dirname(config_file)
    file_base = os.path.basename(config_file)
    config_file_copy_name = 'original.' + file_base
    shutil.copy(
        config_file,
        config['experiment']['experiment_dir'] / os.path.basename(config_file_copy_name)
    )


    # -----------------------------------
    # The major phases of setting up an experiment.

    # DART Simulation
    dart_compile = setup_dart(config)

    # Generate an ensemble of the forcings
    pywatershed_ens_forcing = setup_pywatershed_forcings(config)

    # Generate an ensemble of parameters from a single file for pywatershed 
    pywatershed_ens_param = setup_pywatershed_ens(config)

    # Gererate an ensemble of initial state from a single restart file 
    pywatershed_ens_sim = setup_initial_ens(config)    

    # Observation Preparation
    pywatershed_obs_prep = setup_obs_prep(config)
   
    # Set up the pywatershed ensemble based on single domain, ensemble parameters, and ensemble initial condition 
    #pywatershed_sim = setup_pywatershed(config)

    # ###################################
    # Pywatershed Ensemble Run: establish the run directory.
    config['experiment']['run_dir'].mkdir(parents=True)
    os.chdir(config['experiment']['run_dir'])


    # this is where it creates a folder for each member I guess
    for i in range(config['ensemble']['size']):
        member_name = f"member_{i:03d}"
        member_path = config['experiment']['run_dir'] / f"member_{i:03d}"
        member_path.mkdir(exist_ok=True)
        print(f"Created: {member_path}")

    # Since parameter_restart files are not part of the model, they are not captured by compose.
    # Put them in the run directory now.
    # If there are no parameter restarts in the inital ensemble, they are skipped.
    member_dirs = config['initial_ens']['path'].glob('member_*')
    for mm in member_dirs:
        param_rst_file = list(mm.glob('parameter_restart.*'))
        if len(param_rst_file) == 0:
            continue
        # There will only be one of these by design.
        param_rst_file = param_rst_file[0]
        new_file = config['experiment']['run_dir'] / (mm.name + '/' + param_rst_file.name)
        shutil.copy(str(param_rst_file), str(new_file))
        new_file.chmod(0o755)

    # ###################################
    # Now that the run_dir is established, point the experiment and run dirs to each other.
    link_to_run_dir = config['experiment']['experiment_dir'].joinpath('run_dir')
    if link_to_run_dir.is_symlink():
        link_to_run_dir.unlink()
    link_to_run_dir.symlink_to(config['experiment']['run_dir'])

    link_to_exp_dir = config['experiment']['run_dir'].joinpath('experiment_dir')
    if link_to_exp_dir.is_symlink():
        link_to_exp_dir.unlink()
    link_to_exp_dir.symlink_to(config['experiment']['experiment_dir'])


    # #######################################################
    # With the run directory created, we can stage the ensemble to take the initial ensemble.

    for mm in range(config['ensemble']['size']):

        rst_path = config['initial_ens']['path'] / ('member_' + "%03d" % (mm,))
        run_path = config['experiment']['run_dir'] / ('member_' + "%03d" % (mm,))

        channel_rst_old = list(run_path.glob('*-outflow_ts.nc'))
        _ = [rr.unlink() for rr in channel_rst_old]
        channel_rst = sorted(rst_path.glob('*-outflow_ts.nc'))
        if len(channel_rst) == 0:
                raise ValueError('No *-outflow_ts.nc file present.')
        if len(channel_rst) > 1:
            warnings.warn("Multiple outflow_ts.nc files supplied, using the first in the list.")
        # Do NOT symlink the initial ensemble, copy it as the files are overwritten
        shutil.copy(channel_rst[0], str(run_path / channel_rst[0].name))
        (run_path / channel_rst[0].name).chmod(0o755)

        hru_rst_old = list(run_path.glob("*-gwres_stor.nc"))
        _ = [rr.unlink() for rr in hru_rst_old]
        hru_rst = sorted(rst_path.glob('*-gwres_stor.nc'))
        if len(hru_rst) > 1:
            warnings.warn("Multiple gwres_stor.nc files supplied, using the first in the list.")
        if len(hru_rst) > 0:
            shutil.copy(hru_rst[0], str(run_path / hru_rst[0].name))
            (run_path / hru_rst[0].name).chmod(0o755)

        param_rst = list(rst_path.glob('parameters*nc'))
        if len(param_rst) > 0:
            param_rst = param_rst[0]
            # Do NOT symlink the initial ensemble, copy it as the files are overwritten
            shutil.copy(param_rst, str(run_path / param_rst.name))
            (run_path / param_rst.name).chmod(0o755)


        # AREZOO
        # ket s copy the single restart file that we do not perturb 
        rst_path = pathlib.Path(config['initial_ens']['from_filter']['input_nml']['filter_nml']['input_state_file_list']['seg_inflow_restart_file_list.txt'])
        dest_path = run_path / rst_path.name
        dest_path.symlink_to(rst_path)

    # If the intial_ens directory contains a from_filter/ dir, then stage the inflation files.
    from_filter_path = config['initial_ens']['path'] / 'from_filter'
    if from_filter_path.exists():

        #Opportunistically stage the inflation files: if they are there, copy them to the run_dir.

        prior_inf_files = from_filter_path.glob('output_*priorinf*')
        for pr in prior_inf_files:
            new_name = pr.name.replace('output', 'input')
            _ = shutil.copy(pr, config['experiment']['run_dir'] / new_name)

        post_inf_files = from_filter_path.glob('output_*postinf*')
        for po in post_inf_files:
            new_name = po.name.replace('output', 'input')
            _ = shutil.copy(po, config['experiment']['run_dir'] / new_name)

        # Adding the stage of the hybrid weighting coefficient files
        hybrid_wght_files = from_filter_path.glob('output_*hybridweight*')
        for hyb in hybrid_wght_files:
            new_name = hyb.name.replace('output', 'input')
            _ = shutil.copy(hyb, config['experiment']['run_dir'] / new_name)

    # ###################################
    # Place scripts into the run dir.
    print("Staging scripts.")

    # Various scripts (tied to run_experiment)
    #script_list = ['advance_ensemble.py', 'get_ensemble_time.py', 'set_obs_seq_times.py']
    script_list = ['run_filter_experiment.py']
    for ss in script_list:
        adv_ens_script_src = config['dart']['dart_src'] / \
                             ('models/pywatershed/pws_dart/pwsdart/core/' + ss)
        adv_ens_script_link = config['experiment']['run_dir'] / ss
        adv_ens_script_link.symlink_to(adv_ens_script_src)

    # DART binaries / namelist
    print('Staging DART executables and input.nml.')
    model_work_path = ( config['experiment']['experiment_dir'] / config['dart']['build_dir'] /
                        'models/pywatershed/work' )

    dart_exe_names = config['run_experiment']['dart']['exes']
    for dd in dart_exe_names:
        dart_exe = model_work_path / dd
        dart_exe_link = config['experiment']['run_dir'] / dd
        dart_exe_link.symlink_to(dart_exe)

    input_nml_patches = config['run_experiment']['dart']['input_nml_patches']

    dart_input_nml = model_work_path / 'input.nml'
    dart_input_nml_copy = config['experiment']['run_dir'] / 'input.nml'

    # Also need domain files. Get the entire config dir.
    domain_file = pathlib.Path(config['pywatershed']['domain_dir'] + '/parameters_dis_seg_app.nc')
    config_file_dest = pathlib.Path(config['experiment']['run_dir'] / 'parameters_dis_seg_app.nc')
    config_file_dest.symlink_to(
            domain_file
        )

    domain_file = pathlib.Path(config['pywatershed']['domain_dir'] + '/parameters_PRMSChannel.nc')
    config_file_dest = pathlib.Path(config['experiment']['run_dir'] / 'parameters_PRMSChannel.nc')
    config_file_dest.symlink_to(
            domain_file
        )

    domain_file = pathlib.Path(config['pywatershed']['domain_dir'] + '/parameters_dis_hru.nc')
    config_file_dest = pathlib.Path(config['experiment']['run_dir'] / 'parameters_dis_hru.nc')
    config_file_dest.symlink_to(
            domain_file
        )

    # Move climatology file list also to the rundir
    static_file = model_work_path / 'static_file_list.txt'
    if static_file.exists():
        shutil.copy(static_file, config['experiment']['run_dir'] / 'static_file_list.txt')

    # input.nml patches and checks
    print('Apply namelist patches and checks.')
    input_nml = f90nml.read(dart_input_nml)

    # Ensure consistency with initial ensemble from filter if being used.
    if config['initial_ens']['from_filter']['use']:
        initial_ens_in_nml = f90nml.read(config['initial_ens']['path'] / 'from_filter/input.nml')

        nml_list = ['model_nml',              'model_nml',
                    'filter_nml',             'filter_nml']
        val_list = ['domain_order',           'domain_shapefiles',
                    'output_state_file_list', 'input_state_file_list']

        for nml, val in zip(nml_list, val_list):

            if val not in input_nml_patches[nml].keys():
                # Handle this one exception.
                if nml == 'filter_nml' and val == 'output_state_file_list':
                    input_nml[nml][val] = initial_ens_in_nml[nml]['input_state_file_list']
                else:
                    input_nml[nml][val] = initial_ens_in_nml[nml][val]
            else:
                raise ValueError('Using inital ensemble from filter but supplying '
                                 'filter a ' + nml + ':' + val + ' patch.')

    if input_nml_patches:
        for k0 in input_nml_patches.keys():
            for k1 in input_nml_patches[k0].keys():
                if k1 == 'ens_size':
                    if int(input_nml_patches[k0][k1]) != int(config['ensemble']['size']):
                        raise ValueError("input_nml_patches ens_size do not match config:ensemble:size")
                input_nml[k0][k1] = input_nml_patches[k0][k1]

    # verify all ens_size in the entire namelist... 
    # Since there only two levels, just write a double for loop.
    for k0 in input_nml.keys():
        for k1 in input_nml[k0].keys():
            if k1 in ['ens_size', 'num_output_obs_members', 'num_output_state_members']:
                input_nml[k0][k1] = int(config['ensemble']['size'])

    input_nml.write(dart_input_nml_copy)    


    # lets just make the member_XXX directories 
    pws_job = setup_pywatershed_job(config)
    #######################################################

    # Setup run_filter_experiment.csh_template
    # Copy to run_dir then preprocess the template.

    # job execution options from config
    # Short-hand
    job_exe = config['run_experiment']['job_execution']
    sched = job_exe['scheduler']
    config_time = config['run_experiment']['time']

    # If cheyenne, we know ppn (proc per node) <= 36.
    this_machine = get_machine()
    if this_machine != 'derecho':
        warnings.warn("Dont YET know how to handle schduler on machines other than derecho.")
    # TODO: check that if cheynne the number is 36 and
    # also check this against job_execution:moachine:ppn
    if sched['ppn_use'] > 128:
        raise ValueError("Derecho has maximum ppn of 128.")

    # Write a submission script
    submit_script_template = config['dart']['dart_src'] / \
                             ('models/pywatershed/pws_dart/pwsdart/core/'
                              'submission_scripts/submit_filter_experiment.sh')
    submit_script_specific = config['experiment']['run_dir'] / 'submit_filter_experiment.sh'
    shutil.copy(submit_script_template, submit_script_specific)

    pbs_select_specific = "select={0}:ncpus={1}:mpiprocs={1}".format(
        sched['nnodes'],
        sched['ppn_use']
    )

    walltime_specific = sched['walltime']
    if len(walltime_specific.split(':')) == 2:
        walltime_specific = walltime_specific + ':00'
    walltime_specific = 'walltime=' + walltime_specific

    template_dict = {
        'JOB_NAME_TEMPLATE': sched['job_name'],
        'PBS_SELECT_TEMPLATE': pbs_select_specific,
        'WALLTIME_TEMPLATE': walltime_specific,
        'ACCOUNT_TEMPLATE': sched['account'],
        'QUEUE_TEMPLATE': sched['queue'],
        'EMAIL_WHEN_TEMPLATE': sched['email_when'],
        'EMAIL_WHO_TEMPLATE': sched['email_who']
    }

    for key, value in template_dict.items():
        replace_in_file(submit_script_specific, key, value)


    # ###################################
    # Observation preparation
    _ = setup_obs_prep(config)


    # ###################################
    # Sanity checks

    # 1) there is danger if the filter_nml: domain_order and domain_shapefiles are different for
    # different parts of the experiment. Check the consistency between use in initial ensemble
    # generation and the filter exeperiment.
    # Check these after they are written to disk.
    if config['initial_ens']['from_filter']['use']:

        initial_ens_in_nml = f90nml.read(config['initial_ens']['path'] / 'from_filter/input.nml')
        filter_in_nml = f90nml.read(config['experiment']['run_dir'] / 'input.nml')

        filter_domain_order = filter_in_nml['model_nml']['domain_order']
        ini_ens_domain_order = initial_ens_in_nml['model_nml']['domain_order']

        filter_domain_shapefiles = filter_in_nml['model_nml']['domain_shapefiles']
        ini_ens_domain_shapefiles = initial_ens_in_nml['model_nml']['domain_shapefiles']

        message=''
        if (filter_domain_order != ini_ens_domain_order):
            message += "\nInitial ensemble and filter input.nml model_nml domain_order do not match.\n"
            message += "filter (run dir) domain_order: " + ''.join(filter_domain_order) + '\n'
            message += "inital ensemble  domain_order: " + ''.join(ini_ens_domain_order) + '\n\n'

        if (filter_domain_shapefiles != ini_ens_domain_shapefiles):
            message += "Initial ensemble and filter input.nml model_nml domain_shapefiles do not match.\n"
            message += "filter (run dir) domain_shapefiles: " + ''.join(filter_domain_shapefiles) + '\n'
            message += "inital ensemble  domain_shapefiles: " + ''.join(ini_ens_domain_shapefiles) + '\n\n'

        if message != '':
            message += ('Severe problems likely with domain order and shapefiles. Please check \n'
                        'for consistency in your yaml file. Some settings may be taken from \n'
                        'default input.nml. Use run_dir with EXTREME CAUTION. \n')
            raise ValueError(message)

    # ###################################
    # User instructions/guidance
    print(
    """ 
    
#-------------------------------------------------------

The experiment has been established in 
    experiment_dir = {0}
    run_dir        = {1}

To run a filter experiment:
    cd {1}
    qsub submit_filter_experiment.sh
    """.format(
            config['experiment']['experiment_dir'],
            config['experiment']['run_dir']
        )
    )

    return 0
  
