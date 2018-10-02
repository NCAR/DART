# Setup a "HydroDartRun". For now this is a filter experiment.

from boltons.iterutils import remap
import code
import datetime
import f90nml
import os
import pathlib
from pprint import pprint
import shlex
import shutil
import subprocess
import sys
import warnings
import wrfhydropy

from .setup_dart import setup_dart
from .setup_experiment_tools import establish_config, replace_in_file, get_machine
from .setup_initial_ens import setup_initial_ens
from .setup_obs_prep import setup_obs_prep
from .setup_wrf_hydro import setup_wrf_hydro
from .setup_wrf_hydro_ens import setup_wrf_hydro_ens
from .setup_wrf_hydro_ens_job import setup_wrf_hydro_ens_job

def setup_experiment(config_file):

    config = establish_config(config_file)

    # The following list and function turn strings in the config
    # into pathlib.PosixPath objs.
    config_abs_paths_list = [
        'experiment_dir',
        'run_dir',
        'path',
        'hydro_file_list.txt',
        'param_file_list.txt',
        'dart_src',
        'wrf_hydro_src',
        'domain_src',
        'setup_py',
        'noise_function',
        'output_dir',
        'input_dir'
    ]


    def visit_abs_paths(path, key, value):
        if value is None:
            return True
        # Making a bet that the same keys at different hierarchical levels
        # will both be desired as pathlib.PosixPath objects.
        if key in config_abs_paths_list and type(value) is not dict:
            return key, pathlib.PosixPath(value)
        else:
            return True


    config = remap(config, visit_abs_paths)

    # JLM: Not sure this print is really worth it. I rarely look at it.
    # TODO(JLM): consider a verbose mode for this script.
    # print("The configuration file:")
    # pprint(config)
    # print('')

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
    # TODO(JLM): What to do if the dir exists? Use python mode 'w' for clobber?
    # Setting up and experiment should be a totally automated process,
    # so if anything needs to change (in the yaml), then just resetup and rerun.

    # Copy the config file to the experiment dir
    # TODO(JLM): Not sure the best approach here: Rename this copied file? This file does not match
    # the internal representation, it reflects state of the original setup. One approach:
    # copy it here and rename as a dot file: ".original.<file_name>" then stash that
    # path in the internal config file representation. Maybe also stash the path to the original file?
    file_dir = os.path.dirname(config_file)
    file_base = os.path.basename(config_file)
    config_file_copy_name = file_dir + '.original.' + file_base
    shutil.copy(
        config_file,
        config['experiment']['experiment_dir'] / os.path.basename(config_file_copy_name)
    )


    # -----------------------------------
    # The major phases of setting up an experiment.

    # DART Simulation
    dart_compile = setup_dart(config)

    # WRF-Hydro deterministic
    wrf_hydro_sim = setup_wrf_hydro(config)

    # WRF-Hydro Ensemble Setup
    wrf_hydro_ens_sim = setup_wrf_hydro_ens(config, wrf_hydro_sim)

    # Initial ensemble setup and enemble configuration to it.
    wrf_hydro_ens_sim = setup_initial_ens(config, wrf_hydro_ens_sim)

    # The initial ensemble files will be copied to the run dir after compose.
    # Before compose, we have to change the paths to be used for the restart files.
    # shorten these names...
    nlsm_off = 'noahlsm_offline'
    rst_fnm = 'restart_filename_requested'
    att_tuple = ('base_hrldas_namelist', nlsm_off, rst_fnm)
    values = [pathlib.Path(mm.base_hrldas_namelist[nlsm_off][rst_fnm]).name \
              for mm in wrf_hydro_ens_sim.members]
    wrf_hydro_ens_sim.set_member_diffs(att_tuple, values)

    # Set the ensemble to use the copys in the run dir.
    att_tuple = ('base_hydro_namelist', 'hydro_nlist', 'restart_file')
    values = [pathlib.Path(mm.base_hydro_namelist['hydro_nlist']['restart_file']).name \
              for mm in wrf_hydro_ens_sim.members]
    wrf_hydro_ens_sim.set_member_diffs(att_tuple, values)

    # Add an intial job
    wrf_hydro_ens_sim = setup_wrf_hydro_ens_job(config, wrf_hydro_ens_sim)

    # ###################################
    # WRF-Hydro Ensemble Run: establish the run directory.
    print('Compose wrfhydro.EnsembleSim object.')
    config['experiment']['run_dir'].mkdir(parents=True)
    os.chdir(config['experiment']['run_dir'])
    wrf_hydro_ens_sim.compose(
        rm_members_from_memory=False,
        check_nlst_warn=True
    )


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

    # Delete the members from memory at the end of the experiment setup. (Though it's only needed
    # one place and probably could use the job for that purpose.)

    # # ###################################
    # # HydroDartRun Object
    # print('HydroDartRun object.')
    # # TODO(JLM): filter/run/etc?? I'd be more in favor of calling it a run and specifying
    # # the executable task to the job. But then we are mixing types of jobs.
    # # Perhaps the different executables have different job lists?
    # hydro_dart_run = dartclasses.HydroDartRun(
    #     dart_compile=dart_compile,
    #     wrf_hydro_ens_sim=wrf_hydro_ens_sim,
    #     config=config
    # )
    # hydro_dart_run.pickle()

    # #del dart_compile, wrf_hydro_ens_sim
    # #wrf_hydro_model_pickle.unlink()  ## TODO(JLM): Temporary comment only for testing.

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


    ###########################################################################
    # JLM keep refactoring below here.... 


    # #######################################################
    # With the run directory created, we can stage the ensemble to take the initial ensemble.

    for mm in range(config['ensemble']['size']):

        rst_path = config['initial_ens']['path'] / ('member_' + "%03d" % (mm,))
        member = wrf_hydro_ens_sim.members[mm]
        run_path = config['experiment']['run_dir'] / ('member_' + "%03d" % (mm,))

        hydro_rst_old = list(run_path.glob('HYDRO_RST*'))
        _ = [rr.unlink() for rr in hydro_rst_old]
        hydro_rst = sorted(rst_path.glob('HYDRO_RST*'))
        if len(hydro_rst) == 0:
                raise ValueError('No HYDRO_RST file present.')
        if len(hydro_rst) > 1:
            warnings.warn("Multiple HYDRO_RST files supplied, using the first in the list.")
        # Do NOT symlink the initial ensemble, copy it as the files are overwritten
        shutil.copy(hydro_rst[0], str(run_path / hydro_rst[0].name))
        (run_path / hydro_rst[0].name).chmod(0o755)

        lsm_rst_old = list(run_path.glob("RESTART.*"))
        _ = [rr.unlink() for rr in lsm_rst_old]
        lsm_rst = sorted(rst_path.glob('RESTART*'))
        if len(lsm_rst) == 0:
            forc_typ = member.base_hrldas_namelist['wrf_hydro_offline']['forc_typ']
            if forc_typ not in [9, 10]:
                raise ValueError('No RESTART file present (and forc_typ is neither 9 nor 10).')
        if len(lsm_rst) > 1:
            warnings.warn("Multiple RESTART files supplied, using the first in the list.")
        # Do NOT symlink the initial ensemble, copy it as the files are overwritten
        if len(lsm_rst) > 0:
            shutil.copy(lsm_rst[0], str(run_path / lsm_rst[0].name))
            (run_path / lsm_rst[0].name).chmod(0o755)

        param_rst = list(rst_path.glob('parameters*nc'))
        if len(param_rst) > 0:
            param_rst = param_rst[0]
            # Do NOT symlink the initial ensemble, copy it as the files are overwritten
            shutil.copy(param_rst, str(run_path / param_rst.name))
            (run_path / param_rst.name).chmod(0o755)


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


    # ###################################
    # Place scripts into the run dir.
    # Should it be checke that this file (setup_experiment.py) is in
    # wrf_hydro_dart/models/wrf_hydro/shell_scripts?
    print("Staging scripts.")

    # Various scripts (tied to run_experiment)
    script_list = ['advance_ensemble.py', 'get_ensemble_time.py', 'set_obs_seq_times.py']
    for ss in script_list:
        adv_ens_script_src = config['dart']['dart_src'] / \
                             ('models/wrf_hydro/hydro_dart_py/hydrodartpy/core/' + ss)
        adv_ens_script_link = config['experiment']['run_dir'] / ss
        adv_ens_script_link.symlink_to(adv_ens_script_src)


    # Noise scripts
    perturb_forcing = config['run_experiment']['perturb_forcing']
    if perturb_forcing['perturb']:
        for ff in perturb_forcing['noise_function_files']:
            for mm in range(config['ensemble']['size']):
                the_link = config['experiment']['run_dir'] / \
                           (('member_' + "%03d" % (mm,)) + '/' + os.path.basename(ff))
                the_link.symlink_to(ff)


    # DART binaries / namelist
    print('Staging DART executables and input.nml.')
    model_work_path = ( config['experiment']['experiment_dir'] / config['dart']['build_dir'] /
                        'models/wrf_hydro/work' )

    dart_exe_names = config['run_experiment']['dart']['exes']
    for dd in dart_exe_names:
        dart_exe = model_work_path / dd
        dart_exe_link = config['experiment']['run_dir'] / dd
        dart_exe_link.symlink_to(dart_exe)

    input_nml_patches = config['run_experiment']['dart']['input_nml_patches']

    dart_input_nml = model_work_path / 'input.nml'
    dart_input_nml_copy = config['experiment']['run_dir'] / 'input.nml'

    # write model namelists to the run dir (needed by dart to run filter)
    model_nlsts = ['base_hydro_namelist', 'base_hrldas_namelist'] # 0 for testing if 1 is necessary
    model_output_nlsts = ['hydro.namelist', 'namelist.hrldas']
    for nn,oo in zip(model_nlsts, model_output_nlsts):
        nlst = wrf_hydro_sim.__dict__[nn]
        for k0 in nlst.keys():
            for k1 in nlst[k0].keys():
                if type(nlst[k0][k1]) is str:
                    nlst[k0][k1] = nlst[k0][k1].replace('./','./member_000/')
        f90nml.namelist.Namelist(nlst).write(config['experiment']['run_dir'] / oo)


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
            if k1 in ['ens_size', 'num_output_obs_members']:
                input_nml[k0][k1] = int(config['ensemble']['size'])

    input_nml.write(dart_input_nml_copy)    


    #######################################################

    # Setup run_filter_experiment.csh_template
    # Copy to run_dir then preprocess the template.
    run_experiment_template = \
        config['dart']['dart_src'] / 'models/wrf_hydro/shell_scripts/run_filter_experiment.csh.template'
    run_experiment_specific = config['experiment']['run_dir'] / 'run_filter_experiment.csh'
    shutil.copy(str(run_experiment_template), str(run_experiment_specific))

    replace_in_file(
        run_experiment_specific,
        'EXPERIMENT_DIRECTORY_TEMPLATE',
        str(config['experiment']['run_dir'])
    )

    replace_in_file(
        run_experiment_specific,
        'OBSERVATION_DIR_TEMPLATE',
        str(config['observation_preparation']['all_obs_dir'])
    )

    # Ensemble advance
    advance_cmd = "`python advance_ensemble.py --hold True"
    config_perturb_forcing = config['run_experiment']['perturb_forcing']
    if config_perturb_forcing['perturb']:
        advance_cmd += " --job_entry_cmd '" + config_perturb_forcing['noise_cmd'] + "'`\n"
        advance_cmd += 'set ens_adv_exit_code = $? \n'

    replace_in_file(run_experiment_specific, 'PYTHON_ADVANCE_TEMPLATE', advance_cmd)

    # Call filter again
    if config['run_experiment']['wrf_hydro_ens_advance']['with_filter']:
        filter_again_cmd =  'echo "ens_adv_exit_code: $ens_adv_exit_code" \n'
        filter_again_cmd += 'if ( $ens_adv_exit_code == 1 ) then \n'
        filter_again_cmd += '    echo "Ensemble advance failed, exiting." \n'
        filter_again_cmd += '    rm .filter_not_complete" \n'
        filter_again_cmd += 'exit 1 \n'
        filter_again_cmd += 'endif \n'
        filter_again_cmd += './run_filter_experiment.csh &\n'
    else:
        filter_again_cmd =  'set next_filter_afterok = `echo "$next_filter_afterok"'
        filter_again_cmd += ' | cut -d " " -f1 | tr -d " "`\n'
        filter_again_cmd += 'qsub -W depend=afterok:$next_filter_afterok run_filter_experiment.csh\n'
        filter_again_cmd += 'qrls $next_filter_afterok'

    replace_in_file(run_experiment_specific, 'FILTER_AGAIN_TEMPLATE', filter_again_cmd)


    # PBS tmp dir on cheyenne
    tmp_dir_cmd = 'setenv TMPDIR /glade/scratch/$USER/temp\nmkdir -p $TMPDIR'
    replace_in_file(run_experiment_specific, 'TMP_DIR_TEMPLATE', tmp_dir_cmd)

    # PBS options from config
    # Short-hand
    dart_sched = config['run_experiment']['dart']['scheduler']
    config_time = config['run_experiment']['time']

    # The easy ones.
    replace_in_file(run_experiment_specific, 'JOB_NAME_TEMPLATE', dart_sched['job_name'])
    replace_in_file(run_experiment_specific, 'ACCOUNT_TEMPLATE', dart_sched['account'])
    replace_in_file(run_experiment_specific, 'EMAIL_WHO_TEMPLATE', dart_sched['email_who'])
    replace_in_file(run_experiment_specific, 'EMAIL_WHEN_TEMPLATE', dart_sched['email_when'])

    # Model end time
    end_time_fmt = datetime.datetime.strptime(config_time['end_time'],'%Y-%m-%d_%H:%M')
    end_time_fmt = end_time_fmt.strftime('%Y%m%d%H%M')
    replace_in_file(run_experiment_specific, 'END_DATE_TEMPLATE', end_time_fmt)

    # Model advance time
    replace_in_file(run_experiment_specific, 'ADV_MODEL_HRS_TEMPLATE', config_time['advance_model_hours'])

    # Wall time
    dart_walltime = dart_sched['walltime']
    if len(dart_walltime.split(':')) == 2:
        dart_walltime = dart_walltime + ':00'
    dart_walltime = 'walltime=' + dart_walltime
    replace_in_file(run_experiment_specific, 'WALLTIME_TEMPLATE', dart_walltime)

    # If cheyenne, we know ppn (proc per node) <= 36.
    this_machine = get_machine()
    if this_machine != 'cheyenne':
        warnings.warn("Dont YET know how to handle schduler on machines other than cheyenne.")
    if dart_sched['ppn_max'] > 36:
        raise ValueError("Cheyenne has maximum ppn of 36.")

    # Do the nodes vs total number of processors match/make sense?
    if dart_sched['nnodes']*dart_sched['ppn_max'] < dart_sched['nproc']:
        raise ValueError(
            "run_experiment:dart:scheduler:nproc is the TOTAL number of processors requested."
        )

    # pbs "-l select=" : Distribute the processors over the nodes how to split?
    dart_nproc_last_node = \
        (dart_sched['nproc'] - (dart_sched['nnodes'] * dart_sched['ppn_max'])) % dart_sched['ppn_max']

    if dart_nproc_last_node > 0:
        if dart_nproc_last_node >= dart_sched['ppn_max']:
            raise ValueError('nproc - (nnodes * ppn) = {0} >= ppn'.format(dart_nproc_last_node))

    if dart_nproc_last_node == 0 or int(dart_sched['nnodes']) < 2:

        if config['run_experiment']['wrf_hydro_ens_advance']['with_filter']:
            prcstr = "select=1:ncpus=36:mpiprocs=36"
        else:
            prcstr = "select={0}:ncpus={1}:mpiprocs={1}"
            if int(dart_sched['nnodes']) < 2:
                prcstr = prcstr.format(dart_sched['nnodes'], dart_sched['nproc'])
            else:
                prcstr = prcstr.format(dart_sched['nnodes'], dart_sched['ppn_max'])

    else:
        prcstr = "select={0}:ncpus={1}:mpiprocs={1}+1:ncpus={2}:mpiprocs={2}"
        prcstr = prcstr.format(int(dart_sched['nnodes'])-1,
                               dart_sched['ppn_max'],
                               dart_nproc_last_node)

    pbs_select_cmd = prcstr
    replace_in_file(run_experiment_specific, 'PBS_SELECT_TEMPLATE', pbs_select_cmd)

    # queue
    if dart_sched['nproc'] <= 18:

        if not config['run_experiment']['wrf_hydro_ens_advance']['with_filter']:

            # Shared queue
            if dart_sched['queue'] != 'share':
                warnings.warn(
                    'You have not selected share queue but are requesting <= 18 cores.\n' +
                    'DART jobs will be sent to the shared queue.'
                )
                share_use_array_cmd = 'setenv MPI_USE_ARRAY false'
                launch_cmd = '"mpirun `hostname` -np {0}"'.format(dart_sched['nproc'])
                dart_queue = 'share'

        else:

            share_use_array_cmd = 'setenv MPI_USE_ARRAY false'
            launch_cmd = '"mpirun `hostname` -np {0}"'.format(dart_sched['nproc'])
            dart_queue = dart_sched['queue']

    else:

        share_use_array_cmd = ''
        launch_cmd = 'mpiexec_mpt'
        dart_queue = dart_sched['queue']

    # Apply the choices
    replace_in_file(run_experiment_specific, 'SHARE_USE_ARRAY_TEMPLATE', share_use_array_cmd)
    replace_in_file(run_experiment_specific, 'LAUNCH_CMD_TEMPLATE', launch_cmd)
    replace_in_file(run_experiment_specific, 'QUEUE_TEMPLATE', dart_queue)

    # If with_filter, change the top-level script to put the stdout/err to disk (not PBS buffer).
    if config['run_experiment']['wrf_hydro_ens_advance']['with_filter']:

        ## TODO(JLM): this should be a script in a file that we edit... 
        submit_filter_file = run_experiment_specific.parent / 'submit_filter_experiment.csh'
        subprocess.run(
            (
                'echo "#!/bin/tcsh" >> ' + str(submit_filter_file) + 
                ' && egrep "^#PBS" ' + str(run_experiment_specific) + ' >> ' + str(submit_filter_file) +
                ' && echo "touch .filter_not_complete" >> ' + str(submit_filter_file) +
                ' && echo "./run_filter_experiment.csh >& submit_filter_experiment.stdeo" >> ' +
                str(submit_filter_file) +
                ' && echo "while (\`ls .filter_not_complete | wc -l\`)" >> ' + str(submit_filter_file) +
                ' && echo "    sleep 20" >> ' + str(submit_filter_file) +
                ' && echo "end" >> ' + str(submit_filter_file) +
                ' && echo "exit 0" >> ' + str(submit_filter_file)
            ),
            shell=True,
            executable='/bin/bash'
        )
        submit_filter_file.chmod(0o755)



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
    # Write the EnsembleSimulation object to disk after deleting its members from memory.
    wrf_hydro_ens_sim.rm_members()
    wrf_hydro_ens_sim.pickle(config['experiment']['run_dir'] / 'WrfHydroEnsembleSim.pkl')


    # ###################################
    # User instructions/guidance

    print(
    """

--------------------------------------------------------

The experiment has been established in 
    experiment_dir = {0}
    run_dir        = {1}

To run a filter experiment:
    cd {1}
    qsub {2}_filter_experiment.csh
    """.format(
        config['experiment']['experiment_dir'],
        config['experiment']['run_dir'],
        'submit' if config['run_experiment']['wrf_hydro_ens_advance']['with_filter'] else 'run'
        )
    )

    return 0


