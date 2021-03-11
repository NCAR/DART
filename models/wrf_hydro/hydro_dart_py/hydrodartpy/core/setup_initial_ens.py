from copy import deepcopy
import f90nml
import json
import os
import pathlib
import shlex
import shutil
import subprocess
import warnings
import wrfhydropy

from .create_parameter_restart_file import create_parameter_restart_file
from .obs_seq_dummy import obs_seq_dummy
from .setup_experiment_tools import get_top_level_dir_from_config

# TODO(JLM): this could be further broken into param restart creation
# and "from filter" parts, though threre are some interdependencies.

def setup_initial_ens(config, wrf_hydro_ens_sim):

    # This leverages the same setup object for advancing the model.

    # -------------------------------------------------------
    # This section is established needed variables for either of
    # 1) create param file or
    # 2) create ensemble from filter.
    if config['initial_ens']['param_restart_file']['create'] or \
       config['initial_ens']['from_filter']['create']:

        init_ens_dir = pathlib.PosixPath(config['initial_ens']['path'])
        init_ens_dir.mkdir(parents=True, exist_ok=True)

        # Short-hand
        input_state_file_list = \
            config['initial_ens']['from_filter']['input_nml']['filter_nml']['input_state_file_list']

        # TODO(JLM): Dont allow clobber? Could wait until write, but maybe favor preemption.
        hydro_patches_json = config['wrf_hydro']['domain_src'] / 'hydro_namelist_patches.json'
        json_namelist = json.load(hydro_patches_json.open())
        base_namelist = deepcopy(json_namelist['base'])
        config_patches = deepcopy(json_namelist[config['wrf_hydro']['model_config']])
        #Update the base namelist with the config patches
        config_namelist = wrfhydropy.core.namelist.dict_merge(base_namelist,config_patches)
        domain_config_rst_path = config['wrf_hydro']['domain_src'] / \
            pathlib.Path(config_namelist['hydro_nlist']['restart_file']).parent

        if input_state_file_list['hydro_file_list.txt'] is None:
            input_state_file_list['hydro_file_list.txt'] = \
                wrf_hydro_ens_sim.members[0].base_hydro_namelist['hydro_nlist']['restart_file']

    # -------------------------------------------------------
    # Create a parameter restart file
    if config['initial_ens']['param_restart_file']['create']:

        print('Create parameter restart files.')
        
        if input_state_file_list['param_file_list.txt'] is not None:
            raise ValueError('You are requesting to create a parameter restart file '
                             'and specifying an existing one.')

        create_param_file = config['initial_ens']['param_restart_file']

        param_rst_file = create_parameter_restart_file(
            out_path=init_ens_dir,
            out_mode=create_param_file['mode'],
            hydro_rst_file=input_state_file_list['hydro_file_list.txt'],
            restart_file=None,
            existing_variables=create_param_file['existing_variables'],
            new_variables=create_param_file['new_variables'],
            values=create_param_file['values']
        )

        input_state_file_list['param_file_list.txt'] = param_rst_file
        
        # protect the file
        param_rst_file.chmod(0o444)


    # -------------------------------------------------------
    # Initial ensemble creation from filter
    if config['initial_ens']['from_filter']['create']:

        print('Create initial ensemble using filter.')

        # make a place to filter in the init_ens_filter_dir
        init_ens_filter_dir = pathlib.PosixPath(config['initial_ens']['path'] / 'from_filter')
        init_ens_filter_dir.mkdir()

        # Copy filter and input_nml
        models_dir = [ww for ww in config['dart']['work_dirs'] if 'models' in ww]
        if len(models_dir) != 1:
            raise ValueError("Currently only supporting individual models")
        model_work_path = (config['experiment']['experiment_dir'] /
                           config['dart']['build_dir'] /
                           models_dir[0])
        shutil.copy(model_work_path / 'filter', init_ens_filter_dir / 'filter')
        shutil.copy(model_work_path / 'input.nml', init_ens_filter_dir / 'input.nml')

        # Bring in input_nml, apply patches.
        init_ens_filter_nml = f90nml.read(init_ens_filter_dir / 'input.nml')

        # This one is handled differently in the YAML than in the fortran.
        init_ens_filter_nml_update = deepcopy(
            config['initial_ens']['from_filter']['input_nml'])
        _ = init_ens_filter_nml_update['filter_nml'].pop('input_state_file_list')
        check_keys = init_ens_filter_nml_update['model_nml'].keys()
        
        for bb in ['input_state_file_list', 'domain_order', 'domain_shapefiles']:
            if bb in check_keys:
                warnings.warn(
                    "Key " + bb + " in initial_ens input_nml model_nml is ignored.")
        
        ## Apply patches
        init_ens_filter_nml = init_ens_filter_nml.todict()
        for kk in list(init_ens_filter_nml_update.keys()):
            if init_ens_filter_nml_update[kk] != {}:
                #_ = init_ens_filter_nml_update.pop(kk)
                init_ens_filter_nml[kk].update(init_ens_filter_nml_update[kk])        
        init_ens_filter_nml = f90nml.Namelist(init_ens_filter_nml)
        
        input_file_lists = []
        domain_order_list = []
        domain_shapefile_list = []
        
        if input_state_file_list['hydro_file_list.txt'] is not None:
            input_file_lists.append('hydro_file_list.txt')
            domain_order_list.append('hydro')
            domain_shapefile_list.append(str(input_state_file_list['hydro_file_list.txt']))

        if input_state_file_list['restart_file_list.txt'] is not None:
            input_file_lists.append('restart_file_list.txt')
            domain_order_list.append('lsm')
            domain_shapefile_list.append(str(input_state_file_list['restart_file_list.txt']))

        if input_state_file_list['param_file_list.txt'] is not None:
            if input_state_file_list['param_file_list.txt'].parent != init_ens_dir:
                _ = shutil.copy(
                    str(input_state_file_list['param_file_list.txt']),
                    init_ens_dir / input_state_file_list['param_file_list.txt'].name
                )
            input_file_lists.append('param_file_list.txt')
            domain_order_list.append('parameters')
            domain_shapefile_list.append(str(input_state_file_list['param_file_list.txt']))

        init_ens_filter_nml['filter_nml']['input_state_file_list'] = input_file_lists
        init_ens_filter_nml['model_nml']['domain_shapefiles'] = domain_shapefile_list
        # domain_order is not usedin the noah namelist
        if  input_state_file_list['hydro_file_list.txt'] is not None:
            init_ens_filter_nml['model_nml']['domain_order'] = domain_order_list
        else:
            del init_ens_filter_nml['model_nml']['domain_order']
        
        
        # Create the hydro_file_list.txt and param_file_list.txt files
        for ff in input_file_lists:
            with open(init_ens_filter_dir / ff, 'w') as opened_file:
                for mm in range(config['ensemble']['size']):
                    _ = opened_file.write(str(input_state_file_list[ff]) + '\n')

        output_file_lists = []
        if input_state_file_list['hydro_file_list.txt'] is not None:
            output_file_lists.append('hydro_file_list.out.txt')
        if input_state_file_list['restart_file_list.txt'] is not None:
            output_file_lists.append('restart_file_list.out.txt')
        if input_state_file_list['param_file_list.txt'] is not None:
            output_file_lists.append('param_file_list.out.txt')
        init_ens_filter_nml['filter_nml']['output_state_file_list'] = output_file_lists

        # Create the output_state_file_list
        # In the process copy the deterministic restart file to the ensemble
        # restart files so that all variables are present after "from_filter"
        for ff in input_file_lists:
            ff_out = ff.replace('.txt', '.out.txt')
            restart_basename = os.path.basename(str(input_state_file_list[ff]))
            with open(init_ens_filter_dir / ff_out, 'w') as opened_file:
                for mm in range(config['ensemble']['size']):
                    init_member_dir = config['initial_ens']['path'] / ('member_' + "%03d" % (mm,))
                    init_member_dir.mkdir(parents=True, exist_ok=True)
                    _ = opened_file.write(str(init_member_dir / restart_basename) + '\n')
                    shutil.copy(str(input_state_file_list[ff]), str(init_member_dir / restart_basename))
                    (init_member_dir / restart_basename).chmod(0o644)

        init_ens_filter_nml['filter_nml']['perturb_from_single_instance'] = True
        init_ens_filter_nml['filter_nml']['ens_size'] = config['ensemble']['size']

        # Write the input.nml
        init_ens_filter_nml.write(init_ens_filter_dir / 'input.nml', force=True)

        wrf_hydro_setup_0 = wrf_hydro_ens_sim.members[0]
        nml = f90nml.Namelist(wrf_hydro_setup_0.base_hydro_namelist)
        # TODO(JLM): the following lines should not be necessary:
        _ = nml['hydro_nlist'].pop('chanobs_domain')
        _ = nml['hydro_nlist'].pop('frxst_pts_out')
        _ = nml['hydro_nlist'].pop('io_config_outputs')
        _ = nml['hydro_nlist'].pop('io_form_outputs')
        nml.write(init_ens_filter_dir / 'hydro.namelist', force=True)

        # Have to create a dummy obs_seq.out file to run filter
        if input_state_file_list['hydro_file_list.txt'] is None:
            hydro_rst_deterministic = \
                config['wrf_hydro']['domain_src'] / nml['hydro_nlist']['restart_file']
        else:
            hydro_rst_deterministic = input_state_file_list['hydro_file_list.txt']
        obs_seq_dummy(
            hydro_rst=hydro_rst_deterministic,
            dir=init_ens_filter_dir
        )

        # TODO(JLM): eventually need to do this for the RESTART file namelist

        # Also need domain files. Get the entire config dir.
        top_level_dir = get_top_level_dir_from_config(config, wrf_hydro_setup_0)
        init_ens_filter_config_dir = init_ens_filter_dir / top_level_dir
        init_ens_filter_config_dir.symlink_to(
            wrf_hydro_setup_0.domain.domain_top_dir / top_level_dir
        )

        # Run filter
        if 'cmd' in config['initial_ens']['from_filter'].keys():
            filter_cmd = config['initial_ens']['from_filter']['cmd']
        else:
            filter_cmd = './filter'
            if config['dart']['mpi']:
                filter_cmd = 'mpirun --host `hostname` -np 1 ' + filter_cmd

        spr = subprocess.run(shlex.split(filter_cmd), cwd=init_ens_filter_dir)

        if spr.returncode != 0:
            raise ValueError("Filter failed to create the inital ensemble.")

        # # Advance the initial ens to get a new initial ens?
        # if not config['initial_ens']['advance']['end_time']:
        #     # symlink the filter-created restarts to config['initial_ens']['path']
        #     pass
        # else:
        #     print('Initial Ensemble: Establish WrfHydroEnsembleRun object.')
        #     initial_ens_run_dir = config['initial_ens']['path']
        #     wrf_hydro_ens_run = wrfhydropy.WrfHydroEnsembleRun(
        #         wrf_hydro_ens_sim,
        #         run_dir=initial_ens_run_dir / 'run',
        #         mode=config['initial_ens']['from_filter']['mode']
        #     )
            # TODO(JLM): Roughly...
            # 0) pickle the ens_run
            # 1) Set the restart files to the
            #    config['initial_ens']['path'] / member_iii / file.nc
            # 2) Advance.
            # 3) Then symlink the files back up a dir.

        # Protect the initial ensemble member dirs with chmod 444
        protect = config['initial_ens']['path'].glob('member_*')
        for pp in protect:
            for root, dirs, files in os.walk(pp):
                for d in dirs:
                    os.chmod(os.path.join(root, d), 0o444)
                for f in files:
                    os.chmod(os.path.join(root, f), 0o444)

    # #######################################################
    # Edit the ensemble namelist to take the initial ensemble.

    for mm in range(config['ensemble']['size']):

        rst_path = config['initial_ens']['path'] / ('member_' + "%03d" % (mm,))
        member = wrf_hydro_ens_sim.members[mm]

        forgot = "\n*** Maybe you forgot to create the initial ensemble?? ***"
        
        hydro_rst = list(rst_path.glob('HYDRO_RST*'))
        if len(hydro_rst) > 1:
            raise ValueError('Too many initial HYDRO_RST files in ' + str(rst_path))
        if len(hydro_rst) == 0:
            raise ValueError('No initial HYDRO_RST files provided in ' + str(rst_path) + forgot)
        member.base_hydro_namelist['hydro_nlist']['restart_file'] = str(hydro_rst[0])

        lsm_rst = list(rst_path.glob('RESTART*'))
        if len(lsm_rst) > 1:
            raise ValueError('Too many initial RESTART files in ' + str(rst_path))
        if len(lsm_rst) == 0:
            forc_typ = member.base_hrldas_namelist['wrf_hydro_offline']['forc_typ']
            if forc_typ not in [9, 10]:
                raise ValueError('No inital RESTART files in (and forc_typ is neither 9 nor 10):' +
                                 str(rst_path) + forgot)
        else:
            member.base_hrldas_namelist['noahlsm_offline']['restart_filename_requested'] = lsm_rst[0]

    # update the member_diffs. Would be nice if this were more automated.
    _ = wrf_hydro_ens_sim.member_diffs
    del _

    return wrf_hydro_ens_sim
