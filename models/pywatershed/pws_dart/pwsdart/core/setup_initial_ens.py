from copy import deepcopy
import f90nml
import json
import os
import pathlib
import shlex
import shutil
import subprocess
import warnings

from .obs_seq_dummy import obs_seq_dummy

def setup_initial_ens(config):

    if config['initial_ens']['from_filter']['create']:

        print('Create initial ensemble using filter.')

        init_ens_dir = pathlib.PosixPath(config['initial_ens']['path'])
        init_ens_dir.mkdir(parents=True, exist_ok=True)

        # Short-hand
        input_state_file_list = \
            config['initial_ens']['from_filter']['input_nml']['filter_nml']['input_state_file_list']

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
        
        if input_state_file_list['channel_restart_file_list.txt'] is not None:
            input_file_lists.append('channel_restart_file_list.txt')
            domain_order_list.append('channel')
            domain_shapefile_list.append(str(input_state_file_list['channel_restart_file_list.txt']))

        if input_state_file_list['hru_restart_file_list.txt'] is not None:
            input_file_lists.append('hru_restart_file_list.txt')
            domain_order_list.append('hru')
            domain_shapefile_list.append(str(input_state_file_list['hru_restart_file_list.txt']))

        """ 
            Arezoo : we do not have parameters at this point 

        if input_state_file_list['param_file_list.txt'] is not None:
            if input_state_file_list['param_file_list.txt'].parent != init_ens_dir:
                _ = shutil.copy(
                    str(input_state_file_list['param_file_list.txt']),
                    init_ens_dir / input_state_file_list['param_file_list.txt'].name
                )
            input_file_lists.append('param_file_list.txt')
            domain_order_list.append('parameters')
            domain_shapefile_list.append(str(input_state_file_list['param_file_list.txt']))
        """

        init_ens_filter_nml['filter_nml']['input_state_file_list'] = input_file_lists
        init_ens_filter_nml['model_nml']['domain_shapefiles'] = domain_shapefile_list
        init_ens_filter_nml['model_nml']['domain_order'] = domain_order_list
        
        # Create the channel_restart_file_list.txt and param_file_list.txt files
        for ff in input_file_lists:
            with open(init_ens_filter_dir / ff, 'w') as opened_file:
                for mm in range(config['ensemble']['size']):
                    _ = opened_file.write(str(input_state_file_list[ff]) + '\n')

        output_file_lists = []
        if input_state_file_list['channel_restart_file_list.txt'] is not None:
            output_file_lists.append('channel_restart_file_list.out.txt')
        if input_state_file_list['hru_restart_file_list.txt'] is not None:
            output_file_lists.append('hru_restart_file_list.out.txt')

        """ Arezoo: We do not have parameters at this time
        if input_state_file_list['param_file_list.txt'] is not None:
            output_file_lists.append('param_file_list.out.txt')"""
        init_ens_filter_nml['filter_nml']['output_state_file_list'] = output_file_lists

        # Copy the static file list from the head directory
        static_file = model_work_path / 'static_file_list.txt'
        if static_file.exists():
            shutil.copy(static_file, init_ens_filter_dir / 'static_file_list.txt')

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

        # Arezoo: I do not know why this part does not get overwritten properly, I did it manually for now. 
        init_ens_filter_nml['filter_nml']['obs_sequence_in_name'] = 'obs_seq.out'
        init_ens_filter_nml['filter_nml']['obs_sequence_out_name'] = 'obs_seq.final'

        # Write the input.nml
        init_ens_filter_nml.write(init_ens_filter_dir / 'input.nml', force=True)

        # Have to create a dummy obs_seq.out file to run filter
        if input_state_file_list['channel_restart_file_list.txt'] is None:
             print("Please define the channel_restart_file_list.txt")
             #   hydro_rst_deterministic = \
             #       config['wrf_hydro']['domain_src'] / nml['hydro_nlist']['restart_file']
        else:
            hydro_rst_deterministic = input_state_file_list['channel_restart_file_list.txt']
        obs_seq_dummy(
            hydro_rst=hydro_rst_deterministic,
            dir=init_ens_filter_dir
        )
        obs_seq_dummy

        # TODO : Arezoo: This is right now hardcoded, so fix it later please 

        # Also need domain files. Get the entire config dir.
        domain_file = config['pywatershed']['domain_dir'] + '/parameters_dis_seg_app.nc'
        init_ens_filter_config_dir = init_ens_filter_dir / 'parameters_dis_seg_app.nc'
        print(init_ens_filter_dir)
        init_ens_filter_config_dir.symlink_to(
            domain_file
        )

        domain_file = config['pywatershed']['domain_dir'] + '/parameters_PRMSChannel.nc'
        init_ens_filter_config_dir = init_ens_filter_dir / 'parameters_PRMSChannel.nc'
        print(init_ens_filter_dir)
        init_ens_filter_config_dir.symlink_to(
            domain_file
        )

        domain_file = config['pywatershed']['domain_dir'] + '/parameters_dis_hru.nc'
        init_ens_filter_config_dir = init_ens_filter_dir / 'parameters_dis_hru.nc'
        print(init_ens_filter_dir)
        init_ens_filter_config_dir.symlink_to(
            domain_file
        )

        # Run filter
        if 'cmd' in config['initial_ens']['from_filter'].keys():
            filter_cmd = config['initial_ens']['from_filter']['cmd']
        else:
            filter_cmd = './filter'
            if config['dart']['mpi']:
                filter_cmd = 'mpirun -np 1 ' + filter_cmd

        spr = subprocess.run(shlex.split(filter_cmd), cwd=init_ens_filter_dir)

        if spr.returncode != 0:
            raise ValueError("Filter failed to create the inital ensemble.")

        # Protect the initial ensemble member dirs with chmod 444
        protect = config['initial_ens']['path'].glob('member_*')
        for pp in protect:
            for root, dirs, files in os.walk(pp):
                for d in dirs:
                    os.chmod(os.path.join(root, d), 0o444)
                for f in files:
                    os.chmod(os.path.join(root, f), 0o444)

