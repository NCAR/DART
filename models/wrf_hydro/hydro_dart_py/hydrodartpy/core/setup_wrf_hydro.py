import pickle
import wrfhydropy

def setup_wrf_hydro(config):

    print('wrfhydropy.Domain object.')
    wrf_hydro_domain = wrfhydropy.Domain(
        domain_top_dir=config['wrf_hydro']['domain_src'],
        domain_config=config['wrf_hydro']['model_config']
        #model_version=config['wrf_hydro']['domain_version']
    )

    print('wrfhydropy.Model object.')
    wrf_hydro_build_dir = config['experiment']['experiment_dir'] / config['wrf_hydro']['build_dir']
    wrf_hydro_model_pickle = wrf_hydro_build_dir / 'WrfHydroModel.pkl'

    if config['wrf_hydro']['use_existing_build']:

        if wrf_hydro_model_pickle.exists():
            print('WrfHydroModel object: using existing pickle.')
            wrf_hydro_model = pickle.load(open(wrf_hydro_model_pickle, 'rb'))
        else:
            raise ValueError('Existing WRF_HYDRO model requested but no pickle file found')

    else:

        # TODO(JLM): If fork specified, pull fork. Check out desired commit.
        # TODO(JLM): clone to build_dir? or make a clone_dir in this case?
        #            If not fork, only use local repo state.

        wrf_hydro_model = wrfhydropy.Model(
            source_dir=config['wrf_hydro']['wrf_hydro_src'] / 'trunk/NDHMS/' ,
            model_config=config['wrf_hydro']['model_config'],
            compiler=config['wrf_hydro']['compiler'],
            hydro_namelist_config_file=config['wrf_hydro']['hydro_namelist_config_file'],
            hrldas_namelist_config_file=config['wrf_hydro']['hrldas_namelist_config_file'],
            compile_options_config_file=config['wrf_hydro']['compile_options_config_file']
        )

        # Apply compile option patches.
        if config['wrf_hydro']['compile_options'] is not None:
            wrf_hydro_model.compile_options.update(config['wrf_hydro']['compile_options'])

        wrf_hydro_model.compile(
            compile_dir=config['experiment']['experiment_dir'] / config['wrf_hydro']['build_dir']
        )

    print('wrfhydropy.Job object')
    job = wrfhydropy.Job(
        exe_cmd='mpirun -np 4 ./wrf_hydro.exe',
        job_id='nwm_spinup',
        restart=False
    )

    
    print('wrfhydropy.Simulation object.')
    wrf_hydro_sim = wrfhydropy.Simulation()
    wrf_hydro_sim.add(wrf_hydro_model)
    wrf_hydro_sim.add(wrf_hydro_domain)

    # These were deep copied, delete them.
    del wrf_hydro_model, wrf_hydro_domain

    return wrf_hydro_sim
