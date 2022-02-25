import pickle
from .dartclasses import *

def setup_dart(config):

    print('DART setup object.')

    dart_build_dir = config['experiment']['experiment_dir'] / config['dart']['build_dir']
    dart_compile_pickle = dart_build_dir / 'DartCompile.pkl'

    if config['dart']['use_existing_build']:

        if dart_compile_pickle.exists():
            print('DART setup object: using existing pickle.')
            dart_compile = pickle.load(open(dart_compile_pickle, 'rb'))
        else:
            raise ValueError('Existing DART setup requested but no pickle file found')

    else:

        dart_compile = DartCompile(
            source_dir=config['dart']['dart_src'],
            mkmf_template=config['dart']['mkmf_template'],
            mpi=config['dart']['mpi'],
            build_dir=dart_build_dir,
            work_dirs=config['dart']['work_dirs']
        )

        dart_comp_keys = dart_compile.__dict__.keys()
        models_list = [
            oo for oo in list(dart_compile.__dict__.keys()) if 'model' in oo]
        if len(models_list) != 1:
            raise ValueError("Currently only supporting individual models")
        dart_model = dart_compile.__dict__[models_list[0]]
        dart_model.input_nml['filter_nml']['ens_size'] = config['ensemble']['size']
        dart_model.input_nml['filter_nml']['num_output_state_members']  = config['ensemble']['size']
        dart_model.input_nml['filter_nml']['num_output_obs_members']  = config['ensemble']['size']

        dart_model.input_nml.write(str(dart_model.input_nml_file), force=True)
        dart_compile.pickle()

