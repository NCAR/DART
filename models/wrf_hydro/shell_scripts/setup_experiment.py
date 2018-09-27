# Setup a "HydroDartRun". For now this is a filter experiment.

# Notes:
# * The approach to ensemble construction uses a single "setup" object as
#   basis for the ensemble. This is probably fine for most cases but is
#   not totally general.
#

import argparse
from boltons.iterutils import remap
import code
import copy
from create_parameter_restart_file import create_parameter_restart_file
import datetime
import f90nml
from obs_seq_dummy import obs_seq_dummy
import os
import pathlib
import pickle
from pprint import pprint
import re
import shlex
import shutil
import subprocess
import sys
import warnings
import wrfhydropy

from setup_experiment_tools import establish_config

# ###################################
# Args
# Arguments to this script.

# python setup_experiment.py --help

parser = argparse.ArgumentParser(description='Setup a WRF-Hydro-DART experiment.')
# Single positional argument.
parser.add_argument(
    'config_file',
    metavar='config_file.yaml',
    type=str,
    nargs=1,
    help='The YAML experiment configuration file of arbitrary name.'
)
args = parser.parse_args()

# #################################
# Read the yaml config file
config = establish_config(args.config_file[0])


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
    'noise_function'
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

print("The configuration file:")
pprint(config); print('')


# ###################################
# Setup experiment_dir
print('Create experiment directory.')
config['experiment']['experiment_dir'].mkdir(
    parents=True,
    exist_ok = config['experiment']['experiment_dir_mode'] == 'w'
)
# TODO(JLM): What to do if the dir exists? Use python mode 'w' for clobber?
# Setting up and experiment should be a totally automated process,
# so if anything needs to change (in the yaml), then just resetup and rerun.

# Copy the config file to the experiment dir
# TODO(JLM): Not sure the best approach here: Rename this copied file? This file does not match
# the internal representation, it reflects state of the original setup. One approach:
# copy it here and rename as a dot file: ".original.<file_name>" then stash that
# path in the internal config file representation. Maybe also stash the path to the original file?
file_dir = os.path.dirname(args.config_file[0])
file_base = os.path.basename(args.config_file[0])
config_file_copy_name = file_dir + '.original.' + file_base 
shutil.copy(
    args.config_file[0],
    config['experiment']['experiment_dir'] / os.path.basename(config_file_copy_name)
)


# ###################################
# DART Setup
print('DART setup object')

dart_build_dir = config['experiment']['experiment_dir'] / config['dart']['build_dir']
dart_setup_pickle = dart_build_dir / 'DartSetup.pkl'

if config['dart']['use_existing_build']:

    if dart_setup_pickle.exists():
        print('DART setup object: using existing pickle.')
        dart_setup = pickle.load(open(dart_setup_pickle, 'rb'))
    else:
        raise ValueError('Existing DART setup requested but no pickle file found')

else:

    # todo(JLM): If fork specified, pull fork. Check out desired commit.
    # TODO(JLM): clone to build_dir? or make a clone_dir in this case?
    #            If not fork, only use local repo state.
    dart_setup = wrfhydropy.DartSetup(
        source_dir=config['dart']['dart_src'],
        mkmf_template=config['dart']['mkmf_template'],
        mpi=config['dart']['mpi'],
        build_dir=dart_build_dir
    )

    dart_setup.input_nml['filter_nml']['ens_size'] = config['ensemble']['size']
    dart_setup.input_nml['filter_nml']['num_output_state_members']  = config['ensemble']['size']
    dart_setup.input_nml['filter_nml']['num_output_obs_members']  = config['ensemble']['size']

    dart_setup.input_nml.write(str(dart_setup.input_nml_file), force=True)
    dart_setup.pickle()


# ###################################
# WRF-Hydro deterministic
print('WrfHydroDomain object.')

wrf_hydro_domain = wrfhydropy.WrfHydroDomain(
    domain_top_dir=config['wrf_hydro']['domain_src'],
    domain_config=config['wrf_hydro']['model_config'],
    model_version=config['wrf_hydro']['domain_version']
)

print('WrfHydroModel object.')

wrf_hydro_build_dir = config['experiment']['experiment_dir'] / config['wrf_hydro']['build_dir']
wrf_hydro_model_pickle = wrf_hydro_build_dir / 'WrfHydroModel.pkl'

if config['wrf_hydro']['use_existing_build']:

    if wrf_hydro_model_pickle.exists():
        print('WrfHydroModel object: using existing pickle.')
        wrf_hydro_model = pickle.load(open(wrf_hydro_model_pickle, 'rb'))
    else:
        raise ValueError('Existing WRF_HYDRO model requested but no pickle file found')

else: 

    # todo(JLM): If fork specified, pull fork. Check out desired commit.
    # TODO(JLM): clone to build_dir? or make a clone_dir in this case?
    #            If not fork, only use local repo state.
    wrf_hydro_model = wrfhydropy.WrfHydroModel(
        source_dir=config['wrf_hydro']['wrf_hydro_src'] / 'trunk/NDHMS/' ,
        model_config=config['wrf_hydro']['model_config']
    )
    wrf_hydro_model.compile(
        compile_dir=config['experiment']['experiment_dir'] / config['wrf_hydro']['build_dir'],
        compiler=config['wrf_hydro']['compiler']
    )

# TODO(JLM): defencies above:
# compile should have a commit attribute. consider writing this to file in the compile dir.
# would also like to copy the compile options to the compile dir.

print('WrfHydroSetup object.')

wrf_hydro_setup = wrfhydropy.WrfHydroSetup(
    wrf_hydro_model=wrf_hydro_model,
    wrf_hydro_domain=wrf_hydro_domain
)

# These were deep copied, delet them.
del wrf_hydro_model, wrf_hydro_domain

# WRF-Hydro Ensemble Setup
print('WrfHydroEnsSetup object.')

# TODO(JLM): This is one way of making a fairly homogenous ensemble
# might be desirable to move this to the constructor script for
# more variations

wrf_hydro_ens_setup = wrfhydropy.WrfHydroEnsembleSetup([wrf_hydro_setup])
wrf_hydro_ens_setup.replicate_member(config['ensemble']['size'])

# Note: do not pickle the ensemble setup, only pickle the ensemble run object.

# Set the member run dirs to have the basename of the experiment's run dir.
for mm in wrf_hydro_ens_setup.members:
    mm.run_dir = config['experiment']['run_dir'] / mm.run_dir

# Pick up the inital ensemble after making sure it's created (below).
# We setup the ensemble before creating the initial ensemble in case
# we want to use it to advance the initial ensemble which is created
# from filter (a future feature).

# update the diffs_dict. Would be nice if this were more automated.
_ = wrf_hydro_ens_setup.diffs_dict
del _


# ###################################
# Ensemble construction.
# The file config['ensemble']['setup_py'] contains python actions on the
# WrfHydroEnsembleSetup object.
print('WrfHydroEnsSetup object: ensemble construction.')

if config['ensemble']['constructor'] is None:

    banner = \
        '\n\n *** Interactive Mode *** \n' + \
        '\nBecause no ensemble:constructor file was specified, you are entering\n' + \
        'interactive mode. Here you may compose commands to construct your  \n' + \
        'desired ensemble from the wrf_hydro_ens_setup object. Please save  \n' + \
        'comands to a file you will later specify as ensemble:constructor. Your\n' + \
        'commands are lost upon exiting interactive mode and the process of \n' + \
        'constructing the enemble only completes from a specified file.     \n' + \
        '\n' + \
        'You can restart the experiment setup proces as many times as you   \n' + \
        'When restarting, you can avoid recompiling both DART and WRF-Hydro \n' + \
        'by setting use_existing_build to True for both. This will quickly  \n' + \
        'return you to this ensemble construction step where you can once   \n' + \
        'again edit the "vanilla" ensemble or proceed through the specified \n' + \
        'script. \n\n'

    import readline
    import rlcompleter
    readline.parse_and_bind("tab: complete")
    code.interact(local=locals(), banner=banner)
    sys.exit()

else:

    exec(open(config['ensemble']['constructor']).read())


# ###################################
# Initial ensemble work/setup
# TODO (JLM): make this a separate script.
# This leverages the same setup object for advancing the model.

# This section if either create param file or create ensemble from filter.
if config['initial_ens']['param_restart_file']['create'] or \
   config['initial_ens']['from_filter']['create']:

    init_ens_dir = pathlib.PosixPath(config['initial_ens']['path'])
    init_ens_dir.mkdir(parents=True, exist_ok=True)

    # Short-hand
    input_state_file_list = \
        config['initial_ens']['from_filter']['input_nml']['filter_nml']['input_state_file_list']

    # TODO(JLM): Dont allow clobber? Could wait until write, but maybe favor preemption.
    domain_config_rst_path = \
        config['wrf_hydro']['domain_src'] / (config['wrf_hydro']['model_config'] + '/RESTART/')

    if input_state_file_list['hydro_file_list.txt'] is None:
        input_state_file_list['hydro_file_list.txt'] = \
            sorted(domain_config_rst_path.glob('*HYDRO_RST*'))[0]

    if input_state_file_list['param_file_list.txt'] is None:

        # First check if the file comes with the domain.
        input_state_file_list['param_file_list.txt'] = \
            list(domain_config_rst_path.glob('*param*.nc'))

        if len(input_state_file_list['param_file_list.txt']):
            warnings.warn("Using first parameter restart file found with the model domain.")
            input_state_file_list['param_file_list.txt'] = \
                input_state_file_list['param_file_list.txt'][0]
        else:
            warnings.warn("Setting parameter restart file to " +
                          "config['initial_ens']['param_restart_file'].")
            input_state_file_list['param_file_list.txt'] = \
                pathlib.PosixPath(config['initial_ens']['param_restart_file']['out_file'])


if config['initial_ens']['param_restart_file']['create']:

    print('Create parameter restart files.')
    # Create parameter "restart" file from HYDRO_RST or RESTART files.

    if input_state_file_list['param_file_list.txt'] is not None and \
       not input_state_file_list['param_file_list.txt'].exists():
        create_param_file=config['initial_ens']['param_restart_file']
        create_parameter_restart_file(
            out_file=create_param_file['out_file'],
            out_mode=create_param_file['mode'],
            hydro_rst_file=input_state_file_list['hydro_file_list.txt'],
            restart_file=None,
            existing_variables=create_param_file['existing_variables'],
            new_variables=create_param_file['new_variables'],
            values=create_param_file['values']
        )


# Initial ensemble creation (optional) from filter
if config['initial_ens']['from_filter']['create']:

    print('Create initial ensemble using filter.')

    # make a place to filter in the init_ens_filter_dir
    init_ens_filter_dir = pathlib.PosixPath(config['initial_ens']['path'] / 'from_filter')
    init_ens_filter_dir.mkdir()

    # Copy filter and input_nml
    shutil.copy(
        config['experiment']['experiment_dir'] / config['dart']['build_dir'] / 'filter' ,
        init_ens_filter_dir / 'filter'
    )
    shutil.copy(
        config['experiment']['experiment_dir'] / config['dart']['build_dir'] / 'input.nml' ,
        init_ens_filter_dir / 'input.nml'
    )

    # Bring in input_nml, apply patches.
    init_ens_filter_nml = f90nml.read(init_ens_filter_dir / 'input.nml' )

    input_file_lists = []
    if input_state_file_list['hydro_file_list.txt'] is not None:
        input_file_lists.append('hydro_file_list.txt')
    if input_state_file_list['restart_file_list.txt'] is not None:
        input_file_lists.append('restart_file_list.txt')
    if input_state_file_list['param_file_list.txt'] is not None:
        input_file_lists.append('param_file_list.txt')
    init_ens_filter_nml['filter_nml']['input_state_file_list'] = input_file_lists

    domain_shapefiles = [ str(input_state_file_list[ii]) for ii in input_file_lists]
    init_ens_filter_nml['model_nml']['domain_shapefiles'] = domain_shapefiles

    # Create the hydro_file_list.txt and param_file_list.txt files
    for ff in input_file_lists:
        with open(init_ens_filter_dir / ff, 'w') as opened_file:
            for mm in range(config['ensemble']['size']):
                _ = opened_file.write(str(input_state_file_list[ff])+'\n')

    output_file_lists = []
    if input_state_file_list['hydro_file_list.txt'] is not None:
        output_file_lists.append('hydro_file_list.out.txt')
    if input_state_file_list['restart_file_list.txt'] is not None:
        output_file_lists.append('restart_file_list.out.txt')
    if input_state_file_list['param_file_list.txt'] is not None:
        output_file_lists.append('param_file_list.out.txt')
    init_ens_filter_nml['filter_nml']['output_state_file_list'] = output_file_lists

    # Create the output_state_file_list make these all lke
    for ff in input_file_lists:
        ff_out = ff.replace('.txt','.out.txt')
        restart_basename=os.path.basename(str(input_state_file_list[ff]))
        with open(init_ens_filter_dir / ff_out, 'w') as opened_file:
            for mm in range(config['ensemble']['size']):
                init_member_dir = config['initial_ens']['path'] / ('member_' + "%03d" % (mm,))
                init_member_dir.mkdir(parents=True, exist_ok=True)
                _ = opened_file.write(str(init_member_dir / restart_basename) + '\n')

    init_ens_filter_nml['filter_nml']['perturb_from_single_instance'] = True
    init_ens_filter_nml['filter_nml']['ens_size'] = config['ensemble']['size']
    init_ens_filter_nml['model_nml'].update(
        config['initial_ens']['from_filter']['input_nml']['model_nml']
    )
    # Write the input.nml
    init_ens_filter_nml.write(init_ens_filter_dir / 'input.nml', force=True)


    # TODO(JLM): model_nml: domain_shapefiles is not a great term.

    nml = f90nml.Namelist(wrf_hydro_setup.hydro_namelist)
    # TODO(JLM): the following lines should not be necessary:
    _ = nml['hydro_nlist'].pop('chanobs_domain')
    _ = nml['hydro_nlist'].pop('frxst_pts_out')
    _ = nml['hydro_nlist'].pop('io_config_outputs')
    _ = nml['hydro_nlist'].pop('io_form_outputs')
    nml.write(init_ens_filter_dir / 'hydro.namelist', force=True)

    # Have to create a dummy obs_seq.out file to run filter
    obs_seq_dummy(
        hydro_rst=config['wrf_hydro']['domain_src'] / nml['hydro_nlist']['restart_file'],
        dir=init_ens_filter_dir
    )

    # TODO(JLM): eventually need to do this for the RESTART file namelist

    # Also need domain files. Get the entire config dir.
    init_ens_filter_config_dir = init_ens_filter_dir / 'NWM'
    init_ens_filter_config_dir.symlink_to(wrf_hydro_setup.domain.domain_top_dir / 'NWM')

    # Run filter
    filter_cmd = './filter'
    if config['dart']['mpi']:
        filter_cmd = 'mpirun -np 1 ' + filter_cmd
    subprocess.run(
        shlex.split(filter_cmd),
        cwd=init_ens_filter_dir
    )

    if not config['initial_ens']['advance']['end_time']:
        # symlink the filter-created restarts to config['initial_ens']['path']
        pass
    else:
        print('Initial Ensemble: Establish WrfHydroEnsembleRun object.')
        initial_ens_run_dir = config['initial_ens']['path']
        wrf_hydro_ens_run = wrfhydropy.WrfHydroEnsembleRun(
            wrf_hydro_ens_setup,
            run_dir=initial_ens_run_dir / 'run',
            mode=config['initial_ens']['from_filter']['mode']
        )

        # Set the restart files to the
        # config['initial_ens']['path'] / member_iii / file.nc
        # Advance.
        # Then symlink the files back up a dir.

# #######################################################
# Edit the ensemble namelist to take the initial ensemble.

for mm in range(config['ensemble']['size']):

    rst_path = config['initial_ens']['path'] / ('member_' + "%03d" % (mm,))
    member = wrf_hydro_ens_setup.members[mm]

    hydro_rst = list(rst_path.glob('HYDRO_RST*'))
    if len(hydro_rst) > 1:
        raise ValueError('Too many initial HYDRO_RST files in ' + str(rst_path))
    member.hydro_namelist['hydro_nlist']['restart_file'] = hydro_rst[0]

    lsm_rst = list(rst_path.glob('RESTART*'))
    if len(lsm_rst) > 1:
        raise ValueError('Too many initial RESTART files in ' + str(rst_path))
    if len(lsm_rst) == 0:
        forc_typ = member.namelist_hrldas['wrf_hydro_offline']['forc_typ']
        if forc_typ not in [9,10]:
            raise ValueError('No RESTART file present (and forc_typ is neither 9 nor 10).')
    else:
        member.namelist_hrldas['noahlsm_offline']['restart_filename_requested'] = lsm_rst[0]


# update the diffs_dict. Would be nice if this were more automated.
_ = wrf_hydro_ens_setup.diffs_dict
del _


# ###################################
# Ensemble Run: establish the run directory.
print('WrfHydroEnsembleRun object.')
wrf_hydro_ens_run = wrfhydropy.WrfHydroEnsembleRun(
    wrf_hydro_ens_setup,
    run_dir=config['experiment']['run_dir']
#    mode=config['experiment']['run_dir_mode']
)

#This was deep copied into the ensemble run object.
del wrf_hydro_ens_setup
#dart_setup_pickle.unlink() ## TODO(JLM): Temporary comment only for testing.

# Now that the run_dir is established, point the experiment and run dirs to each other.
link_to_run_dir = config['experiment']['experiment_dir'].joinpath('run_dir')
if link_to_run_dir.exists():
    link_to_run_dir.unlink()
link_to_run_dir.symlink_to(config['experiment']['run_dir'])

link_to_exp_dir = config['experiment']['run_dir'].joinpath('experiment_dir')
if link_to_exp_dir.exists():
    link_to_exp_dir.unlink()
link_to_exp_dir.symlink_to(config['experiment']['experiment_dir'])


# #######################################################
# Stage the ensemble to take the initial ensemble.

for mm in range(config['ensemble']['size']):

    rst_path = config['initial_ens']['path'] / ('member_' + "%03d" % (mm,))
    member = wrf_hydro_ens_run.ens_setup.members[mm]
    run_path = config['experiment']['run_dir'] / ('member_' + "%03d" % (mm,))

    hydro_rst_old = list(run_path.glob('HYDRO_RST*'))
    _ = [rr.unlink() for rr in hydro_rst_old]
    hydro_rst = list(rst_path.glob('HYDRO_RST*'))
    # Copy these incase they get clobbered
    shutil.copy(str(hydro_rst[0]), str(run_path / hydro_rst[0].name))

    lsm_rst_old = list(run_path.glob("RESTART.*"))
    _ = [rr.unlink() for rr in lsm_rst_old]
    lsm_rst = list(rst_path.glob('RESTART*'))
    if len(lsm_rst) == 0:
        forc_typ = member.namelist_hrldas['wrf_hydro_offline']['forc_typ']
        if forc_typ not in [9,10]:
            raise ValueError('No RESTART file present (and forc_typ is neither 9 nor 10).')
    else:
        shutil.copy(str(lsm_rst), str(run_path / lsm_rst[0].name))

    param_rst = list(rst_path.glob('parameters*nc'))[0]
    shutil.copy(str(param_rst), str(run_path / param_rst.name))


# ###################################
# HydroDartRun Object
print('HydroDartRun object.')
# TODO(JLM): filter/run/etc?? I'd be more in favor of calling it a run and specifying
# the executable task to the job. But then we are mixing types of jobs.
# Perhaps the different executables have different job lists?
hydro_dart_run = wrfhydropy.HydroDartRun(
    dart_setup=dart_setup,
    wrf_hydro_ens_run=wrf_hydro_ens_run,
    config=config
)

hydro_dart_run.pickle()

#del dart_setup, wrf_hydro_ens_run
#wrf_hydro_model_pickle.unlink()  ## TODO(JLM): Temporary comment only for testing.

# ###################################
# Place scripts into the run dir.
# Should it be checke that this file (setup_experiment.py) is in
# wrf_hydro_dart/models/wrf_hydro/shell_scripts?
print("Staging scripts.")

# Various scripts (tied to run_filter)
script_list = ['advance_ensemble.py', 'get_ensemble_time.py', 'set_obs_seq_times.py']
for ss in script_list:
    adv_ens_script_src = config['dart']['dart_src'] / ('models/wrf_hydro/shell_scripts/' + ss)
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

dart_exe_names = config['run_experiment']['dart']['exes']
for dd in dart_exe_names:
    dart_exe = config['experiment']['experiment_dir'] / \
               (config['dart']['build_dir'] + '/' + dd)
    dart_exe_link = config['experiment']['run_dir'] / dd
    dart_exe_link.symlink_to(dart_exe)

input_nml_patches = config['run_experiment']['dart']['input_nml_patches']
dart_input_nml = dart_exe = config['experiment']['experiment_dir'] / \
           (config['dart']['build_dir'] + '/input.nml')
dart_input_nml_copy = config['experiment']['run_dir'] / 'input.nml'

# write model namelists to the run dir (needed by dart in some generic way -- why exactly?
# do the namelist details matter?
model_nlsts = [['hydro_namelist', 'namelist_hrldas'][0]] # 0 for testing if 1 is necessary
for nn in model_nlsts:
    nlst = wrf_hydro_setup.__dict__[nn]
    for k0 in nlst.keys():
        for k1 in nlst[k0].keys():
            if type(nlst[k0][k1]) is str:
                nlst[k0][k1] = nlst[k0][k1].replace('./','./member_000/')
    f90nml.namelist.Namelist(nlst).write(config['experiment']['run_dir'] / nn.replace('_','.'))


# input.nml patches and checks
print('Apply namelist patches and checks.')
input_nml = f90nml.read(dart_input_nml)

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
    for k1 in input_nml.keys():
        if k1 in ['ens_size', 'num_output_obs_members']:
            input_nml[k0][k1] = int(config['ensemble']['size'])

init_shapefiles = [str(ff) for ff in list((config['initial_ens']['path'] / 'member_000').glob('*'))]
input_nml['model_nml']['domain_shapefiles'] = init_shapefiles
# TODO(JLM): revisit when there's some kind of order.... 

input_nml.write(dart_input_nml_copy)    


# Setup run_filter.csh_template
# Copy to run_dir then preprocess the template.
run_filter_template = \
    config['dart']['dart_src'] / 'models/wrf_hydro/shell_scripts/run_filter.csh.template'
run_filter_specific = config['experiment']['run_dir'] / 'run_filter.csh'
shutil.copy(str(run_filter_template), str(run_filter_specific))

def replace_in_file(file, str_find, str_replace):
    with open(file, "r") as sources:
        lines = sources.readlines()
    with open(file, "w") as sources:
        for line in lines:
            _ = sources.write(re.sub(str_find, str_replace, line))

replace_in_file(
    run_filter_specific,
    'EXPERIMENT_DIRECTORY',
    str(config['experiment']['run_dir'])
)

replace_in_file(
    run_filter_specific,
    'OBSERVATION_DIR',
    str(config['observation_preparation']['all_obs_dir'])
)

advance_cmd = "python advance_ensemble.py"
config_perturb_forcing = config['run_experiment']['perturb_forcing']
if config_perturb_forcing['perturb']:
    advance_cmd += " --job_entry_cmd '" + config_perturb_forcing['noise_cmd'] + "'"

replace_in_file(
    run_filter_specific,
    'PYTHON_ADVANCE',
    advance_cmd
)

# Interactive use if desired.
import readline
import rlcompleter
readline.parse_and_bind
code.interact(local=locals())

sys.exit()
