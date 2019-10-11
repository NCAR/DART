import code
import os
import sys
import wrfhydropy


def setup_wrf_hydro_ens(config, wrf_hydro_sim):

    print('wrfhydropy.EnsembleSimulation object.')

    # TODO(JLM): This is one way of making a fairly homogenous ensemble
    # might be desirable to move this to the constructor script for
    # more variations

    wrf_hydro_ens_sim = wrfhydropy.EnsembleSimulation(ncores=config['ensemble']['ncores_setup'])
    wrf_hydro_ens_sim.add([wrf_hydro_sim])
    wrf_hydro_ens_sim.replicate_member(config['ensemble']['size'])

    # Note: do not pickle the ensemble setup, only pickle the ensemble run object.

    # Set the member run dirs to have the basename of the experiment's run dir.
    #for mm in wrf_hydro_ens_sim.members:
    #    mm.run_dir = config['experiment']['run_dir'] / mm.run_dir

    # Pick up the inital ensemble after making sure it's created (below).
    # We setup the ensemble before creating the initial ensemble in case
    # we want to use it to advance the initial ensemble which is created
    # from filter (a future feature).

    # update the members_diffs. Would be nice if this were more automated.
    _ = wrf_hydro_ens_sim.member_diffs
    del _

    # ###################################
    # Ensemble construction.
    # The file config['ensemble']['setup_py'] contains python actions on the
    # WrfHydroEnsembleSetup object.
    print('wrfhydropy.EnsembleSimulation object: ensemble construction.')

    if config['ensemble']['constructor'] is None:

        banner = \
            '\n\n *** Interactive Mode *** \n' + \
            '\nBecause no ensemble:constructor file was specified, you are entering\n' + \
            'interactive mode. Here you may compose commands to construct your  \n' + \
            'desired ensemble from the wrf_hydro_ens_sim object. Please save  \n' + \
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
    
    return wrf_hydro_ens_sim
