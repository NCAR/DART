import pathlib

# Must turn off nudging in the namelist, or the requested files will be validated.
# There's a wrfhydropy issue about this being a pain.
nudge_keys = wrf_hydro_ens_sim.members[0].base_hydro_namelist['nudging_nlist'].keys()
for kk in nudge_keys:
    att_tuple = ('base_hydro_namelist', 'nudging_nlist', kk)
    values = ''
    wrf_hydro_ens_sim.set_member_diffs(att_tuple, values)

    
# Set the perturbed forcing dir.
# Setting a directory, like this, which does not exist at simulation compose time
# (when the ensemble is established on disk) requires the argument:
# wrfhydropy.Simulaton.compose(check_nlst_warn=True)
att_tuple = ('base_hrldas_namelist', 'noahlsm_offline', 'indir')
values = './FORCING_perturbed'
wrf_hydro_ens_sim.set_member_diffs(att_tuple, values)


# Forc-type 10 (channel+bucket).
att_tuple = ('base_hrldas_namelist', 'wrf_hydro_offline', 'forc_typ')
values = 10
wrf_hydro_ens_sim.set_member_diffs(att_tuple, values)


# This should be the default, but I will leave it here.
att_tuple = ('base_hydro_namelist', 'hydro_nlist', 'compound_channel')
values = True
wrf_hydro_ens_sim.set_member_diffs(att_tuple, values)


# Change restart frequency to hourly in hydro namelist.
att_tuple = ('base_hydro_namelist', 'hydro_nlist', 'rst_dt')
values = 60
wrf_hydro_ens_sim.set_member_diffs(att_tuple, values)


# Change output frequency. We are
# currently only interested in prognostic variables and are doing hourly restarts.
# Must turn off channel_bucket_influx to do this.
att_tuple = ('base_hydro_namelist', 'hydro_nlist', 'out_dt')
values = 999999999
wrf_hydro_ens_sim.set_member_diffs(att_tuple, values)

att_tuple = ('base_hydro_namelist', 'hydro_nlist', 'output_channelbucket_influx')
values = 0
wrf_hydro_ens_sim.set_member_diffs(att_tuple, values)

# Change model advance to 1 hour in hrldas namelist
# This is governed by the configuration namelist setting:
# run_experiment:  time:  advance_model_hours: 


# =============================================================================
# Multiphysics
# -----------------------------------------------------------------------------
# The run dir does not yet exist. Put the ensemble domain files into the experiment directory.
domain_dir = wrf_hydro_ens_sim.members[0].domain.domain_top_dir
src_dir = domain_dir / 'NWM/DOMAIN/ensemble_param_files/jtti_streamflow_da_multiphysics_channel6_bucket2-low-var'
ens_domain_dir = config['experiment']['experiment_dir'] / 'ensemble_domain'

import os
import shutil

# Copy them for posterity into the experiment
if ens_domain_dir.exists():
    os.chmod(ens_domain_dir, 0o777)
    [pp.chmod(0o777) for pp in ens_domain_dir.glob('*')]
    shutil.rmtree(str(ens_domain_dir))

shutil.copytree(str(src_dir), str(ens_domain_dir))

# Mannings N
# Create a set of routelink files for the ensemble
# Code is in this dir
rl_ens_files = sorted(ens_domain_dir.glob('Route_Link_edit*'))
rl_ens_files = [str(ff) for ff in rl_ens_files]
att_tuple = ('base_hydro_namelist', 'hydro_nlist', 'route_link_f')
values = rl_ens_files
wrf_hydro_ens_sim.set_member_diffs(att_tuple, values)

# GWBUCK PARM
att_tuple = ('base_hydro_namelist', 'hydro_nlist', 'gwbuckparm_file')
values = '/glade/work/jamesmcc/domains/private/florence_933020089/NWM/DOMAIN/GWBUCKPARM_uniform_v2.1_stats.nc'
wrf_hydro_ens_sim.set_member_diffs(att_tuple, values)
