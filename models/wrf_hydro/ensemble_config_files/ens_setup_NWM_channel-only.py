# The benefit of using set_member_diffs is that it can parallelized in the future.

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


# The nwm_ana default for channel-only is forc-type 10 (channel+bucket).
# We want just the channel.
att_tuple = ('base_hrldas_namelist', 'wrf_hydro_offline', 'forc_typ')
values = 9
wrf_hydro_ens_sim.set_member_diffs(att_tuple, values)


# Sixmile does not currently have a routelink with support for compound-channel.
att_tuple = ('base_hydro_namelist', 'hydro_nlist', 'compound_channel')
values = False
wrf_hydro_ens_sim.set_member_diffs(att_tuple, values)


# Change restart frequency to hourly in hydro namelist.
att_tuple = ('base_hydro_namelist', 'hydro_nlist', 'rst_dt')
values = 60
wrf_hydro_ens_sim.set_member_diffs(att_tuple, values)


# Change output frequency to hourly in hydro namelist. We are
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
