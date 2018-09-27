# Change restart frequency to hourly in hydro namelist
att_tuple = ('hydro_namelist', 'hydro_nlist', 'rst_dt')
values = 60
wrf_hydro_ens_setup.set_diffs_dict(att_tuple, values)
[mm.hydro_namelist['hydro_nlist']['rst_dt'] for mm in wrf_hydro_ens_setup.members]

# Change restart frequency to hourly in hrldas namelist
att_tuple = ('namelist_hrldas', 'noahlsm_offline', 'restart_frequency_hours')
values = 1
wrf_hydro_ens_setup.set_diffs_dict(att_tuple, values)
[mm.namelist_hrldas['noahlsm_offline']['restart_frequency_hours'] for mm in wrf_hydro_ens_setup.members]

# Change model advance to 1 hour in hrldas namelist
# This is governed by the configuration namelist setting:
# run_experiment:  time:  advance_model_hours: 


# No other differences across the ensemble, only the FORCING dir for each
# will be set at run time by the noise_model.

# We could to parameter differences here.
