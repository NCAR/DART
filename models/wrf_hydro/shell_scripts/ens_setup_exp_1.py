# Change restart frequency to hourly
## Hydro restart
att_tuple = ('hydro_namelist', 'hydro_nlist', 'rst_dt')
values = 60
wrf_hydro_ens_setup.set_diffs_dict(att_tuple, values)

[mm.hydro_namelist['hydro_nlist']['rst_dt'] for mm in wrf_hydro_ens_setup.members]

att_tuple = ('namelist_hrldas', 'noahlsm_offline', 'restart_frequency_hours')
values = 1
wrf_hydro_ens_setup.set_diffs_dict(att_tuple, values)
[mm.namelist_hrldas['noahlsm_offline']['restart_frequency_hours'] for mm in wrf_hydro_ens_setup.members]


# Channel-only is forc_typ = 8
att_tuple = ('namelist_hrldas', 'wrf_hydro_offline', 'forc_typ')
# Only set_diffs_dict values of length 1 and length of the ensemble are allowed.
# This is a scalar -> homogenous change to the ensemble
values = 8
#wrf_hydro_ens_setup.set_diffs_dict(att_tuple, values)
#[mm.namelist_hrldas['wrf_hydro_offline']['forc_typ'] for mm in wrf_hydro_ens_setup.members]
# The above is a homogenous change, will not be seen in the diffs_dict
#pprint(wrf_hydro_ens_setup.diffs_dict)

# Point to channel-only forcing
att_tuple = ('domain', 'forcing_dir')
# This is a non-homogenous change, will show up in the diffs_dict
values = [pathlib.PosixPath('FORCING_' + mm.number) for mm in wrf_hydro_ens_setup.members]
#wrf_hydro_ens_setup.set_diffs_dict(att_tuple, values)
#pprint(wrf_hydro_ens_setup.diffs_dict)

# TODO(JLM): domain forcing_data should be calculted by setter actions on forcing_dir.

# Does the previous need a namelist change?
# Are these two coordinate when the run is setup?
