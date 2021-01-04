
<span id="TOP" class="anchor"></span>

![DARTlogo](https://github.com/NCAR/DART/blob/Manhattan/docs/images/Dartboard7.png)

### Further description of the scripts

Scripts call a fortran tool to post-process the FESOM and DART outputs. They include a **tool** variable
which is inserted into the namelist.config which is read by the fortran program. Please see
descriptions in ```diagnostics/src/fesom_post_main.F90```. Other variables can be set in
```dart.postproc.env``` and ```fesom.postproc.env``` source files.

All of the scripts call visualization scripts written in GMT except ```dart_obs_seq_diag``` uses
FERRET to visualize DART diagnostic outputs.

| script name|tool code | plotting tool |description|
|------------|----------|---------------|------------|
|```compute_ensemble_mean``` | 7   |```plot_ensemble_mean.gmt```    | computes ensemble mean and extracts a transect or level|
|```compute_increment```     |13   |```plot_increment.gmt```      | computes increment using DART diagnostic output |
|```compute_NR_diff```       | 9   |```plot_NR_diff.gmt```     | computes the difference between a nature run and the ensemble prior mean |
|```observe_nature_run```    | 8,11|    | creates synthetic observations from a nature run |
|```transect_daily_mean```   | 2   |```transect_daily_mean.gmt```    | extracts and plots a transect of an individual ensemble member|
|```zlevel_daily_mean```     | 3   |```zlevel_yearly_mean.gmt```    | extracts and plots a level of an individual ensemble member|
|```dart.postproc.env```     |     |    | DART environment variables |
|```fesom.postproc.env```    |     |    | FESOM environment variables|
|```dart_obs_seq_diag```     |     | ```frt.obs_diag_TeMPLaTe.jnl``` and ```frt.obs_epoch_TeMPLaTe.jnl```   | DART observation-space statistics from ```obs_epoch.nc``` and ```obs_diag.nc```|

Similar scripts and plotting tools can be written to use the tools below.

|tool code | called routine |
|-------|---------------------|
|tool=1 | basin_mean_evolution|
|tool=2 | read_thalweg_from_nc|
|tool=3 | read_section_from_netcdf|
|tool=4 | marmara_mean_evolution|
|tool=5 | calc_section_monthly_mean|
|tool=6 | calc_thalweg_monthly_mean|
|tool=7 | read_ensemble_from_netcdf|
|tool=8 | synthetic_ferrybox_from_nr|
|tool=9 | read_section_from_NR_diff|
|tool=10| read_ctd_data|
|tool=11| profile_from_netcdf|
|tool=12| velocity_at_the_exit|
|tool=13| read_section_from_inc|
|tool=14| dardanelles_for_MFS|
|tool=15| total_kinetic_energy|
|tool=16| surface_kinetic_energy|
|tool=17| calc_section_annual_mean|
|tool=18| calc_thalweg_annual_mean|
|tool=19| compute_vorticity|
|tool=20| compute_wind_stress_curl|
|tool=21| compute_net_flux|
|tool=22| compute_surface_buoyancy|
|tool=23| compute_forcing_monthly_timeseries|
|tool=24| compute_wind_work|
|tool=25| read_ship_track|
|tool=26| bosphorus_for_blk_mfs|
|tool=27| compute_volume_transport|
