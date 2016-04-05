#!/bin/ksh -x
#########################################################################
#
# Purpose: Script to create DART/WRF Namelist 
#
#########################################################################
# Exp 1
#   inf_damping                 = 1.00,                   1.00,
#   inf_sd_initial              = 0.00,                   0.00,
#   inf_sd_lower_bound          = 0.00,                   0.00,
# Exp 2
#   inf_damping                 = 1.00,                   1.00,
#   inf_sd_initial              = 0.00,                   0.00,
#   inf_sd_lower_bound          = 0.00,                   0.00,
# Exp 3
#   inf_damping                 = 0.90,                   1.00,
#   inf_sd_initial              = 0.60,                   0.00,
#   inf_sd_lower_bound          = 0.60,                   0.00,
# Exp 4
#   inf_damping                 = 0.90,                   1.00,
#   inf_sd_initial              = 0.60,                   0.50,
#   inf_sd_lower_bound          = 0.60,                   0.10,
#   inf_initial                 = 1.00,                   1.00,
#   inf_lower_bound             = 1.0,                    1.0,
#   inf_upper_bound             = 1000000.0,              1000000.0
#
# CREATE DART/WRF NAMELIST FILE
rm -f input.nml
touch input.nml
cat > input.nml << EOF
 &perfect_model_obs_nml
   start_from_restart    = .true.,
   output_restart        = .true.,
   async                 = 2,
   init_time_days        = -1,
   init_time_seconds     = -1,
   first_obs_days           = 148816,
   first_obs_seconds        = 75601,
   last_obs_days            = 148817,
   last_obs_seconds         = 10800,
   output_interval       = 1,
   restart_in_file_name  = "perfect_ics",
   restart_out_file_name = "perfect_restart",
   obs_seq_in_file_name  = "obs_seq.in",
   obs_seq_out_file_name = "obs_seq.out",
   adv_ens_command       = "./advance_model.csh",
   output_timestamps     = .false.,
   trace_execution       = .false.,
   output_forward_op_errors = .false.,
   print_every_nth_obs   = -1,
   silence               = .false.,
/
 &filter_nml
   async                    = 2,
   adv_ens_command          = "./advance_model.csh",
   ens_size                 = ${NL_ENS_SIZE},
   start_from_restart       = .true.,
   output_restart           = ${NL_OUTPUT_RESTART},
   obs_sequence_in_name     = "obs_seq.out",
   obs_sequence_out_name    = "obs_seq.final",
   restart_in_file_name     = "filter_ic_prior",
   restart_out_file_name    = "filter_ic_post",
   init_time_days           = -1,
   init_time_seconds        = -1,
   first_obs_days           = ${NL_FIRST_OBS_DAYS},
   first_obs_seconds        = ${NL_FIRST_OBS_SECONDS},
   last_obs_days            = ${NL_LAST_OBS_DAYS},
   last_obs_seconds         = ${NL_LAST_OBS_SECONDS},
   num_output_state_members = 0,
   num_output_obs_members   = 0,
   output_interval          = 1,
   num_groups               = 1,
   input_qc_threshold       = 4.0,
   outlier_threshold        = 3.0,
   output_forward_op_errors = .false.,
   output_timestamps        = .true.,
   output_inflation         = .true.,
   trace_execution          = .true.,
   inf_flavor                  = ${NL_INF_FLAVOR_PRIOR}, ${NL_INF_FLAVOR_POST},
   inf_initial_from_restart    = ${NL_INF_INITIAL_FROM_RESTART_PRIOR}, ${NL_INF_INITIAL_FROM_RESTART_POST},
   inf_sd_initial_from_restart = ${NL_INF_SD_INITIAL_FROM_RESTART_PRIOR}, ${NL_INF_SD_INITIAL_FROM_RESTART_POST},
   inf_output_restart          = .true.,                 .true.,
   inf_deterministic           = .true.,                 .true.,
   inf_in_file_name            = 'prior_inf_ic_old',     'post_inf_ic_old',
   inf_out_file_name           = 'prior_inf_ic_new',     'post_inf_ic_new',
   inf_diag_file_name          = 'prior_inf_diag',       'post_inf_diag',
   inf_initial                 = 1.00,                   1.00,
   inf_sd_initial              = 0.60,                   0.0,
   inf_damping                 = 0.90,                   1.00,
   inf_lower_bound             = 1.0,                    1.0,
   inf_upper_bound             = 1000000.0,              1000000.0
   inf_sd_lower_bound          = 0.60,                   0.0,
/
 &smoother_nml
   num_lags              = 0
   start_from_restart    = .false.
   output_restart        = .false.
   restart_in_file_name  = 'smoother_ics'
   restart_out_file_name = 'smoother_restart',
/
 &assim_tools_nml
   filter_kind                     = 1,
   cutoff                          = 0.1,
   sort_obs_inc                    = .false.,
   spread_restoration              = .false.,
   sampling_error_correction       = .false.,
   adaptive_localization_threshold = -1,
   print_every_nth_obs             = 1000,
   special_localization_obs_types  = ${NL_SPECIAL_LOCALIZATION_OBS_TYPES},
   special_localization_cutoffs    = ${NL_SPECIAL_LOCALIZATION_CUTOFFS},
/
 &cov_cutoff_nml
   select_localization             = 1,
/
 &obs_sequence_nml
   write_binary_obs_sequence   = .false.,
/
 &preprocess_nml
   input_obs_kind_mod_file  = '${DART_DIR}/obs_kind/DEFAULT_obs_kind_mod.F90',
   output_obs_kind_mod_file = '${DART_DIR}/obs_kind/obs_kind_mod.f90',
   input_obs_def_mod_file   = '${DART_DIR}/obs_def/DEFAULT_obs_def_mod.F90',
   output_obs_def_mod_file  = '${DART_DIR}/obs_def/obs_def_mod.f90',
   input_files              = '${DART_DIR}/obs_def/obs_def_reanalysis_bufr_mod.f90',
                              '${DART_DIR}/obs_def/obs_def_radar_mod.f90',
                              '${DART_DIR}/obs_def/obs_def_metar_mod.f90',
                              '${DART_DIR}/obs_def/obs_def_dew_point_mod.f90',
                              '${DART_DIR}/obs_def/obs_def_altimeter_mod.f90',
                              '${DART_DIR}/obs_def/obs_def_gps_mod.f90',
                              '${DART_DIR}/obs_def/obs_def_gts_mod.f90',
                              '${DART_DIR}/obs_def/obs_def_vortex_mod.f90',
                              '${DART_DIR}/obs_def/obs_def_IASI_O3_mod.f90',
/  
 &obs_kind_nml
   assimilate_these_obs_types = 'RADIOSONDE_TEMPERATURE',
                                'RADIOSONDE_U_WIND_COMPONENT',
                                'RADIOSONDE_V_WIND_COMPONENT',
                                'ACARS_U_WIND_COMPONENT',
                                'ACARS_V_WIND_COMPONENT',
                                'ACARS_TEMPERATURE',
                                'AIRCRAFT_U_WIND_COMPONENT',
                                'AIRCRAFT_V_WIND_COMPONENT',
                                'AIRCRAFT_TEMPERATURE',
                                'SAT_U_WIND_COMPONENT',
                                'SAT_V_WIND_COMPONENT',
                                'IASI_O3_RETRIEVAL',
/
 &wrf_obs_preproc_nml
   obs_boundary             = 0.0,
   increase_bdy_error       = .false.,
   maxobsfac                = 2.5,
   obsdistbdy               = 15.0,
   sfc_elevation_check      = .false.,
   sfc_elevation_tol        = 300.0,
   obs_pressure_top         = 0.0,
   obs_height_top           = 200000000000000.0,
   include_sig_data         = .true.,
   tc_sonde_radii           = -1.0,
   superob_aircraft         = .false.,
   aircraft_horiz_int       = 36.0,
   aircraft_pres_int        = 2500.0,
   superob_sat_winds        = .false.,
   sat_wind_horiz_int       = 100.0,
   sat_wind_pres_int        = 2500.0,
/ 
 &obs_def_radar_mod_nml
   apply_ref_limit_to_obs      =   .false.,
   reflectivity_limit_obs      =  -10.0,
   lowest_reflectivity_obs     =  -10.0,
   apply_ref_limit_to_fwd_op   =   .false.,
   reflectivity_limit_fwd_op   =  -10.0,
   lowest_reflectivity_fwd_op  =  -10.0,
   max_radial_vel_obs          =   1000000,
   allow_wet_graupel           =   .false.,
   microphysics_type           =       2  ,
   allow_dbztowt_conv          =   .false.,
   dielectric_factor           =  0.224,
   n0_rain                     =  8.0e6,
   n0_graupel                  =  4.0e6,
   n0_snow                     =  3.0e6,
   rho_rain                    = 1000.0,
   rho_graupel                 =  400.0,
   rho_snow                    =  100.0,
/
 &assim_model_nml
   write_binary_restart_files      = .false.,
   netCDF_large_file_support       = .false.,
/
 &model_nml
   default_state_variables     = .false.,
   wrf_state_variables         = ${NL_MY_STATE_VARIABLES}
   wrf_state_bounds            = ${NL_MY_STATE_BOUNDS}
   output_state_vector         = .false.,
   num_domains                 = ${NL_NUM_DOMAINS},
   calendar_type               = 3,
   assimilation_period_seconds = 21600,
   vert_localization_coord     = 3,
   center_search_half_length   = 500000.0,
   center_spline_grid_scale    = 10   
   sfc_elev_max_diff           = ${NL_SFC_ELEV_MAX_DIFF},
   circulation_pres_level      = 80000.0,
   circulation_radius          = 108000.0,
/
 &location_nml
   horiz_dist_only                 = .false.,
   vert_normalization_pressure     = 187500.0,
   vert_normalization_height       = 5000000.0,
   vert_normalization_level        = 2666.7,
   approximate_distance            = .false.,
   nlon                            = 141,
   nlat                            = 72,
   output_box_info                 = .false.,
/
 &dart_to_wrf_nml
   model_advance_file           = ${NL_MODEL_ADVANCE_FILE},
   adv_mod_command              = '${NL_ADV_MOD_COMMAND}',
   dart_restart_name            = '${NL_DART_RESTART_NAME}',
   print_data_ranges            = .false.,
   debug                        = .false.,
/
 &wrf_to_dart_nml
   dart_restart_name   = '${NL_DART_RESTART_NAME}',
   print_data_ranges   = .false.,
   debug               = .false.,
/ 
 &utilities_nml
   TERMLEVEL                   = 1,
   module_details              = .false.,
   logfilename                 = 'dart_log.out',
   nmlfilename                 = 'dart_log.nml',
   write_nml                   = 'file',
/
 &reg_factor_nml
   select_regression           = 1,
   input_reg_file              = "time_mean_reg",
   save_reg_diagnostics        = .false.,
   reg_diagnostics_file        = 'reg_diagnostics',
/
 &ensemble_manager_nml
   single_restart_file_in  = ${NL_SINGLE_RESTART_FILE_IN},
   single_restart_file_out = ${NL_SINGLE_RESTART_FILE_OUT},
   perturbation_amplitude  = 0.2,
/
 &obs_diag_nml
   obs_sequence_name          = 'obs_seq.final',
   first_bin_center           = $YEAR_INIT, $MONTH_INIT, $DAY_INIT, $HOUR_INIT, 0, 0,
   last_bin_center            = $YEAR_END,  $MONTH_END,  $DAY_END,  $HOUR_END,  0, 0,
   bin_separation             = 0, 0, 0, 6, 0, 0 ,
   bin_width                  = 0, 0, 0, 6, 0, 0 ,
   time_to_skip               = 0, 0, 0, 0, 0, 0 ,
   max_num_bins               = 1000,
   rat_cri                    = 5000.0,
   input_qc_threshold         = 4.0,
   Nregions                   = 1,
   lonlim1                    = 235.0,
   lonlim2                    = 295.0,
   latlim1                    = 25.0,
   latlim2                    = 55.0,
   reg_names                  = 'North America',
   print_mismatched_locs      = .false.,
   print_obs_locations        = .false.,
   verbose                    = .false.,
/
 &schedule_nml
   calendar        = 'Gregorian',  
   first_bin_start =  2008, 9, 4, 3, 0, 0 ,
   first_bin_end   =  2008, 9, 4, 9, 0, 0 ,
   last_bin_end    =  2008, 9, 15, 0, 0, 0 ,
   bin_interval_days    = 0,
   bin_interval_seconds = 21600,   
   max_num_bins         = 1000,
   print_table          = .true.,
/   
 &obs_sequence_tool_nml
   num_input_files   = 1, 
   filename_seq      = 'obs_seq.out',
   filename_out      = 'obs_seq.processed', 
   first_obs_days           = 148816,
   first_obs_seconds        = 75601,
   last_obs_days            = 148817,
   last_obs_seconds         = 10800,
   obs_types         = '', 
   keep_types        = .false., 
   print_only        = .false., 
   min_lat           = -90.0, 
   max_lat           =  90.0, 
   min_lon           =   0.0, 
   max_lon           = 360.0,
/
 &restart_file_utility_nml
   input_file_name              = "filter_ic_post",
   output_file_name             = "assim_model_state_ic",
   ens_size                     = 10,
   single_restart_file_in       = ${NL_SINGLE_RESTART_FILE_IN},
   single_restart_file_out      = ${NL_SINGLE_RESTART_FILE_OUT},
   write_binary_restart_files   = .false.,
   overwrite_data_time          = .false.,
   new_data_days                = -1,
   new_data_secs                = -1,
   input_is_model_advance_file  = .false.,
   output_is_model_advance_file = .true.,
   overwrite_advance_time       = .true.,
   new_advance_days             = 148821,
   new_advance_secs             = 43200, 
/
 &restart_file_tool_nml
   input_file_name              = "filter_ic_post",
   output_file_name             = "assim_model_state_ic",
   ens_size                     = 10,
   single_restart_file_in       = ${NL_SINGLE_RESTART_FILE_IN},
   single_restart_file_out      = ${NL_SINGLE_RESTART_FILE_OUT},
   write_binary_restart_files   = .false.,
   overwrite_data_time          = .false.,
   new_data_days                = -1,
   new_data_secs                = -1,
   input_is_model_advance_file  = .false.,
   output_is_model_advance_file = .true.,
   overwrite_advance_time       = .true.,
   new_advance_days             = 148821,
   new_advance_secs             = 43200, 
/
 &replace_wrf_fields_nml
   debug = .false.,
   fail_on_missing_field = .false.,
   fieldnames = "SNOWC",
                "ALBBCK",
                "TMN",
                "TSK",
                "SH2O",
                "SMOIS",
                "SEAICE",
                "HGT_d01",
                "TSLB",
                "SST",
                "SNOWH",
                "SNOW",
   fieldlist_file = '',
/
 &obs_seq_to_netcdf_nml
   obs_sequence_name = 'obs_seq.final',
   obs_sequence_list = '',
   append_to_netcdf  = .false.,
   lonlim1    =    0.0,
   lonlim2    =  360.0,
   latlim1    =  -90.0,
   latlim2    =   90.0,
   verbose    = .false.,
/
 &obs_seq_coverage_nml
   obs_sequence_list = 'obs_coverage_list.txt',
   obs_sequence_name = '',
   obs_of_interest   = 'METAR_U_10_METER_WIND',
   textfile_out      = 'METAR_U_10_METER_WIND_obsdef_mask.txt',
   netcdf_out        = 'METAR_U_10_METER_WIND_obsdef_mask.nc',
   lonlim1    =    0.0,
   lonlim2    =  360.0,
   latlim1    =  -90.0,
   latlim2    =   90.0,
   nTmin      =      8,
   nTmax      =      8,
   verbose    = .false.,
/
 &obs_selection_nml
   filename_seq          = 'obs_seq.out',
   filename_seq_list     = '',
   filename_out          = 'obs_seq.processed',
   selections_file       = 'obsdef_mask.txt',
   selections_is_obs_seq = .false.,
   print_only            = .false.,
   calendar              = 'gregorian', 
/ 
EOF
