#!/bin/ksh 
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Purpose: run real.exe and generate perturbed IC/BC's using WRFVAR
# and prepares the run directory by linking all of the needed files
# output from this routine.
#
# This script was originally written by Yongsheng Chen and modified
# by Hui Liu, Glen Romine and Nancy Collins. Please contact 
# dart@ucar.edu with any bugs or suggested revisions. 
#
# There is a section labeled 'PATHS' below which needs to be set to
# point to your local directory trees.
# 
# You should have already built and tested your DART and wrf model
# builds prior to trying this script.
#
# Options that you commonly modify can be placed in the 'OPTIONS' 
# block below, while other options are less frequently tuned in the 
# are accessible under the 'NAMELISTS' block below. The idea here 
# is to reduce the risk of inconsistent namelist options by having 
# only one place where options are defined.
#
# Things you are assumed to have for this script to work:
# 1. real.exe - compiled for serial execution
# 2. da_wrfvar.exe - compiled for serial execution
# 3. dart_to_wrf  |
# 4. wrf_to_dart  | - these should all be in $DART/models/wrf/work
# 5. pert_wrf_bc  |
# 6. da_advance_time.exe - tool for converting time formats that
#                        - should be in the WRFDA build
# 7. be.dat - background error stats data, placed or linked into
#             the ICBC directory
# 8. met_em* files - output from WPS
#
# If you have nests, you also need the file below
# 9. ndown.exe - compiled for mpi execution
#
# This script also uses the ncks utility within the NCO distribution
#
# Note that the number of processors needs to equal the number of 
# ensemble members for running this script on bluefire.
# ########
# CAUTION - this script will remove the OUTPUT and ASSIM_DIR
# directories during execution. Change to different paths to
# avoid wiping out old directories.
# ########
#-----------------------------------------------------------------------
#BSUB -n 32
#BSUB -a poe
#BSUB -R "span[ptile=16]"
#BSUB -x
#BSUB -J icbc
#BSUB -o output.icbc
#BSUB -e err.icbc
#BSUB -q regular
# EDIT for a valid project number if applicable
#BSUB -P 111111111
#BSUB -W 0:30

set echo

# SYSTEM SPECIFIC SETTINGS
# For the system you are running on, indicate the number of procs per node
export      PROCS_PER_NODE=32

# Initial Condition files
# Set to '1' if you want a single IC file, any other # if you want separate files (the 
# latter is suggested if you have large grids and lots of members)
export      single_file=1

#Time info - 4 digit year, no leading zeroes for others as this is handled later:
export      YEAR_INIT=2008
export      MONTH_INIT=5
export      DAY_INIT=22
export      HOUR_INIT=12
export      YEAR_END=2008
export      MONTH_END=5
export      DAY_END=24
export      HOUR_END=12
# the time window for the assimilation run, in hours
export      DA_TIME_WINDOW=48
# the frequency of lateral boundary condition updates, in hours
export      LBC_FREQ=6
# for the assimilation run, pre-generated BC files(0), or on-the-fly(1)
export      OTF=1
#--------------------------------------------------------------------------------------
# PATHS
#--------------------------------------------------------------------------------------
# parent tree where all of the wrf and wrfda files reside
#EDIT lower case paths - do not edit upper case path parts
export      CODE_DIR=/path/to/wrf
# where the wrf model executable dir is within this
export      WRF_DIR=$CODE_DIR/WRFV3
# where the wrfvar tree is
export      VAR_DIR=$CODE_DIR/WRFDA
# where the DART home directory is
export      DART_DIR=/path/to/dart/DART
# where the da_advance_time executable is located
export      TOOL_DIR=$VAR_DIR/var/da
# where the namelists will be written for 3dvar, DART filter (input.nml), and wrf,
# as well as the location for the be.dat file
export      ICBC_DIR=/path/to/here
# where the WPS processed files reside to feed to real.exe (met_em*)
export      WPS_DATA_DIR=/$CODE_DIR/WPS
# temporary home where the IC and BC files will be written by this script
export      OUTPUT=/path/to/output
# work directory for temporary files 
export      RUN_DIR=/path/to/work/dir
# assimilation directory where you expect to run the assimilation system....
# where you plan to run filter
export      ASSIM_DIR=/path/to/assim/run/dir
# location of the obs_seq.out file 
export      OBS_DIR=$DART_DIR/ncep_obs/work
# location of background error covariance file for da_wrfvar
export      BE_file=$VAR_DIR/var/run/be.dat.cv3

#--------------------------------------------------------------------------------------
# OPTIONS
#--------------------------------------------------------------------------------------
# 1. Common input.nml options (used by wrf_to_dart and dart_to_wrf and a copy is placed 
#                              in the assimilation directory) 
#
# how many members in your ensemble
  export ENS_SIZE=32
#Domain info:  only for nested domains, otherwise, = 1.
  export MAX_DOM=1
# List of WRF state variables to pass to DART - needs to make sense for the model
# parameter definitions - especially the moist variables related to the microphysics
# of your choice.  Note that any states that you want to update with filter, or
# monitor, needs to be on this list. Note the lists below are in double quotes.
  export  my_state_variables="'U','QTY_U_WIND_COMPONENT','TYPE_U','UPDATE','999',
                              'V','QTY_V_WIND_COMPONENT','TYPE_V','UPDATE','999',
                              'W','QTY_VERTICAL_VELOCITY','TYPE_W','UPDATE','999',
                              'PH','QTY_GEOPOTENTIAL_HEIGHT','TYPE_GZ','UPDATE','999',
                              'T','QTY_POTENTIAL_TEMPERATURE','TYPE_T','UPDATE','999',
                              'MU','QTY_PRESSURE','TYPE_MU','UPDATE','999',
                              'QVAPOR','QTY_VAPOR_MIXING_RATIO','TYPE_QV','UPDATE','999',
                              'QCLOUD','QTY_CLOUD_LIQUID_WATER','TYPE_QC','UPDATE','999',
                              'QRAIN','QTY_RAINWATER_MIXING_RATIO','TYPE_QR','UPDATE','999',
                              'QICE','QTY_CLOUD_ICE','TYPE_QI','UPDATE','999',
                              'QSNOW','QTY_SNOW_MIXING_RATIO','TYPE_QS','UPDATE','999',
                              'U10','QTY_U_WIND_COMPONENT','TYPE_U10','UPDATE','999',
                              'V10','QTY_V_WIND_COMPONENT','TYPE_V10','UPDATE','999',
                              'T2','QTY_TEMPERATURE','TYPE_T2','UPDATE','999',
                              'TH2','QTY_POTENTIAL_TEMPERATURE','TYPE_TH2','UPDATE','999',
                              'Q2','QTY_SPECIFIC_HUMIDITY','TYPE_Q2','UPDATE','999',
                              'PSFC','QTY_PRESSURE','TYPE_PS','UPDATE','999',"
# List of state variables that should be 'positive definite' - where negative values that
# emerge from filter are reset to the value indicated below (zero). Microphysical variables
# are logically added here.
  export       my_state_bounds="'QVAPOR','0.0','NULL','CLAMP',
                                'QCLOUD','0.0','NULL','CLAMP',
                                'QRAIN','0.0','NULL','CLAMP',
                                'QICE','0.0','NULL','CLAMP',
                                'QSNOW','0.0','NULL','CLAMP',"
  export my_inf_initial_from_restart=.false.
  export my_inf_sd_initial_from_restart=.false.

# 2. Common wrf model options (used by real.exe, and placed in assimilation directory)
#  num. grid pts E-W in D01 - definition used in both wrfvar and real
  export    WE_D01=45
#  num. grid pts N-S in D01 - definition used in both wrfvar and real
  export    SN_D01=45
#  Vertical levels in D01 (make sure you have matching eta entries)- definition used in both wrfvar and real
  export    VL_D01=36
#  resolution E-W in D01 in m - definition used in both wrfvar and real
  export    DXM_D01=120000.
#  resolution N-S in D01 in m - definition used in both wrfvar and real
  export    DYM_D01=120000.
# explicit microphysics - definition used in both wrfvar and real. Don't forget
# to match up the state variable list above based on the micro package you selected
  export    MP_PHYS=4    
# surface layer physics (drag) - definition used in both wrfvar and real
  export    SFC_PHYS=1    
# land-surface layer physics (thermal) - definition used in both wrfvar and real
  export    LSF_PHYS=1    
# land-surface layer physics (thermal) - definition used in both wrfvar and real
  export    SOIL_LYR=5    

# 3. Common wrfvar options
# scaling factors for WRFDA NCEP covariances. If you increase the horizontal scale, you
# should also increase the perturbation scale
  export      DA_VAR_SCALING=0.25
  export      PSCALE=$DA_VAR_SCALING
#  export      PSCALE=0.25
  export      HSCALE=1.0
  export      VSCALE=1.5

#--------------------------------------------------------------------------------------
# NAMELISTS
#--------------------------------------------------------------------------------------
cat > input.nml.tmp << EOF
 &filter_nml
   async                    =  2,
   adv_ens_command          = "./advance_model.csh",
   ens_size                 =  $ENS_SIZE,
   start_from_restart       = .true.,
   output_restart           = .false.,
   obs_sequence_in_name     = "obs_seq.out",
   obs_sequence_out_name    = "obs_seq.final",
   restart_in_file_name     = "filter_ics",
   restart_out_file_name    = "filter_restart",
   init_time_days           = -1,
   init_time_seconds        = -1,
   first_obs_days           = -1,
   first_obs_seconds        = -1,
   last_obs_days            = -1,
   last_obs_seconds         = -1,
   num_output_state_members = 0,
   num_output_obs_members   = 0,
   output_interval          = 1,
   num_groups               = 1,
   input_qc_threshold       = 4.0,
   outlier_threshold        = 3.0,
   output_forward_op_errors = .false.,
   output_timestamps        = .false.,
   output_inflation         = .true.,
   trace_execution          = .false.,
   inf_flavor                  = 2,                      0,
   inf_initial_from_restart    = $my_inf_initial_from_restart,    .false.,
   inf_sd_initial_from_restart = $my_inf_sd_initial_from_restart, .false.,
   inf_output_restart          = .true.,                 .true.,
   inf_deterministic           = .true.,                 .true.,
   inf_in_file_name            = 'prior_inf_ic_old',     'post_inf_ic_old',
   inf_out_file_name           = 'prior_inf_ic_new',     'post_inf_ic_new',
   inf_diag_file_name          = 'prior_inf_diag',       'post_inf_diag',
   inf_initial                 = 1.00,                   1.00,
   inf_sd_initial              = 0.60,                   0.50,
   inf_damping                 = 0.90,                   1.00,
   inf_lower_bound             = 1.0,                    1.0,
   inf_upper_bound             = 1000000.0,              1000000.0,
   inf_sd_lower_bound          = 0.60,                   0.10
/
 &ensemble_manager_nml
   single_restart_file_in  = .true.,
   single_restart_file_out = .true.,
   perturbation_amplitude  = 0.2 /

 &smoother_nml
   num_lags              = 0
   start_from_restart    = .false.
   output_restart        = .false.
   restart_in_file_name  = 'smoother_ics'
   restart_out_file_name = 'smoother_restart' /

 &assim_tools_nml
   filter_kind                     = 1,
   cutoff                          = 0.16,
   sort_obs_inc                    = .false.,
   spread_restoration              = .false.,
   sampling_error_correction       = .false.,
   print_every_nth_obs             = 1000,
   adaptive_localization_threshold = -1 /

 &cov_cutoff_nml
   select_localization             = 1  /

 &assim_model_nml
   write_binary_restart_files      = .true.  /

 &location_nml
   horiz_dist_only                 = .false.,
   vert_normalization_pressure     = 187500.0,
   vert_normalization_height       = 5000000.0,
   vert_normalization_level        = 2666.7,
   approximate_distance            = .false.,
   nlon                            = 141,
   nlat                            = 72,
   output_box_info                 = .false.  /

 &model_nml
   output_state_vector         = .false.,
   default_state_variables     = .false.,
   wrf_state_variables         = $my_state_variables
   wrf_state_bounds            = $my_state_bounds
   num_domains                 = $MAX_DOM,
   surf_obs                    = .false.,
   calendar_type               = 3,
   assimilation_period_seconds = 21600,
   adv_mod_command             = "./wrf.exe",
   vert_localization_coord     = 2,
   center_search_half_length   = 400000.0,
   center_spline_grid_scale    = 10  /

 &utilities_nml
   TERMLEVEL                   = 1,
   logfilename                 = 'dart_log.out',
   nmlfilename                 = 'dart_log.nml',
   write_nml                   = 'file',
   module_details              = .false.  /

 &reg_factor_nml
   select_regression           = 1,
   input_reg_file              = "time_mean_reg",
   save_reg_diagnostics        = .false.,
   reg_diagnostics_file        = 'reg_diagnostics'  /

 &obs_sequence_nml
   write_binary_obs_sequence   = .false.  /

 &preprocess_nml
   input_obs_kind_mod_file  = '../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90',
   output_obs_kind_mod_file = '../../../assimilation_code/modules/observations/obs_kind_mod.f90',
   input_obs_def_mod_file   = '../../../observations/forward_operators/DEFAULT_obs_def_mod.F90',
   output_obs_def_mod_file  = '../../../observations/forward_operators/obs_def_mod.f90',
   input_files              = '../../../observations/forward_operators/obs_def_reanalysis_bufr_mod.f90',
                              '../../../observations/forward_operators/obs_def_altimeter_mod.f90',
                              '../../../observations/forward_operators/obs_def_radar_mod.f90',
                              '../../../observations/forward_operators/obs_def_metar_mod.f90',
                              '../../../observations/forward_operators/obs_def_dew_point_mod.f90',
                              '../../../observations/forward_operators/obs_def_gps_mod.f90',
                              '../../../observations/forward_operators/obs_def_gts_mod.f90',
                              '../../../observations/forward_operators/obs_def_QuikSCAT_mod.f90',
                              '../../../observations/forward_operators/obs_def_vortex_mod.f90'  /

 &obs_kind_nml
   assimilate_these_obs_types = 'RADIOSONDE_TEMPERATURE',
                                'METAR_U_10_METER_WIND',
                                'METAR_V_10_METER_WIND',
                                'METAR_TEMPERATURE_2_METER',
                                'RADIOSONDE_U_WIND_COMPONENT',
                                'RADIOSONDE_V_WIND_COMPONENT',
                                'SAT_U_WIND_COMPONENT',
                                'SAT_V_WIND_COMPONENT',
   evaluate_these_obs_types   = 'RADIOSONDE_SPECIFIC_HUMIDITY'  /

 &obs_diag_nml
   obs_sequence_name          = 'obs_seq.final',
   first_bin_center           =  $YEAR_INIT, $MONTH_INIT, $DAY_INIT, $HOUR_INIT, 0, 0 ,
   last_bin_center            =  $YEAR_END , $MONTH_END , $DAY_END , $HOUR_END , 0, 0 ,
   bin_separation             =     0, 0, 0, 6, 0, 0 ,
   bin_width                  =     0, 0, 0, 6, 0, 0 ,
   time_to_skip               =     0, 0, 0, 0, 0, 0 ,
   max_num_bins               = 1000,
   trusted_obs                = 'null',
   Nregions                   = 1,
   lonlim1                    =  235.0,
   lonlim2                    =  295.0,
   latlim1                    =   25.0,
   latlim2                    =   55.0,
   reg_names                  = 'North America',
   print_mismatched_locs      = .false.,
   create_rank_histogram      = .true.,
   outliers_in_histogram      = .true.,
   use_zero_error_obs         = .false.,
   verbose                    = .false.  /

 &restart_file_utility_nml
   input_file_name              = "restart_file_input",
   output_file_name             = "restart_file_output",
   ens_size                     = 1,
   single_restart_file_in       = .true.,
   single_restart_file_out      = .true.,
   write_binary_restart_files   = .true.,
   overwrite_data_time          = .false.,
   new_data_days                = -1,
   new_data_secs                = -1,
   input_is_model_advance_file  = .false.,
   output_is_model_advance_file = .true.,
   overwrite_advance_time       = .true.,
   new_advance_days             = -1,
   new_advance_secs             = -1 /

 &dart_to_wrf_nml
   restart_file = .FALSE. /
EOF

#--------------------------------------------------------------------------------------
# Need a template version for advance_model.csh as well
#
cat > namelist.wrf << EOF
 &time_control
 run_days                            = 0,
 run_hours                           = 0,
 run_minutes                         = 0,
 run_seconds                         = 0,
 start_year                          = _MAX_DOM_*_START_YEAR_,
 start_month                         = _MAX_DOM_*_START_MONTH_,
 start_day                           = _MAX_DOM_*_START_DAY_,
 start_hour                          = _MAX_DOM_*_START_HOUR_,
 start_minute                        = 00, 00, 00,
 start_second                        = 00, 00, 00,
 end_year                            = _MAX_DOM_*_END_YEAR_,
 end_month                           = _MAX_DOM_*_END_MONTH_,
 end_day                             = _MAX_DOM_*_END_DAY_,
 end_hour                            = _MAX_DOM_*_END_HOUR_,
 end_minute                          = 00,   00,   00,
 end_second                          = 00,   00,   00,
 interval_seconds                    = 21600
 input_from_file                     = .true., .true.,
 history_interval                    = 360,  360,
 frames_per_outfile                  = 1,  1,
 restart                             = .false.,
 restart_interval                    = 5000,
 io_form_history                     = 2
 io_form_restart                     = 2
 io_form_input                       = 2
 io_form_boundary                    = 2
 all_ic_times                        = .true.,
 debug_level                         = 0
 /

&domains
 time_step                           = 720,
 time_step_fract_num                 = 0,
 time_step_fract_den                 = 1,
 max_dom                             = _MAX_DOM_,
 e_we                                = $WE_D01,  90,
 e_sn                                = $SN_D01,  90,
 e_vert                              = $VL_D01,   36,
 p_top_requested                     = 5000,
 num_metgrid_levels                  = 27,
 num_metgrid_soil_levels             = 4,
 dx                                  = $DXM_D01,  40000,
 dy                                  = $DYM_D01,  40000,
 grid_id                             = 1,     2,
 parent_id                           = 0,     1,
 i_parent_start                      = 1,     25,
 j_parent_start                      = 1,     10,
 parent_grid_ratio                   = 1,     3,
 parent_time_step_ratio              = 1,     3,
 feedback                            = 1,
 smooth_option                       = 0,
 lagrange_order                      = 2
 interp_type                         = 2
 lowest_lev_from_sfc                 = .false.,
 force_sfc_in_vinterp                = 0,
 zap_close_levels                    = 500,
 sfcp_to_sfcp                        = .true.,
 use_levels_below_ground             = .true.,
 eta_levels                          = 1.00000, 0.99293, 0.98531, 0.97658, 0.96619, 0.95361, 0.93667, 0.91446, 0.88731, 0.85558, 0.81974, 0.78026, 0.73768, 0.69258, 0.64553, 0.59712, 0.54431, 0.48851, 0.43114, 0.37362, 0.31726, 0.26747, 0.22363, 0.18516, 0.15153, 0.12223, 0.09683, 0.07490, 0.05605, 0.04100, 0.03100, 0.02200, 0.01300, 0.00700, 0.00300, 0.00000,
 /

 &physics
 mp_physics                          = $MP_PHYS,     $MP_PHYS,
 ra_lw_physics                       = 1,     1,
 ra_sw_physics                       = 1,     1,
 radt                                = 30,    10,
 sf_sfclay_physics                   = $SFC_PHYS,    $SFC_PHYS,
 sf_surface_physics                  = $LSF_PHYS,    $LSF_PHYS,
 bl_pbl_physics                      = 1,     1,
 bldt                                = 0,     0,
 cu_physics                          = 1,     1,
 cudt                                = 0,     0,
 isfflx                              = 1,
 ifsnow                              = 0,
 icloud                              = 1,
 surface_input_source                = 1,
 num_soil_layers                     = $SOIL_LYR,
 sf_urban_physics                    = 0,
 mp_zero_out                         = 2,
 mp_zero_out_thresh                  = 1.e-8,
 maxiens                             = 1,
 maxens                              = 3,
 maxens2                             = 3,
 maxens3                             = 16,
 ensdim                              = 144,
 /

 &fdda
 /
 &dynamics
 w_damping                           = 1,
 diff_opt                            = 1,
 km_opt                              = 4,
 diff_6th_opt                        = 0,      0,
 diff_6th_factor                     = 0.12,   0.12,
 damp_opt                            = 0,
 zdamp                               = 5000.,  5000.,
 dampcoef                            = 0.2,    0.2,
 khdif                               = 0,      0,
 kvdif                               = 0,      0,
 non_hydrostatic                     = .true., .true.,
 moist_adv_opt                       = 1,      1,
 scalar_adv_opt                      = 1,      1,
 iso_temp                            = 200.,
 /

 &bdy_control
 spec_bdy_width                      = 5,
 spec_zone                           = 1,
 relax_zone                          = 4,
 specified                           = .true., .false.,
 nested                              = .false., .true.,
 /

 &grib2
 /

 &namelist_quilt
 nio_tasks_per_group = 0,
 nio_groups = 1,
 /
EOF

#--------------------------------------------------------------------------------------
# Note the file below is a template - which is updated for individual calls to 
# WRFVAR
cat > namelist.3dvar << EOF
 &wrfvar1
 check_max_iv_print                  = .false.,
 write_increments                    = .false.,
 /

 &wrfvar2
 /

 &wrfvar3
 /

 &wrfvar4
 use_synopobs                        = .false.,
 use_shipsobs                        = .false.,
 use_metarobs                        = .false.,
 use_soundobs                        = .false.,
 use_pilotobs                        = .false.,
 use_airepobs                        = .false.,
 use_geoamvobs                       = .false.,
 use_polaramvobs                     = .false.,
 use_bogusobs                        = .false.,
 use_buoyobs                         = .false.,
 use_profilerobs                     = .false.,
 use_satemobs                        = .false.,
 use_gpspwobs                        = .false.,
 use_gpsrefobs                       = .false.,
 use_qscatobs                        = .false.,
 use_radarobs                        = .false.,
 use_radar_rv                        = .false.,
 use_radar_rf                        = .false.,
 use_airsretobs                      = .false.,
 /

 &wrfvar5
 put_rand_seed                       = .true.,
 check_max_iv                        = .false.,
 /
 &wrfvar6
 max_ext_its                         = 1
 ntmax                               = 200,
 eps                                 = 0.01,
 /

 &wrfvar7
 cv_options                          = 3,
 as1                                 = $PSCALE, $HSCALE, $VSCALE,
 as2                                 = $PSCALE, $HSCALE, $VSCALE,
 as3                                 = $PSCALE, $HSCALE, $VSCALE,
 as4                                 = $PSCALE, $HSCALE, $VSCALE,
 as5                                 = $PSCALE, $HSCALE, $VSCALE,
 rf_passes                           = 6,
 var_scaling1                        = $DA_VAR_SCALING,
 var_scaling2                        = $DA_VAR_SCALING,
 var_scaling3                        = $DA_VAR_SCALING,
 var_scaling4                        = $DA_VAR_SCALING,
 var_scaling5                        = $DA_VAR_SCALING,
 je_factor                           = 1.0,
 /

 &wrfvar8
 /

 &wrfvar9
 trace_csv                           = .false.,
 use_html                            = .false.,
 /

 &wrfvar10
 /
 &wrfvar11
 cv_options_hum                      = 1,
 check_rh                            = 1,
 set_omb_rand_fac                    = 1.0,
 seed_array1                         = _SEED_ARRAY1_,
 seed_array2                         = _SEED_ARRAY2_,
 /

 &wrfvar12
 balance_type   = 1,
 /

 &wrfvar13
 vert_corr                           = 2,
 vertical_ip                         = 0,
 vert_evalue                         = 1,
 max_vert_var1                       = 99.0,
 max_vert_var2                       = 99.0,
 max_vert_var3                       = 99.0,
 max_vert_var4                       = 99.0,
 max_vert_var5                       = 0.0,
 /

 &wrfvar14
 /

 &wrfvar15
 num_pseudo                          = 0,
 /

 &wrfvar16
 /

 &wrfvar17
 analysis_type                       = 'RANDOMCV',
 /
 &wrfvar18
 analysis_date                       = '_ANALYSIS_DATE_.0000',
 /

 &wrfvar19
 /

 &wrfvar20
 /

 &wrfvar21
 /

 &wrfvar22
 /

 &wrfvar23
 /
 &time_control
 start_year=_START_YEAR_,
 start_month=_START_MONTH_,
 start_day=_START_DAY_,
 start_hour=_START_HOUR_,
 end_year=_END_YEAR_,
 end_month=_END_MONTH_,
 end_day=_END_DAY_,
 end_hour=_END_HOUR_,
 debug_level                         = 0
 /
 
 &domains
 e_we=$WE_D01,
 e_sn=$SN_D01,
 e_vert=$VL_D01,
 dx=$DXM_D01,
 dy=$DYM_D01,
 /
 
 &physics
 mp_physics=$MP_PHYS,
 sf_sfclay_physics=$SFC_PHYS,
 sf_surface_physics=$LSF_PHYS,
 num_soil_layers=$SOIL_LYR,/
EOF

#--------------------------------------------------------------------------------------
# You hopefully will not need to modify options below this line
#--------------------------------------------------------------------------------------

ln -sf $BE_file be.dat
#Date string needs to be a fixed format
typeset -Z2 HOUR_INIT2=$HOUR_INIT
typeset -Z2 DAY_INIT2=$DAY_INIT
typeset -Z2 MONTH_INIT2=$MONTH_INIT
export      INITIAL_DATE=$YEAR_INIT$MONTH_INIT2$DAY_INIT2$HOUR_INIT2

if [[ -d $OUTPUT ]]; then \rm -rf $OUTPUT; fi
mkdir -p $OUTPUT

\rm output.*

if [[ -d $RUN_DIR ]]; then \rm -rf $RUN_DIR; fi
mkdir -p $RUN_DIR
cd $RUN_DIR

# ----------------
# copy/link necessary files
# ----------------
ln -sf $WRF_DIR/run/* .
ln -sf $WRF_DIR/main/real.exe real.exe
ln -sf $WRF_DIR/main/ndown.exe ndown.exe
\rm namelist.input output*
ln -sf $ICBC_DIR/namelist.input.wrf namelist.input

# ----------------
# run real.exe 
# ----------------
   echo "starting real section..."
   start_date=$INITIAL_DATE
   end_date=`$TOOL_DIR/da_advance_time.exe $INITIAL_DATE $DA_TIME_WINDOW`
   fcst_hour=$DA_TIME_WINDOW
yyyy=`echo $start_date | cut -c 1-4`
  mm=`echo $start_date | cut -c 5-6`
  dd=`echo $start_date | cut -c 7-8`
  hh=`echo $start_date | cut -c 9-10`
yyyy_end=`echo $end_date | cut -c 1-4`
  mm_end=`echo $end_date | cut -c 5-6`
  dd_end=`echo $end_date | cut -c 7-8`
  hh_end=`echo $end_date | cut -c 9-10`

m4 -D_FCST_=$fcst_hour -D_MAX_DOM_=$MAX_DOM \
   -D_START_YEAR_=$yyyy -D_START_MONTH_=$mm -D_START_DAY_=$dd -D_START_HOUR_=$hh \
   -D_END_YEAR_=$yyyy_end -D_END_MONTH_=$mm_end -D_END_DAY_=$dd_end -D_END_HOUR_=$hh_end \
   $ICBC_DIR/namelist.wrf  > namelist.input

   ln -sf $WPS_DATA_DIR/met_em.d0* .

./real.exe

# ----------------
# rename files 
# ----------------

  export NANALYSIS=`expr $fcst_hour \/ $LBC_FREQ \+ 1`
##   mv wrflowinp_d01 $OUTPUT/.

  it=1
  this_date=$start_date
  echo !!!! $this_date and $g_date[0] $g_date[1]  and ${wdate}
  while [[ $it -le $NANALYSIS ]] ; do
    set -A g_date `$TOOL_DIR/da_advance_time.exe $this_date 0 -g`
    set -A w_date `$TOOL_DIR/da_advance_time.exe $this_date 0 -w`


    if [[ $it -eq 1 ]] ; then
     id=1
     while [[ $id -le $MAX_DOM ]] ; do
        echo move named mean file
        mv wrfinput_d0$id $OUTPUT/wrfinput_d0${id}_${g_date[0]}_${g_date[1]}_mean
        (( id = id + 1 ))
     done
    else
     (( itp1 = it - 1 ))
     echo move dated mean file
     mv wrfinput_d01.${w_date}  $OUTPUT/wrfinput_d01_${g_date[0]}_${g_date[1]}_mean
# Extract the individual bdy conditions for each time. The file time is the target time,
# and so should be coupled with a wrfinput file one analysis time earlier
     ncks -F -O -d Time,${itp1} wrfbdy_d01 $OUTPUT/wrfbdy_d01_${g_date[0]}_${g_date[1]}_mean
    fi
    (( it = it + 1 ))
    this_date=`$TOOL_DIR/da_advance_time.exe $this_date $LBC_FREQ`
  done
# File with all of the bdy times
  mv wrfbdy_d01 $OUTPUT/wrfbdy_d01

# ----------------------------------------------------------
# prepare data/scripts to run wrfvar (v3) to perturb icbc
# ----------------------------------------------------------
  mpmd_cmdfile=$RUN_DIR/run_wrfvar_mpmd_cmdfile
  ie=1
  while [[ $ie -le $ENS_SIZE ]]; do
    if [[ ! -d $ie ]]; then mkdir $ie ; fi 
    cd $ie
    ln -sf $VAR_DIR/run/LANDUSE.TBL .
    ln -sf $VAR_DIR/run/gribmap.txt .
#     compile wrfvar v3 in single thread
    ln -sf $VAR_DIR/var/build/da_wrfvar.exe    da_wrfvar.exe
    ln -sf $ICBC_DIR/be.dat .
    ln -sf $DART_DIR/models/wrf/work/pert_wrf_bc .

    cp $ICBC_DIR/input.nml.tmp  input.nml 

    cat > run_wrfvar_script.ksh << EOF
#!/bin/ksh -v
   # ksh run_wrfvar_script.ksh element start_date
   hostname
   ie=\$1
   this_date=\$2
   cd \$RUN_DIR/\$ie
   it=1
   while [[ \$it -le \$NANALYSIS ]] ; do
      set -A g_date \`\$TOOL_DIR/da_advance_time.exe \$this_date 0 -g\`
      this_date_wrf=\`\$TOOL_DIR/da_advance_time.exe \$this_date 0 -w\`

      seed_array1=\$this_date
      (( seed_array2 = ie * 10000 ))

      yyyy=\`echo \$this_date | cut -c 1-4\`
        mm=\`echo \$this_date | cut -c 5-6\`
        dd=\`echo \$this_date | cut -c 7-8\`
        hh=\`echo \$this_date | cut -c 9-10\`

      m4 -D_ANALYSIS_DATE_=\$this_date_wrf \
         -D_SEED_ARRAY1_=\$seed_array1 -D_SEED_ARRAY2_=\$seed_array2 \
         -D_START_YEAR_=\$yyyy -D_START_MONTH_=\$mm -D_START_DAY_=\$dd -D_START_HOUR_=\$hh \
         -D_END_YEAR_=\$yyyy -D_END_MONTH_=\$mm -D_END_DAY_=\$dd -D_END_HOUR_=\$hh \
         -D_VAR_SCALING1_=\$var_scaling1 -D_VAR_SCALING2_=\$var_scaling2 \
         -D_VAR_SCALING3_=\$var_scaling3 -D_VAR_SCALING4_=\$var_scaling4 \
         -D_VAR_SCALING5_=\$var_scaling5 \
         $ICBC_DIR/namelist.3dvar > namelist.input

      ln -sf $OUTPUT/wrfinput_d01_\${g_date[0]}_\${g_date[1]}_mean fg
# If using an mpi version of da_wrfvar, swap the commented lines below
# and swap the mpi call on the script below - see MPIDA
#      mpirun.lsf da_wrfvar.exe >> output.wrfvar 2>&1
      ./da_wrfvar.exe >> output.wrfvar 2>&1

      mv wrfvar_output \$OUTPUT/wrfinput_d01_\${g_date[0]}_\${g_date[1]}_\$ie

      # -------------------------
      # update wrfbdy
      # -------------------------

      if [[ it -eq 1 ]] ; then
         ln -sf \$OUTPUT/wrfinput_d01_\${g_date[0]}_\${g_date[1]}_\$ie wrfinput_this
      else
         ln -sf \$OUTPUT/wrfinput_d01_\${g_date[0]}_\${g_date[1]}_\$ie wrfinput_next
         cp \$OUTPUT/wrfbdy_d01_\${g_date[0]}_\${g_date[1]}_mean wrfbdy_this

#  -------   perturb the wrfbdy files  --------
         pert_wrf_bc > output.pert_wrf_bc.\${g_date[0]}_\${g_date[1]} 2>&1

         mv wrfbdy_this \$OUTPUT/wrfbdy_d01_\${g_date[0]}_\${g_date[1]}_\$ie
         ln -sf \$OUTPUT/wrfinput_d01_\${g_date[0]}_\${g_date[1]}_\$ie wrfinput_this
      fi

      (( it = it + 1 ))
      this_date=\`\$TOOL_DIR/da_advance_time.exe \$this_date \$LBC_FREQ\`
   done
EOF

   chmod +x run_wrfvar_script.ksh

   cd ..

   # ------------------------------
   # create a command-file for mpmd
   # ------------------------------

   if [[ -f $mpmd_cmdfile && $ie -eq 1 ]]; then \rm $mpmd_cmdfile ; fi
   echo $RUN_DIR/$ie/run_wrfvar_script.ksh $ie $start_date \& >> $mpmd_cmdfile
   (( ie = ie + 1 ))
done

# --------------------------------
# run mpmd job to perturb ic/bc
# --------------------------------
#   using mpi wrfvar.exe
# MPIDA - uncomment the 3 lines below for MPI wrf_var
#chmod +x $mpmd_cmdfile
#$mpmd_cmdfile
#rm $mpmd_cmdfile
#   using single-thread wrfvar.exe
# MPIDA - the command below works on bluefire with Serial wrfvar, but
# needs to have the number of processors assigned to the job equal to
# the number of ensemble members. So - the loop below is recommended
# when not running on an IBM system.
#export MP_PGMMODEL=mpmd
#mpirun.lsf -cmdfile $mpmd_cmdfile
  num_loops=`expr $ENS_SIZE / $PROCS_PER_NODE`
  num_extra=`expr $ENS_SIZE % $PROCS_PER_NODE`
  if [[ $num_extra -ne 0 ]]; then ((num_loops=num_loops+1)) ; fi
  echo $num_loops
# step through the command file in chunks the size of PROCS_PER_NODE:
  it=1
  while [[ $it -le $num_loops ]] ; do
   let step_size=$PROCS_PER_NODE*$it
   if [[ $it -le $num_loops && $num_extra -eq 0 ]]; then
     short_cmd=$RUN_DIR/short_cmdfile
     job=`head -\$step_size $mpmd_cmdfile | tail -\$PROCS_PER_NODE`
     echo $job >> $short_cmd
     echo wait >> $short_cmd
     chmod +x $short_cmd
     $short_cmd
     rm $short_cmd
   else
     short_cmd=$RUN_DIR/short_cmdfile
     job=`tail -\$num_extra $mpmd_cmdfile`
     echo $job >> $short_cmd
     echo wait >> $short_cmd
     chmod +x $short_cmd
     $short_cmd
     rm $short_cmd
   fi
   ((it=it+1))
  done


# --------------------------------

#  for use as template of WRF/DART filter only
   cp  $OUTPUT/wrfinput_d01_${g_date[0]}_${g_date[1]}_mean $OUTPUT/wrfinput_d01

# -------------------------
# create filter_ics if desired
# -------------------------

   start_date_wrf=`$TOOL_DIR/da_advance_time.exe $start_date 0 -w`
   set -A g_date `$TOOL_DIR/da_advance_time.exe $start_date 0 -g`
   yyyy=`echo $start_date | cut -c 1-4`
     mm=`echo $start_date | cut -c 5-6`
     dd=`echo $start_date | cut -c 7-8`
     hh=`echo $start_date | cut -c 9-10`

   m4 -D_FCST_=0 -D_MAX_DOM_=$MAX_DOM \
      -D_START_YEAR_=$yyyy -D_START_MONTH_=$mm -D_START_DAY_=$dd -D_START_HOUR_=$hh \
      -D_END_YEAR_=$yyyy -D_END_MONTH_=$mm -D_END_DAY_=$dd -D_END_HOUR_=$hh \
      $ICBC_DIR/namelist.wrf  > namelist.input

    cp $ICBC_DIR/input.nml.tmp  input.nml 

   ln -sf $DART_DIR/models/wrf/work/dart_to_wrf .
   ln -sf $DART_DIR/models/wrf/work/wrf_to_dart .

   if [[ MAX_DOM -gt 2 ]]; then
      echo MAX_DOM=$MAX_DOM is not supported. Stop!
      exit
   fi

   \rm filter_ic*

   ie=1
   while [[ $ie -le $ENS_SIZE ]]; do
      ln -sf $OUTPUT/wrfinput_d01_${g_date[0]}_${g_date[1]}_$ie wrfinput_d01
      if [[ MAX_DOM -eq 2 ]]; then
         ln -sf $OUTPUT/wrfinput_d01_${g_date[0]}_${g_date[1]}_$ie wrfout_d01_${start_date_wrf}
         ln -sf $OUTPUT/wrfinput_d02_${g_date[0]}_${g_date[1]}_mean wrfinput_d02
         mpirun.lsf ndown.exe
      fi

      ./wrf_to_dart 
      if [[ $single_file -eq 1 ]]; then
         cat dart_wrf_vector >> filter_ics
      else 
         typeset -Z4 tail_num=$ie
         mv dart_wrf_vector filter_ic.$tail_num
      fi
      (( ie = ie + 1 ))
   done
 \mv  filter_ic*  $OUTPUT/. 

#############################################################################
# Prepare for running filter by creating an assimilation directory and the
# needed files and links
#############################################################################

if [[ -d $ASSIM_DIR ]]; then \rm -rf $ASSIM_DIR; fi
   mkdir -p $ASSIM_DIR
   cd $ASSIM_DIR
# Copy in the input.nml file for DART
   cp $ICBC_DIR/input.nml.tmp  $ASSIM_DIR/input.nml 
# WRFDA and WRF namelists
   head -123 $ICBC_DIR/namelist.3dvar > $ICBC_DIR/namelist.tmp
   cat $ICBC_DIR/namelist.tmp $ICBC_DIR/namelist.wrf > $ASSIM_DIR/namelist.input

# File with the value of variance for 3dvar
cat > bc_pert_scale << EOF
   $PSCALE
   $HSCALE
   $VSCALE
EOF
# Link all of the executables to the run directory
# DART first
   ln -sf $DART_DIR/models/wrf/work/wrf_to_dart    .
   ln -sf $DART_DIR/models/wrf/work/dart_to_wrf    .
   ln -sf $DART_DIR/models/wrf/work/pert_wrf_bc    .
   ln -sf $DART_DIR/models/wrf/work/wakeup_filter  .
   ln -sf $DART_DIR/models/wrf/work/obs_diag       .
   ln -sf $DART_DIR/models/wrf/work/filter         .
   cp $DART_DIR/models/wrf/shell_scripts/advance_model.csh  .
   cp $DART_DIR/models/wrf/work/runme_filter       .
   cp $DART_DIR/models/wrf/work/advance_time       .
# WRF and WRFDA executables 
   mkdir -p $ASSIM_DIR/WRF_RUN
   ln -sf $WRF_DIR/run/*                           ./WRF_RUN/
   rm ./WRF_RUN/namelist.input
   ln -sf $TOOL_DIR/da_wrfvar.exe                  ./WRF_RUN/
   ln -sf $TOOL_DIR/da_advance_time.exe            ./WRF_RUN/
# Background covariance
   ln -sf $ICBC_DIR/be.dat                         ./WRF_RUN/
# Initial condition files for filter
   ln -sf $OUTPUT/filter_ic*                       .
   ln -sf $OUTPUT/wrfinput_d01                     .
# Link the observation sequence file
   ln -sf $OBS_DIR/obs_seq.out                     .
# Create directory where wrfout files are written
   mkdir -p $ASSIM_DIR/WRFOUT
# Making WRF directory and linking the IC and BC files to there
   mkdir -p $ASSIM_DIR/WRF
   cd $ASSIM_DIR/WRF
   ln -sf $OUTPUT/wrfbdy_d01                       .
   ln -sf $OUTPUT/wrfinput_d01                     .
   if [[ OTF -eq 1 ]]; then
     ln -sf $OUTPUT/wrfbdy_*mean                   .
     ln -sf $OUTPUT/wrfinput_*mean                 .
   else
     ln -sf $OUTPUT/wrfbdy_*                       .
   fi
   ie=1
   while [[ $ie -le $ENS_SIZE ]]; do
      ln -sf $OUTPUT/wrfinput_d01_${g_date[0]}_${g_date[1]}_$ie wrfinput_d01_$ie
      if [[ MAX_DOM -eq 2 ]]; then
      ln -sf $OUTPUT/wrfinput_d02_${g_date[0]}_${g_date[1]}_$ie wrfinput_d02_$ie
      fi
      (( ie = ie + 1 ))
   done

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

