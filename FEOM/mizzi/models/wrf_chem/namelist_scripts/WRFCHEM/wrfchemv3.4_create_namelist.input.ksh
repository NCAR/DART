#!/bin/ksh -x
#########################################################################
#
# Purpose: Script to create WRFCHEM Namelist 
#
#########################################################################
#
# CREATE WRF NAMELIST FILE
rm -f namelist.input
touch namelist.input
cat > namelist.input << EOF
&time_control
run_days                            = ${NL_RUN_DAYS:-1},
run_hours                           = ${NL_RUN_HOURS:-0},
run_minutes                         = ${NL_RUN_MINUTES:-0},
run_seconds                         = ${NL_RUN_SECONDS:-0},
start_year                          = ${NL_START_YEAR:-2001},
start_month                         = ${NL_START_MONTH:-06},
start_day                           = ${NL_START_DAY:-11},
start_hour                          = ${NL_START_HOUR:-12},
start_minute                        = ${NL_START_MINUTE:-00},
start_second                        = ${NL_START_SECOND:-00},
end_year                            = ${NL_END_YEAR:-2001},
end_month                           = ${NL_END_MONTH:-06},
end_day                             = ${NL_END_DAY:-12},
end_hour                            = ${NL_END_HOUR:-12},
end_minute                          = ${NL_END_MINUTE:-00},
end_second                          = ${NL_END_SECOND:-00},
interval_seconds                    = ${NL_INTERVAL_SECONDS:-10800},
input_from_file                     = ${NL_INPUT_FROM_FILE:-.true.},
history_interval                    = ${NL_HISTORY_INTERVAL:-60},
frames_per_outfile                  = ${NL_FRAMES_PER_OUTFILE:-1},
restart                             = ${NL_RESTART:-.false.},
restart_interval                    = ${NL_RESTART_INTERVAL:-1440},
cycling                             = ${NL_CYCLING:-.false.},
io_form_history                     = ${NL_IO_FORM_HISTORY:-2},
io_form_restart                     = ${NL_IO_FORM_RESTART:-2},
io_form_input                       = ${NL_IO_FORM_INPUT:-2},
io_form_boundary                    = ${NL_IO_FORM_BOUNDARY:-2},
auxinput5_inname                    = ${NL_IO_AUXINPUT5_INNAME:-'null'},
auxinput6_inname                    = ${NL_IO_AUXINPUT6_INNAME:-'null'},
auxinput7_inname                    = ${NL_IO_AUXINPUT7_INNAME:-'null'},
auxinput5_interval_m                = ${NL_AUXINPUT5_INTERVAL_M:-60},
auxinput6_interval_m                = ${NL_AUXINPUT6_INTERVAL_M:-60},
auxinput7_interval_m                = ${NL_AUXINPUT7_INTERVAL_M:-60},
frames_per_auxinput5                = ${NL_FRAMES_PER_AUXINPUT5:-1},
frames_per_auxinput6                = ${NL_FRAMES_PER_AUXINPUT6:-1},
frames_per_auxinput7                = ${NL_FRAMES_PER_AUXINPUT7:-1},
io_form_auxinput5                   = ${NL_IO_FORM_AUXINPUT5:-2},
io_form_auxinput6                   = ${NL_IO_FORM_AUXINPUT6:-2},
io_form_auxinput7                   = ${NL_IO_FORM_AUXINPUT7:-2},
iofields_filename                   = ${NL_IOFIELDS_FILENAME:-'null'},
write_input                         = ${NL_WRITE_INPUT:-.true},
inputout_interval                   = ${NL_INPUTOUT_INTERVAL:-180},
input_outname                       = ${NL_INPUT_OUTNAME:-'wrf_3dvar_input_d<domain>_<date>'},
debug_level                         = ${NL_DEBUG_LEVEL:-0},
/
&domains
time_step                           = ${NL_TIME_STEP:-60},
time_step_fract_num                 = ${NL_TIME_STEP_FRACT_NUM:-0},
time_step_fract_den                 = ${NL_TIME_STEP_FRACT_DEN:-1},
max_dom                             = ${NL_MAX_DOM:-1},
e_we                                = ${NL_E_WE:-91},
e_sn                                = ${NL_E_SN:-82},
e_vert                              = ${NL_E_VERT:-28},
p_top_requested                     = ${NL_P_TOP_REQUESTED:-50000},
interp_type                         = ${NL_INTERP_TYPE:-2}, 
t_extrap_type                       = ${NL_T_EXTRAP_TYPE:-2},
num_metgrid_levels                  = ${NL_NUM_METGRID_LEVELS:-40},
num_metgrid_soil_levels             = ${NL_NUM_METGRID_SOIL_LEVELS:-4},
dx                                  = ${NL_DX:-10000},
dy                                  = ${NL_DY:-10000},
grid_id                             = ${NL_GRID_ID:-1},
parent_id                           = ${NL_PARENT_ID:-0},
i_parent_start                      = ${NL_I_PARENT_START:-1},
j_parent_start                      = ${NL_J_PARENT_START:-1},
parent_grid_ratio                   = ${NL_PARENT_GRID_RATIO:-1},
parent_time_step_ratio              = ${NL_PARENT_TIME_STEP_RATIO:-1},
feedback                            = ${NL_FEEDBACK:-1},
smooth_option                       = ${NL_SMOOTH_OPTION:-0},
eta_levels                          = ${NL_ETA_LEVELS:--1},
starting_time_step                  = ${NL_STARTING_TIME_STEP:--1},
use_adaptive_time_step              = ${NL_USE_ADAPTIVE_TIME_STEP:-false},
force_sfc_in_vinterp                = ${NL_FORCE_SFC_IN_VINTERP:-1},
max_step_increase_pct               = ${NL_MAX_STEP_INCREASE_PCT:-5},
/
&physics
mp_physics                          = ${NL_MP_PHYSICS:-2},
ra_lw_physics                       = ${NL_RA_LW_PHYSICS:-1},
ra_sw_physics                       = ${NL_RA_SW_PHYSICS:-2},
radt                                = ${NL_RADT:-40},
sf_sfclay_physics                   = ${NL_SF_SFCLAY_PHYSICS:-2},
sf_surface_physics                  = ${NL_SF_SURFACE_PHYSICS:-2},
bl_pbl_physics                      = ${NL_BL_PBL_PHYSICS:-2},
bldt                                = ${NL_BLDT:-0},
cu_physics                          = ${NL_CU_PHYSICS:-3},
cudt                                = ${NL_CUDT:-0},
isfflx                              = ${NL_ISFFLX:-1},
ifsnow                              = ${NL_IFSNOW:-0},
icloud                              = ${NL_ICLOUD:-0},
surface_input_source                = ${NL_SURFACE_INPUT_SOURCE:-1},
num_soil_layers                     = ${NL_NUM_SOIL_LAYERS:-4},
sf_urban_physics                    = ${NL_SF_URBAN_PHYSICS:-0},
maxiens                             = ${NL_MAXIENS:-1},
maxens                              = ${NL_MAXENS:-3},
maxens2                             = ${NL_MAXENS2:-3},
maxens3                             = ${NL_MAXENS3:-16},
ensdim                              = ${NL_ENSDIM:-144},
mp_zero_out                         = ${NL_MP_ZERO_OUT:-2},
cu_rad_feedback                     = ${NL_CU_RAD_FEEDBACK:-.false.},
cu_diag                             = ${NL_CU_DIAG:-1},
progn                               = ${NL_PROGN:-0},
cugd_avedx                          = ${NL_CUGD_AVEDX:-1},
num_land_cat                        = ${NL_NUM_LAND_CAT:-24},
/
&fdda
/
&dfi_control
/
&tc
/
&scm
/
&dynamics
use_baseparam_fr_nml                = ${NL_USE_BASEPARAM_FR_NML:-.false.},
w_damping                           = ${NL_W_DAMPING:-1},
diff_opt                            = ${NL_DIFF_OPT:-1},
km_opt                              = ${NL_KM_OPT:-4},
diff_6th_opt                        = ${NL_DIFF_6TH_OPT:-0},
diff_6th_factor                     = ${NL_DIFF_6TH_FACTOR:-0.12},
base_temp                           = ${NL_BASE_TEMP:-290.0},
damp_opt                            = ${NL_DAMP_OPT:-3},
zdamp                               = ${NL_ZDAMP:-5000},
dampcoef                            = ${NL_DAMPCOEF:-0.01},
khdif                               = ${NL_KHDIF:-0},
kvdif                               = ${NL_KVDIF:-0},
non_hydrostatic                     = ${NL_NON_HYDROSTATIC:-.true.},
moist_adv_opt                       = ${NL_MOIST_ADV_OPT:-2},
scalar_adv_opt                      = ${NL_SCALAR_ADV_OPT:-2},
time_step_sound                     = ${NL_TIME_STEP_SOUND:-0},
rk_ord                              = ${NL_RK_ORD:-3},
moist_adv_opt                       = ${NL_MOIST_ADV_OPT:-2},
scalar_adv_opt                      = ${NL_SCALAR_ADV_OPT:-2},
chem_adv_opt                        = ${NL_CHEM_ADV_OPT:-2},
tke_adv_opt                         = ${NL_TKE_ADV_OPT:-2},
/
&bdy_control
spec_bdy_width                      = ${NL_SPEC_BDY_WIDTH:-5},
spec_zone                           = ${NL_SPEC_ZONE:-1},
relax_zone                          = ${NL_RELAX_ZONE:-4},
specified                           = ${NL_SPECIFIED:-.false.},
nested                              = ${NL_NESTED:-.false.},
real_data_init_type                 = ${NL_REAL_DATA_INIT_TYPE:-3},
/
&grib2
/
&namelist_quilt
nio_tasks_per_group                 = ${NL_NIO_TASKS_PER_GROUP:-0},
nio_groups                          = ${NL_NIO_GROUPS:-1},
/
&chem
kemit                              = ${NL_KEMIT:-10},
chem_opt                           = ${NL_CHEM_OPT:-0},
bioemdt                            = ${NL_BIOEMDT:-0},
photdt                             = ${NL_PHOTDT:-0},
chemdt                             = ${NL_CHEMDT:-0},
io_style_emissions                 = ${NL_IO_STYLE_EMISSIONS:-2},
emiss_inpt_opt                     = ${NL_EMISS_INPT_OPT:-0},
emiss_opt                          = ${NL_EMISS_OPT:-0},
chem_in_opt                        = ${NL_CHEM_IN_OPT:-0},
phot_opt                           = ${NL_PHOT_OPT:-0},
gas_drydep_opt                     = ${NL_GAS_DRYDEP_OPT:-0},
aer_drydep_opt                     = ${NL_AER_DRYDEP_OPT:-0},
bio_emiss_opt                      = ${NL_BIO_EMISS_OPT:-0},
gas_bc_opt                         = ${NL_GAS_BC_OPT:-1},
gas_ic_opt                         = ${NL_GAS_IC_OPT:-1},
aer_bc_opt                         = ${NL_AER_BC_OPT:-1},
aer_ic_opt                         = ${NL_AER_IC_OPT:-1},
gaschem_onoff                      = ${NL_GASCHEM_ONOFF:-0},
aerchem_onoff                      = ${NL_AERCHEM_ONOFF:-0},
wetscav_onoff                      = ${NL_WETSCAV_ONOFF:-0},
cldchem_onoff                      = ${NL_CLDCHEM_ONOFF:-0},
vertmix_onoff                      = ${NL_VERTMIX_ONOFF:-0},
chem_conv_tr                       = ${NL_CHEM_CONV_TR:-0},
seas_opt                           = ${NL_SEAS_OPT:-0}, 
dust_opt                           = ${NL_DUST_OPT:-0}, 
dmsemis_opt                        = ${NL_DMSEMIS_OPT:-0}, 
biomass_burn_opt                   = ${NL_BIOMASS_BURN_OPT:-0},
plumerisefire_frq                  = ${NL_PLUMERISEFIRE_FRQ:-60},
have_bcs_chem                      = ${NL_HAVE_BCS_CHEM:-.false.},
aer_ra_feedback                    = ${NL_AER_RA_FEEDBACK:-0},
ne_area                            = ${NL_NE_AREA:-118}, 
opt_pars_out                       = ${NL_OPT_PARS_OUT:-0}, 
scale_fire_emiss                   = ${NL_SCALE_FIRE_EMISS:-.false.},
have_bcs_upper                     = ${NL_HAVE_BCS_UPPER:-.false.},
fixed_ubc_inname                   = ${NL_FIXED_UBC_INNAME:-'null'},
chemdiag                           = ${NL_CHEMDIAG:-0},
/
&wrfvar1
print_detail_grad                   = ${NL_PRINT_DETAIL_GRAD:-false},
var4d                               = ${NL_VAR4D:-false},
multi_inc                           = ${NL_MULTI_INC:-0},
/
&wrfvar2
/
&wrfvar3
ob_format                           = ${NL_OB_FORMAT:-2},
/
&wrfvar4
use_synopobs                        = ${NL_USE_SYNOPOBS:-true},
use_shipsobs                        = ${NL_USE_SHIPOBS:-true},
use_metarobs                        = ${NL_USE_METAROBS:-true},
use_soundobs                        = ${NL_USE_SOUNDOBS:-true},
use_pilotobs                        = ${NL_USE_PILOTOBS:-true},
use_airepobs                        = ${NL_USE_AIREOBS:-true},
use_geoamvobs                       = ${NL_USE_GEOAMVOBS:-true},
use_polaramvobs                     = ${NL_USE_POLARAMVOBS:-true},
use_bogusobs                        = ${NL_USE_BOGUSOBS:-true},
use_buoyobs                         = ${NL_USE_BUOYOBS:-true},
use_profilerobs                     = ${NL_USE_PROFILEROBS:-true},
use_satemobs                        = ${NL_USE_SATEMOBS:-true},
use_gpspwobs                        = ${NL_USE_GPSPWOBS:-true},
use_gpsrefobs                       = ${NL_USE_GPSREFOBS:-true},
use_ssmisobs                        = ${NL_USE_SSMISOBS:-false},
use_qscatobs                        = ${NL_USE_QSCATOBS:-true},
use_airsretobs                      = ${NL_USE_AIRSRETOBS:-true},
/
&wrfvar5
check_max_iv                        = ${NL_CHECK_MAX_IV:-true},
put_rand_seed                       = ${NL_PUT_RAND_SEED:-false},
/
&wrfvar6
ntmax                               = ${NL_NTMAX:-200},
/
&wrfvar7
cv_options                          = ${NL_CV_OPTIONS:-5},
as1                                 = ${NL_AS1_1:-0.25}, ${NL_AS1_2:-1.0}, ${NL_AS1_3:-1.5}
as2                                 = ${NL_AS2_1:-0.25}, ${NL_AS2_2:-1.0}, ${NL_AS2_3:-1.5}
as3                                 = ${NL_AS3_1:-0.25}, ${NL_AS3_2:-1.0}, ${NL_AS3_3:-1.5}
as4                                 = ${NL_AS4_1:-0.25}, ${NL_AS4_2:-1.0}, ${NL_AS4_3:-1.5}
as5                                 = ${NL_AS5_1:-0.25}, ${NL_AS5_2:-1.0}, ${NL_AS5_3:-1.5}
var_scaling4                        = ${NL_VAR_SCALING4:-1.0},
je_factor                           = ${NL_JE_FACTOR:-1.0},
/
&wrfvar8
/
&wrfvar9
/
trace_use                           = ${NL_TRACE_USE:-false},
/
&wrfvar10
/
&wrfvar11
cv_options_hum                      = ${NL_CV_OPTIONS_HUM:-1},
check_rh                            = ${NL_CHECK_RH:-0},
seed_array1                         = ${NL_SEED_ARRAY1:-2012120100},
seed_array2                         = ${NL_SEED_ARRAY2:-2012120100},
calculate_cg_cost_fn                = ${NL_CALCULATE_CG_COST_FN:-false},
lat_stats_option                    = ${NL_LAT_STATS_OPTION:-false},
/
&wrfvar12
/
&wrfvar13
/
&wrfvar14
/
&wrfvar15
num_pseudo                          = ${NL_NUM_PSEUDO:-0},
pseudo_x                            = ${NL_PSEUDO_X:-1.0},
pseudo_y                            = ${NL_PSEUDO_Y:-1.0},
pseudo_z                            = ${NL_PSEUDO_Z:-1.0},
pseudo_err                          = ${NL_PSEUDO_ERR:-1.0},
pseudo_val                          = ${NL_PSEUDO_VAL:-1.0}
/
&wrfvar16
ensdim_alpha                        = ${NL_ENSDIM_ALPHA:-0},
alphacv_method                      = ${NL_ALPHACV_METHOD:-2},
alpha_corr_type                     = ${NL_ALPHA_CORR_TYPE:-3},
alpha_corr_scale                    = ${NL_ALPHA_CORR_SCALE:-1500.0},
alpha_std_dev                       = ${NL_ALPHA_STD_DEV:-1.0},
alpha_vertloc                       = ${NL_ALPHA_VERTLOC:-false},
alpha_truncation                    = ${NL_ALPHA_TRUNCATION:-0},
/
&wrfvar17
analysis_type                       = ${NL_ANALYSIS_TYPE:-'3DVAR'},
/
&wrfvar18
analysis_date                       = ${NL_ANALYSIS_DATE:-'2012-12-01_00:00:00'},
/
&wrfvar19
pseudo_var                          = ${NL_PSEUDO_VAR:-'t'},
/
&wrfvar20
/
&wrfvar21
time_window_min                     = ${NL_TIME_WINDOW_MIN:-'2012-11-30_23:00:00'},
/
&wrfvar22
time_window_max                     = ${NL_TIME_WINDOW_MAX:-'2012-13-01_01:00:00'},
/
&wrfvar23
jcdfi_use                           = ${NL_JCDFI_USE:-false},
jcdfi_io                            = ${NL_JCDFI_IO:-false},
/
EOF
