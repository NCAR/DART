#! /bin/csh -f
#-----------------------------------------------------------------------
# Script run_wrfreal_remote.csh
#
# Purpose: Run WRFV1 real.exe.
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# [1] Set defaults for required environment variables:
#-----------------------------------------------------------------------

setenv NC                 $1
setenv START_DATE         $2
setenv FCST_RANGE         $3
setenv INTERVAL           $4
setenv WRF_DT             $5
setenv OUT_FREQ           $6
setenv WRF_DIR            $7
setenv WEST_EAST_GRIDS	  $8
setenv SOUTH_NORTH_GRIDS  $9
setenv VERTICAL_GRIDS	  $10
setenv GRID_DISTANCE	  $11

echo ${START_DATE}

source /ocotillo/users/${USER}/bin/get_date_range.csh $START_DATE $FCST_RANGE

setenv INTERVAL_SS `expr $INTERVAL \* 3600`

#-----------------------------------------------------------------------
# [2] Create WRF namelist.input:
#-----------------------------------------------------------------------

cat >! ${WRF_DIR}/test/em_real/namelist.input << EOF
 &time_control
 run_days                            = 0,
 run_hours                           = ${FCST_RANGE},
 run_minutes                         = 0,
 run_seconds                         = 0,
 start_year                          = ${START_YEAR},  ${START_YEAR},  ${START_YEAR},
 start_month                         = ${START_MONTH}, ${START_MONTH}, ${START_MONTH},
 start_day                           = ${START_DAY},   ${START_DAY},   ${START_DAY},
 start_hour                          = ${START_HOUR},  ${START_HOUR},  ${START_HOUR},
 start_minute                        = 00,   00,   00,
 start_second                        = 00,   00,   00,
 end_year                            = ${END_YEAR},  ${END_YEAR},  ${END_YEAR},
 end_month                           = ${END_MONTH}, ${END_MONTH}, ${END_MONTH},
 end_day                             = ${END_DAY},   ${END_DAY},   ${END_DAY},
 end_hour                            = ${END_HOUR},  ${END_HOUR},  ${END_HOUR},
 end_minute                          = 00,   00,   00,
 end_second                          = 00,   00,   00,
 interval_seconds                    = ${INTERVAL_SS},
 input_from_file                     = .true.,.false.,.false.,
 history_interval                    = 360,  60,   60,
 frames_per_outfile                  = 1000, 1000, 1000,
 restart                             = .false.,
 restart_interval                    = 5000,
 io_form_history                     = 2
 io_form_restart                     = 2
 io_form_input                       = 2
 io_form_boundary                    = 2
 debug_level                         = 0
 /

 &domains
 time_step                           = ${WRF_DT},
 time_step_fract_num                 = 0,
 time_step_fract_den                 = 1,
 max_dom                             = 1,
 s_we                                = 1,     1,     1,
 e_we                                = ${WEST_EAST_GRIDS},    112,   94,
 s_sn                                = 1,     1,     1,
 e_sn                                = ${SOUTH_NORTH_GRIDS},    97,    91,
 s_vert                              = 1,     1,     1,
 e_vert                              = ${VERTICAL_GRIDS},    28,    28,
 dx                                  = ${GRID_DISTANCE}, 10000,  3333,
 dy                                  = ${GRID_DISTANCE}, 10000,  3333,
 grid_id                             = 1,     2,     3,
 level                               = 1,     2,     3,
 parent_id                           = 0,     1,     2,
 i_parent_start                      = 0,     31,    30,
 j_parent_start                      = 0,     17,    30,
 parent_grid_ratio                   = 1,     3,     3,
 parent_time_step_ratio              = 1,     3,     3,
 feedback                            = 1,
 smooth_option                       = 0
 /

 &physics
 mp_physics                          = 3,     3,     3,
 ra_lw_physics                       = 1,     1,     1,
 ra_sw_physics                       = 1,     1,     1,
 radt                                = 30,    30,    30,
 sf_sfclay_physics                   = 1,     1,     1,
 sf_surface_physics                  = 1,     1,     1,
 bl_pbl_physics                      = 1,     1,     1,
 bldt                                = 0,     0,     0,
 cu_physics                          = 1,     1,     0,
 cudt                                = 5,     5,     5,
 isfflx                              = 1,
 ifsnow                              = 0,
 icloud                              = 1,
 surface_input_source                = 1,
 num_soil_layers                     = 5,
 maxiens                             = 1,
 maxens                              = 3,
 maxens2                             = 3,
 maxens3                             = 16,
 ensdim                              = 144,
 /

 &dynamics
 dyn_opt                             = 2,
 rk_ord                              = 3,
 w_damping                           = 0,
 diff_opt                            = 0,
 km_opt                              = 1,
 damp_opt                            = 0,
 zdamp                               = 5000.,  5000.,  5000.,
 dampcoef                            = 0.2,    0.2,    0.2
 khdif                               = 0,      0,      0,
 kvdif                               = 0,      0,      0,
 smdiv                               = 0.1,    0.1,    0.1,
 emdiv                               = 0.01,   0.01,   0.01,
 epssm                               = 0.1,    0.1,    0.1
 non_hydrostatic                     = .true., .true., .true.,
 time_step_sound                     = 4,      4,      4,
 h_mom_adv_order                     = 5,      5,      5,
 v_mom_adv_order                     = 3,      3,      3,
 h_sca_adv_order                     = 5,      5,      5,
 v_sca_adv_order                     = 3,      3,      3,
 /

 &bdy_control
 spec_bdy_width                      = 5,
 spec_zone                           = 1,
 relax_zone                          = 4,
 specified                           = .true., .false.,.false.,
 periodic_x                          = .false.,.false.,.false.,
 symmetric_xs                        = .false.,.false.,.false.,
 symmetric_xe                        = .false.,.false.,.false.,
 open_xs                             = .false.,.false.,.false.,
 open_xe                             = .false.,.false.,.false.,
 periodic_y                          = .false.,.false.,.false.,
 symmetric_ys                        = .false.,.false.,.false.,
 symmetric_ye                        = .false.,.false.,.false.,
 open_ys                             = .false.,.false.,.false.,
 open_ye                             = .false.,.false.,.false.,
 nested                              = .false., .true., .true.,
 /

 &namelist_quilt
 nio_tasks_per_group = 0,
 nio_groups = 1,
 /

EOF

#-----------------------------------------------------------------------
# [3.0] Run real.exe:
#-----------------------------------------------------------------------

cd ${WRF_DIR}/test/em_real

echo "   Running real.exe" 
real.exe >>& ./run_real_${NC}.out

echo ""

exit (0)
