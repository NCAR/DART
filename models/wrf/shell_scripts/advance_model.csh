#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

# This script copies the necessary files into the temporary directory
# for a model run. It assumes that there is ${PBS_O_WORKDIR}/WRF directory
# where boundary conditions files reside.

set PBS_O_WORKDIR = $1
set element = $2
set temp_dir = $3

# Shell script to run the WRF model from DART input.

set verbose

rm -rf $temp_dir
mkdir  $temp_dir
cd     $temp_dir

# Copy the initial condition file to the temp directory

cp ${PBS_O_WORKDIR}/wrfinput .
mv ${PBS_O_WORKDIR}/assim_model_state_ic$element dart_wrf_vector
ln -s ${PBS_O_WORKDIR}/input.nml .

ln -s  ${PBS_O_WORKDIR}/RRTM_DATA .
ln -s  ${PBS_O_WORKDIR}/LANDUSE.TBL .

# Convert DART to wrfinput

echo ".true." | ${PBS_O_WORKDIR}/dart_tf_wrf >& out.dart_to_wrf

set time = `head -1 wrf.info`

set secs = $time[1]
set days = $time[2]

# Copy the boundary condition file to the temp directory.
cp ${PBS_O_WORKDIR}/WRF/wrfbdy_${days}_${secs}_$element wrfbdy_d01

set date              = `head -2 wrf.info | tail -1`
set START_YEAR  = $date[1]
set START_MONTH = $date[2]
set START_DAY   = $date[3]
set START_HOUR  = $date[4]
set START_MIN   = $date[5]
set START_SEC   = $date[6]
set date              = `head -3 wrf.info | tail -1`
set END_YEAR    = $date[1]
set END_MONTH   = $date[2]
set END_DAY     = $date[3]
set END_HOUR    = $date[4]
set END_MIN     = $date[5]
set END_SEC     = $date[6]
set INTERVAL_SS       = `head -4 wrf.info | tail -1`
set RUN_HOURS         = `expr $INTERVAL_SS \/ 3600`
set REMAIN            = `expr $INTERVAL_SS \% 3600`
set RUN_MINUTES       = `expr $REMAIN \/ 60`
set RUN_SECONDS       = `expr $REMAIN \% 60`
set WRF_DT            = `head -5 wrf.info | tail -1`
set GRID_DISTANCE     = `head -6 wrf.info | tail -1`
set WEST_EAST_GRIDS   = `head -7 wrf.info | tail -1`
set SOUTH_NORTH_GRIDS = `head -8 wrf.info | tail -1`
set VERTICAL_GRIDS    = `head -9 wrf.info | tail -1`

#-----------------------------------------------------------------------
# Create WRF namelist.input:
#-----------------------------------------------------------------------

cat >! namelist.input << EOF
 &time_control
 run_days                            = 0,
 run_hours                           = ${RUN_HOURS},
 run_minutes                         = ${RUN_MINUTES},
 run_seconds                         = ${RUN_SECONDS},
 start_year                          = ${START_YEAR},  ${START_YEAR},  ${START_YEAR},
 start_month                         = ${START_MONTH}, ${START_MONTH}, ${START_MONTH},
 start_day                           = ${START_DAY},   ${START_DAY},   ${START_DAY},
 start_hour                          = ${START_HOUR},  ${START_HOUR},  ${START_HOUR},
 start_minute                        = ${START_MIN},   ${START_MIN},   ${START_MIN},
 start_second                        = ${START_SEC},   ${START_SEC},   ${START_SEC},
 end_year                            = ${END_YEAR},  ${END_YEAR},  ${END_YEAR},
 end_month                           = ${END_MONTH}, ${END_MONTH}, ${END_MONTH},
 end_day                             = ${END_DAY},   ${END_DAY},   ${END_DAY},
 end_hour                            = ${END_HOUR},  ${END_HOUR},  ${END_HOUR},
 end_minute                          = ${END_MIN},   ${END_MIN},   ${END_MIN},
 end_second                          = ${END_SEC},   ${END_SEC},   ${END_SEC},
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

mv wrfinput wrfinput_d01

# Update boundary conditions

${PBS_O_WORKDIR}/update_wrf_bc >& out.update_wrf_bc

${PBS_O_WORKDIR}/wrf.exe >>& out_wrf_integration

set SUCCESS = `grep "wrf: SUCCESS COMPLETE WRF" out_wrf_integration | cat | wc -l`
if ($SUCCESS != 1) then
   echo $element >> $PBS_O_WORKDIR/blown_${days}_${secs}.out
endif

mv wrfout_d01_* wrfinput

mv dart_wrf_vector dart_wrf_vector.input

# create new input to DART (taken from "wrfinput")
echo ".false." | ${PBS_O_WORKDIR}/dart_tf_wrf >& out.wrf_to_dart

mv dart_wrf_vector $PBS_O_WORKDIR/assim_model_state_ud$element

cd $PBS_O_WORKDIR
#rm -rf $temp_dir

exit
