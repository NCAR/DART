#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &filter_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &filter_nml
   async                        = ${NL_ASYNC:-2},
   adv_ens_command              = ${NL_ADV_ENS_COMMAND:-'./advance_model.csh'},
   ens_size                     = ${NL_ENS_SIZE:-3},
   start_from_restart           = ${NL_START_FROM_RESTART:-.true.},
   output_restart               = ${NL_OUTPUT_RESTART:-'restart.out'},
   obs_sequence_in_name         = ${NL_OBS_SEQUENCE_IN_NAME:-'obs_seq.out'},
   obs_sequence_out_name        = ${NL_OBS_SEQUENCE_OUT_NAME:-'obs_seq.final'},
   restart_in_file_name         = ${NL_RESTART_IN_FILE_NAME:-'filter_ic_prior'},
   restart_out_file_name        = ${NL_RESTART_OUT_FILE_NAME:-'filter_ic_post'},
   init_time_days               = ${NL_INIT_TIME_DAYS:--1},
   init_time_seconds            = ${NL_INIT_TIME_SECONDS:--1},
   first_obs_days               = ${NL_FIRST_OBS_DAYS:-148816},
   first_obs_seconds            = ${NL_FIRST_OBS_SECONDS:-75601},
   last_obs_days                = ${NL_LAST_OBS_DAYS:-148817},
   last_obs_seconds             = ${NL_LAST_OBS_SECONDS:-10800},
   num_output_state_members     = ${NL_NUM_OUTPUT_STATE_MEMBERS:-0},
   num_output_obs_members       = ${NL_NUM_OUTPUT_OBS_MEMBERS:-0},
   output_interval              = ${NL_OUTPUT_INTERVAL:-1},
   num_groups                   = ${NL_NUM_GROUPS:-1},
   input_qc_threshold           = ${NL_INPUT_QC_THRESHOLD:-4.0},
   outlier_threshold            = ${NL_OUTLIER_THRESHOLD:-3.0},
   enable_special_outlier_code  = ${NL_ENABLE_SPECIAL_OUTLIER_CODE:-.false.},
   output_forward_op_errors     = ${NL_OUTPUT_FORWARD_OP_ERRORS:-.false.},
   output_timestamps            = ${NL_PUTPUT_TIMESTAMPS:-.true.},
   output_inflation             = ${NL_OUTPUT_INFLATION:-.true.},
   trace_execution              = ${NL_TRACE_EXECUTION:-.true.},
   inf_flavor                   = ${NL_INF_FLAVOR_PRIOR:-0}, ${NL_INF_FLAVOR_POST:-0},
   inf_initial_from_restart     = ${NL_INF_INITIAL_FROM_RESTART_PRIOR:-.false.}, ${NL_INF_INITIAL_FROM_RESTART_POST:-.false.},
   inf_sd_initial_from_restart  = ${NL_INF_SD_INITIAL_FROM_RESTART_PRIOR:-.false.}, ${NL_INF_SD_INITIAL_FROM_RESTART_POST:-.false.},
   inf_output_restart           = ${NL_INF_OUTPUT_RESTART_PRIOR:-.true.}, ${NL_INF_OUTPUT_RESTART_POST:-.true.},
   inf_deterministic            = ${NL_INF_DETERMINISTIC_PRIOR:-.true.}, ${NL_INF_DETERMINISTIC_POST:-.true.},
   inf_in_file_name             = ${NL_INF_IN_FILE_NAME_PRIOR:-'prior_inf_ic_old'}, ${NL_INF_IN_FILE_NAME_POST:-'post_inf_ic_old'},
   inf_out_file_name            = ${NL_INF_OUT_FILE_NAME_PRIOR:-'prior_inf_ic_new'}, ${NL_INF_OUT_FILE_NAME_POST:-'post_inf_ic_new'},
   inf_diag_file_name           = ${NL_INF_DIAG_FILE_NAME_PRIOR:-'prior_inf_diag'}, ${NL_DIAG_FILE_NAME_POST:-'post_inf_diag'},
   inf_initial                  = ${NL_INF_INITIAL_PRIOR:-1.00}, ${NL_INF_INITIAL_POST:-1.00},
   inf_sd_initial               = ${NL_INF_SD_INITIAL_PRIOR:-0.60}, ${NL_INF_SD_INITIAL_POST:-0.0},
   inf_damping                  = ${NL_INF_DAMPING_PRIOR:-0.90}, ${NL_INF_DAMPING_POST:-1.00},
   inf_lower_bound              = ${NL_INF_LOWER_BOUND_PRIOR:-1.0}, ${NL_INF_LOWER_BOUND_POST:-1.0},
   inf_upper_bound              = ${NL_INF_UPPER_BOUND_PRIOR:-1000000.0}, ${NL_INF_UPPER_BOUND_POST:-1000000.0},
   inf_sd_lower_bound           = ${NL_INF_SD_LOWER_BOUND_PRIOR:-0.60}, ${NL_INF_SD_LOWER_BOUND_POST:-0.0},
/
EOF
#
# Append namelist section to input.nml
if [[ -f input.nml ]]; then
   cat input.nml_temp >> input.nml
   rm input.nml_temp
else
   mv input.nml_temp input.nml
fi


