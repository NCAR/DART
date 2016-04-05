#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &perfect_model_obs_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
&perfect_model_obs_nml
   start_from_restart           = ${NL_START_FROM_RESTART:-.true.},
   output_restart               = ${NL_OUTPUT_RESTART:-.true.},
   async                        = ${NL_ASYNC:-2},
   init_time_days               = ${NL_INIT_TIME_DAYS:--1},
   init_time_seconds            = ${NL_INIT_TIME_SECONDS:--1},
   first_obs_days               = ${NL_FIRST_OBS_DAYS:-148816},
   first_obs_seconds            = ${NL_FIRST_OBS_SECONDS:-75601},
   last_obs_days                = ${NL_LAST_OBS_DAYS:-148817},
   last_obs_seconds             = ${NL_LAST_OBS_SECONDS:-10800},
   output_interval              = ${NL_OUTPUT_INTERVAL:-1},
   restart_in_file_name         = ${NL_RESTART_IN_FILE_NAME:-'perfect_ics'},
   restart_out_file_name        = ${NL_RESTART_OUT_FILE_NAME:-'perfect_restart'},
   obs_seq_in_file_name         = ${NL_OBS_SEQ_IN_FILE_NAME:-'obs_seq.in'},
   obs_seq_out_file_name        = ${NL_OBS_SEQ_OUT_FILE_NAME:-'obs_seq.out'},
   adv_ens_command              = ${NL_ADV_ENS_COMMAND:-'./advance_model.csh'},
   output_timestamps            = ${NL_OUTPUT_TIMESTAMPS:-.false.},
   trace_execution              = ${NL_TRACE_EXECUTION:-.false.},
   output_forward_op_errors     = ${NL_OUTPUT_FORWARD_OP_ERRORS:-.false.},
   print_every_nth_obs          = ${NL_PRINT_EVERY_NTH_OBS:--1},
   silence                      = ${NL_SILENCE:-.false.},
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


