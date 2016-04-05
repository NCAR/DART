#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &restart_file_utility_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &restart_file_utility_nml
   input_file_name              = ${NL_INPUT_FILE_NAME:-'filter_restart_in'},
   output_file_name             = ${NL_OUTPUT_FILE_NAME:-'filter_restart_out'},
   ens_size                     = ${NL_ENS_SIZE:-1},
   single_restart_file_in       = ${NL_SINGLE_RESTART_FILE_IN:-.true.},
   single_restart_file_out      = ${NL_SINGLE_RESTART_FILE_OUT:-.true.},
   write_binary_restart_files   = ${NL_WRITE_BINARY_RESTART_FILE:-.false.},
   overwrite_data_time          = ${NL_OVERWRITE_DATA_TIME:-.false.},
   new_data_days                = ${NL_NEW_DATA_DAYS:--1},
   new_data_secs                = ${NL_NEW_DATA_SECS:--1},
   input_is_model_advance_file  = ${NL_INPUT_IS_MODEL_ADVANCE_FILE:-.false.},
   output_is_model_advance_file = ${NL_OUTPUT_IS_MODEL_ADVANCE_FILE:-.false.},,
   overwrite_advance_time       = ${NL_OVERWRITE_ADVANCE_TIME:-.false.},
   new_advance_days             = ${NL_NEW_ADVANCE_DAYS:--1},
   new_advance_secs             = ${NL_NEW_ADVANCE_SECS:--1}, 
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
