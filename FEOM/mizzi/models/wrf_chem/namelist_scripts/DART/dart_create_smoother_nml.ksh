#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &smoother_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &smoother_nml
   num_lags                  = ${NL_NUM_LAGS:-0},
   start_from_restart        = ${NL_START_FROM_RESTART:-.false.},
   output_restart            = ${NL_OUTPUT_RESTART:-.false.},
   restart_in_file_name      = ${NL_RESTART_IN_FILE_NAME:-'smoother_ics'},
   restart_out_file_name     = ${NL_RESTART_OUT_FILE_NAME:-'smoother_restart'},
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
