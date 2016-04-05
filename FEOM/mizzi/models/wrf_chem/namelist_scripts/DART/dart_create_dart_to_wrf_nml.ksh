#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &dart_to_wrf_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &dart_to_wrf_nml
   model_advance_file           = ${NL_MODEL_ADVANCE_FILE:-.true.},
   adv_mod_command              = ${NL_ADV_MOD_COMMAND:-'wrf.exe'},
   dart_restart_name            = ${NL_DART_RESTART_NAME:-'dart_wrf_vector'},
   print_data_ranges            = ${NL_PRINT_DATA_RANGES:-.false.},
   debug                        = ${NL_DEBUG:-.false.},
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
