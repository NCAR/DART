#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &wrf_to_dart_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &wrf_to_dart_nml
   dart_restart_name   = ${NL_DART_RESTART_NAME:-'null'},
   print_data_ranges   = ${NL_PRINT_DATA_RANGES:-.false.},
   debug               = ${NL_DEBUG:-.false.},
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


