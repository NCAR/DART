#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &ensemble_manager_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &ensemble_manager_nml
   single_restart_file_in  = ${NL_SINGLE_RESTART_FILE_IN:-.true.},
   single_restart_file_out = ${NL_SINGLE_RESTART_FILE_OUT:-.true.},
   perturbation_amplitude  = ${NL_PERTURBATION_AMPLITUDE:-0.2},
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

