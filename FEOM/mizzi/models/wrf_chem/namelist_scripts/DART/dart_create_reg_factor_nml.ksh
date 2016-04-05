#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &reg_factor_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &reg_factor_nml
   select_regression           = ${NL_SELECT_REGRESSION:-1},
   input_reg_file              = ${NL_INPUT_REG_FILE:-'time_mean_reg'},
   save_reg_diagnostics        = ${NL_SAVE_REG_DIAGNOSTICS:-.false.},
   reg_diagnostics_file        = ${NL_REG_DIAGNOSTICS_FILE:-'reg_diagnostics'},
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
