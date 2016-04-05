#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &cov_cutoff_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &cov_cutoff_nml
   select_localization             = ${NL_SELECT_LOCALIZATION:-1},
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
