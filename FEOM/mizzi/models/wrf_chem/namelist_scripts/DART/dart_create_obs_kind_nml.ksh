#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &obs_kind_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &obs_kind_nml
   assimilate_these_obs_types = ${NL_ASSIMILATE_THESE_OBS_TYPES:-'null'},
   evaluate_these_obs_types = ${NL_EVALUATE_THESE_OBS_TYPES:-'null'},
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
