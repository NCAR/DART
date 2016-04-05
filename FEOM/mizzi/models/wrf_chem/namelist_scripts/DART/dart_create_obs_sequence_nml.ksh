#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &obs_sequence_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &obs_sequence_nml
   write_binary_obs_sequence   = ${NL_WRITE_BINARY_OBS_SEQUENCE:-.false.},
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
