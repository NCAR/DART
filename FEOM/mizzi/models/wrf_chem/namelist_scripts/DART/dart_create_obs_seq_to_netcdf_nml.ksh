#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &obs_seq_to_netcdf_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &obs_seq_to_netcdf_nml
   obs_sequence_name = ${NL_OBS_SEQUENCE_NAME:-'null'},
   obs_sequence_list = ${NL_OBS_SEQUENCE_LIST:-'obs_coverage_list.txt'},
   append_to_netcdf  = ${NL_APPEND_TO_NETCDF:-.false.},
   lonlim1           = ${NL_LONLIM1:-0.0},
   lonlim2           = ${NL_LONLIM2:-360.0},
   latlim1           = ${NL_LATLIM1:--90.0},
   latlim2           = ${NL_LATLIM2:-90.0},
   verbose           = ${NL_VERBOSE:-.false.},
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
