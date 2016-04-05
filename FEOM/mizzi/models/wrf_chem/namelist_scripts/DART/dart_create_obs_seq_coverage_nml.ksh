#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &obs_seq_coverage_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &obs_seq_coverage_nml
   obs_sequence_list = ${NL_OBS_SEQUENCE_LIST:-'obs_coverage_list.txt'},
   obs_sequence_name = ${NL_OBS_SEQEUNCE_NAME:-'null'},
   obs_of_interest   = ${NL_OBS_OF_INTEREST:-'null'},
   textfile_out      = ${NL_TEXTFILE_OUT:-'obsdef_mask.txt'},
   netcdf_out        = ${NL_NETCDF_OUT:-'obsdef_mask.nc'},
   lonlim1           = ${NL_LONLIM1:-0.0},
   lonlim2           = ${NL_LONLIM2:-360.0},
   latlim1           = ${NL_LATLIM1:--90.0},
   latlim2           = ${NL_LATLIM2:-90.0},
   nTmin             = ${NL_NTMIN:-8},
   nTmax             = ${NL_NTMAX:-8},
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
