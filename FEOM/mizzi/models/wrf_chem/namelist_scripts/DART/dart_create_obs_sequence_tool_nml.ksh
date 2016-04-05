#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &obs_sequence_tool_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &obs_sequence_tool_nml
   num_input_files          = ${NL_NUM_INPUT_FILES:-1}, 
   filename_seq             = ${NL_FILENAME_SEQ:-'obs_seq.out'},
   filename_out             = ${NL_FILENAME_OUT:-'obs_seq.processed'}, 
   first_obs_days           = ${NL_FIRST_OBS_DAYS:-148816},
   first_obs_seconds        = ${NL_FIRST_OBS_SECONDS:-75601},
   last_obs_days            = ${NL_LAST_OBS_DAYS:-148817},
   last_obs_seconds         = ${NL_LAST_OBS_SECONDS:-10800},
   obs_types                = ${NL_OBS_TYPES:-''}, 
   keep_types               = ${NL_KEEP_TYPES:-.false.}, 
   print_only               = ${NL_PRINT_ONLY:-.false.}, 
   min_lat                  = ${NL_MIN_LAT:--90.0}, 
   max_lat                  = ${NL_MAX_LAT:-90.0}, 
   min_lon                  = ${NL_MIN_LON:-0.0}, 
   max_lon                  = ${NL_MAX_LON:-360.0},
   qc_metadata              = ${NL_QC_METADATA:-''},
   min_qc                   = ${NL_MIN_QC:-0},
   max_qc                   = ${NL_MAX_QC:-5},
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
