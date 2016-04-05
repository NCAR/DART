#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &preprocess_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &preprocess_nml
   input_obs_kind_mod_file  = ${NL_INPUT_OBS_KIND_MOD_FILE:-"null"},
   output_obs_kind_mod_file = ${NL_OUTPUT_OBS_KIND_MOD_FILE:-"null"},
   input_obs_def_mod_file   = ${NL_INPUT_OBS_DEF_MOD_FILE:-"null"},
   output_obs_def_mod_file  = ${NL_OUTPUT_OBS_DEF_MOD_FILE:-"null"},
   input_files              = ${NL_INPUT_FILES:-"null"},
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
