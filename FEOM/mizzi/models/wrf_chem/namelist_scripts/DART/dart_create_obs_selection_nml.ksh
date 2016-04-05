#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &obs_selection_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &obs_selection_nml
   filename_seq          = ${NL_FILENAME_SEQ:-'obs_seq.out'},
   filename_seq_list     = ${NL_FILENAME_SEQ_LIST:-'null'},
   filename_out          = ${NL_FILENAME_OUT:-'obs_seq.processed'},
   selections_file       = ${NL_SELECTIONS_FILE:-'obsdef_mask.txt'},
   selections_is_obs_seq = ${NL_SELECTIONS_IS_OBS_SEQ:-.false.},
   print_only            = ${NL_PRINT_ONLY:-.false.},
   calendar              = ${NL_CALENDAR:-'gregorian'}, 
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
