#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &replace_wrf_fields_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &replace_wrf_fields_nml
   debug = ${NL_DEBUG:-.false.},
   fail_on_missing_field = ${NL_FAIL_ON_MISSING_FIELD:-.false.},
   fieldnames = ${NL_FIELDNAMES:-'null'},
   fieldlist_file = ${NL_FIELDLIST_FILE:-'null'},
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
