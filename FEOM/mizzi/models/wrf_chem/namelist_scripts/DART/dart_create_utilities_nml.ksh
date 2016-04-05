#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &utilities_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &utilities_nml
   TERMLEVEL                   = ${NL_TERMLEVEL:-1},
   module_details              = ${NL_MODULE_DETAILS:-.false.},
   logfilename                 = ${NL_LOGFILENAME:-'dart_log.out'},
   nmlfilename                 = ${NL_NMLFILENAME:-'dart_log.nml'},
   write_nml                   = ${NL_WRITE_NML:-'file'},
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
