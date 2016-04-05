#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &assim_model_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &assim_model_nml
   write_binary_restart_files      = ${NL_WRITE_BINARY_RESTART_FILE:-.false.},
   netCDF_large_file_support       = ${NL_NETCDF_LARGE_FILE_SUPPORT:-.true.},
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
