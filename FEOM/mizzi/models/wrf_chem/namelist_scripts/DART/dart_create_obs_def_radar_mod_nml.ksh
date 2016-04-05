#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &obs_def_radar_mod_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &obs_def_radar_mod_nml
   apply_ref_limit_to_obs      = ${NL_APPLY_REF_LIMIT_TO_OBS:-.false.},
   reflectivity_limit_obs      = ${NL_REFLECTIVITY_LIMIT_OBS:--10.0},
   lowest_reflectivity_obs     = ${NL_LOWEST_REFLECTIVITY_OBS:--10.0},
   apply_ref_limit_to_fwd_op   = ${NL_APPLY_REF_LIMIT_TO_FWD_OP:-.false.},
   reflectivity_limit_fwd_op   = ${NL_REFLECTIVITY_LIMIT_FWD_OP:--10.0},
   lowest_reflectivity_fwd_op  = ${NL_LOWEST_REFLECTIVITY_FWD_OP:--10.0},
   max_radial_vel_obs          = ${NL_MAX_RADIAL_VEL_OBS:-1000000},
   allow_wet_graupel           = ${NL_ALLOW_WET_GRAUPEL:-.false.},
   microphysics_type           = ${NL_MICROPHYSICS_TYPE:-2},
   allow_dbztowt_conv          = ${NL_ALLOW_DBZTOWT_CONV:-.false.},
   dielectric_factor           = ${NL_DIELECTRIC_FACTOR:-0.224},
   n0_rain                     = ${NL_N0_RAIN:-8.0e6},
   n0_graupel                  = ${NL_N0_GRAUPEL:-4.0e6},
   n0_snow                     = ${NL_N0_SNOW:-3.0e6},
   rho_rain                    = ${NL_RHO_RAIN:-1000.0},
   rho_graupel                 = ${NL_RHO_GRAUPEL:-400.0},
   rho_snow                    = ${NL_RHO_SNOW:-100.0},
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
