#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &model_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &model_nml
   default_state_variables     = ${NL_DEFAULT_STATE_VARIABLES:-.false.},
   wrf_state_variables         = ${NL_WRF_STATE_VARIABLES:-"null"}
   wrf_state_bounds            = ${NL_WRF_STATE_BOUNDS:-"null"}
   output_state_vector         = ${NL_OUTPUT_STATE_VECTOR:-.false.},
   num_domains                 = ${NL_NUM_DOMAINS:-1},
   calendar_type               = ${NL_CALENDAR_TYPE:-3},
   assimilation_period_seconds = ${NL_ASSIMILATION_PERIOD_SECONDS:-21600},
   vert_localization_coord     = ${NL_VERT_LOCALIZATION_COORD:-3},
   center_search_half_length   = ${NL_CENTER_SEARCH_HALF_LENGTH:-500000.0},
   center_spline_grid_scale    = ${NL_CENTER_SPLINE_GRID_SCALE:-10},   
   sfc_elev_max_diff           = ${NL_SFC_ELEV_MAX_DIFF:-100.0},
   circulation_pres_level      = ${NL_CIRCULATION_PRES_LEVEL:-80000.0},
   circulation_radius          = ${NL_CIRCULATION_RADIUS:-108000.0},
   allow_obs_below_vol         = ${NL_ALLOW_OBS_BELOW_VOL:-.false.},
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
