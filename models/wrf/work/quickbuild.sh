#!/bin/bash

main() {
set -e

[ -z "$DART" ] && echo "ERROR: Must set DART environment variable" && exit 9

MODEL=wrf
LOCATION=threed_sphere
source $DART/build_templates/buildfunctions.sh

# clean the directory
\rm -f *.o *.mod Makefile .cppdefs

programs=( \
closest_member_tool \
filter \
model_mod_check \
perfect_model_obs \
perturb_single_instance \
wakeup_filter
)

serial_programs=( \
advance_time \
create_fixed_network_seq \
create_obs_sequence \
fill_inflation_restart \
obs_common_subset \
obs_diag \
obs_selection \
obs_seq_coverage \
obs_seq_to_netcdf \
obs_seq_verify \
obs_sequence_tool
)


#radiance_obs_to_netcdf \  # needs rttov

model_serial_programs=(
add_pert_where_high_refl \
advance_cymdh \
convertdate \
ensemble_init \
pert_wrf_bc \
replace_wrf_fields \
select \
update_wrf_bc \
wrf_dart_obs_preprocess
)

arguments "$@"

# clean the directory
\rm -f *.o *.mod Makefile .cppdefs

# build and run preprocess before making any other DART executables
buildpreprocess

# build DART
buildit

# clean up
\rm -f *.o *.mod

}

main "$@"
