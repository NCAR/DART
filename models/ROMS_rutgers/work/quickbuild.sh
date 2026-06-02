#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

main() {

export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildfunctions.sh

MODEL=ROMS_rutgers
LOCATION=threed_sphere


programs=(
closest_member_tool
filter
model_mod_check
perfect_model_obs
perturb_single_instance
wakeup_filter
)

serial_programs=(
advance_time
create_fixed_network_seq
create_obs_sequence
fill_inflation_restart
obs_common_subset
obs_diag
obs_selection
obs_seq_coverage
obs_seq_to_netcdf
obs_seq_verify
obs_sequence_tool
)

# Parse reg_grid/irreg_grid option and remove it from the argument list.
# Usage: quickbuild.sh [reg_grid|irreg_grid] [mpi/nompi/mpif08] [program]
#   reg_grid   - link model_mod_reg_grid   to ../chosen_model_mod.f90 before building
#   irreg_grid - link model_mod_irreg_grid to ../chosen_model_mod.f90 before building
grid_type=""
remaining_args=()
for arg in "$@"; do
    case $arg in
        reg_grid|irreg_grid)
            grid_type=$arg
            ;;
        *)
            remaining_args+=("$arg")
            ;;
    esac
done

if [[ "$grid_type" == "reg_grid" ]]; then
    cp -f "$DART/models/ROMS_rutgers/model_mod_reg_grid" "$DART/models/ROMS_rutgers/model_mod.f90"
elif [[ "$grid_type" == "irreg_grid" ]]; then
    cp -f "$DART/models/ROMS_rutgers/model_mod_irreg_grid" "$DART/models/ROMS_rutgers/model_mod.f90"
else
    echo "ERROR: grid type required."
    echo "Usage: quickbuild.sh reg_grid|irreg_grid [mpi/nompi/mpif08] [program]"
    exit 1
fi

arguments "${remaining_args[@]}"

# clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# build and run preprocess before making any other DART executables
buildpreprocess

# build DART
buildit

# clean up
\rm -f -- *.o *.mod

}

main "$@"
