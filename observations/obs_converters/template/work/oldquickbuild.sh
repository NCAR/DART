#!/usr/bin/env bash
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

main() {

  export DART=$(git rev-parse --show-toplevel)

  # Source build functions
  source "$DART/build_templates/buildconvfunctions.sh"

  # Define converter and location
  CONVERTER=$(basename "$(dirname "$PWD")")
  LOCATION=threed_sphere  # Change this as needed

  # Clean directory
  \rm -f -- *.o *.mod Makefile .cppdefs

  # Build and run preprocess before making any other DART executables
  buildpreprocess

  buildconv
  
  # Find the Fortran source file in the parent directory
  src_file="../${CONVERTER}_to_obs.f90"
  if [[ ! -f "$src_file" ]]; then
    echo "No Fortran source file found at $src_file."
    exit 1
  fi

  MOD_PATH="$DART/assimilation_code/modules"

  sed -i '/if ( num_types_evaluate > 0 ) use_precomputed_FOs_these_obs_types(j:num_types_use_precomputed_FOs) = evaluate_these_obs_types(1:num_types_evaluate)/{
  s/if ( num_types_evaluate > 0 ) use_precomputed_FOs_these_obs_types(j:num_types_use_precomputed_FOs) = evaluate_these_obs_types(1:num_types_evaluate)/if ( num_types_evaluate > 0 ) then/
  a\
  use_precomputed_FOs_these_obs_types(j:num_types_use_precomputed_FOs) = \&\
    evaluate_these_obs_types(1:num_types_evaluate)\
  endif
}' "$MOD_PATH/observations/obs_kind_mod.f90"

  # Set paths
  OBS_PATH="$DART/observations/forward_operators"
  UTIL_PATH="$DART/assimilation_code/modules/utilities"
  LOC_PATH="$DART/assimilation_code/location/$LOCATION"
  OBS_CONV_UTIL_PATH="$DART/observations/obs_converters/utilities"
  OBS_KIND_PATH="$MOD_PATH/observations"

  # Create basic Makefile
  echo "FC = gfortran" > Makefile
  echo "FFLAGS = -O2 -g -Wno-error=line-truncation -I$MOD_PATH -I$OBS_PATH -I$UTIL_PATH -I$LOC_PATH -I$OBS_CONV_UTIL_PATH -I$MOD_PATH/observations" >> Makefile
  echo "" >> Makefile
  echo "OBJS = types_mod.o utilities_mod.o time_manager_mod.o random_seq_mod.o location_mod.o obs_sequence_mod.o obs_utilities_mod.o obs_kind_mod.o netcdf_utilities_mod.o mpi_utilities_mod.o ../${CONVERTER}_to_obs.o" >> Makefile
  echo "" >> Makefile
  echo "all: $CONVERTER" >> Makefile
  echo "" >> Makefile
  echo "$CONVERTER: \$(OBJS)" >> Makefile
  echo "	\$(FC) \$(FFLAGS) -o $CONVERTER \$(OBJS)" >> Makefile
  echo "" >> Makefile
  echo "types_mod.o: $UTIL_PATH/types_mod.f90" >> Makefile
  echo "	\$(FC) \$(FFLAGS) -c $UTIL_PATH/types_mod.f90" >> Makefile
  echo "" >> Makefile
  echo "utilities_mod.o: $UTIL_PATH/utilities_mod.f90 types_mod.o" >> Makefile
  echo "	\$(FC) \$(FFLAGS) -c $UTIL_PATH/utilities_mod.f90" >> Makefile
  echo "" >> Makefile
  echo "time_manager_mod.o: $UTIL_PATH/time_manager_mod.f90 types_mod.o utilities_mod.o" >> Makefile
  echo "	\$(FC) \$(FFLAGS) -c $UTIL_PATH/time_manager_mod.f90" >> Makefile
  echo "" >> Makefile
  echo "random_seq_mod.o: $UTIL_PATH/random_seq_mod.f90 types_mod.o utilities_mod.o" >> Makefile
  echo "	\$(FC) \$(FFLAGS) -c $UTIL_PATH/random_seq_mod.f90" >> Makefile
  echo "" >> Makefile
  echo "mpi_utilities_mod.o: $UTIL_PATH/mpi_utilities_mod.f90 types_mod.o utilities_mod.o" >> Makefile
  echo "	\$(FC) \$(FFLAGS) -c $UTIL_PATH/mpi_utilities_mod.f90" >> Makefile
  echo "" >> Makefile
  echo "location_mod.o: $LOC_PATH/location_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o random_seq_mod.o mpi_utilities_mod.o" >> Makefile
  echo "	\$(FC) \$(FFLAGS) -c $LOC_PATH/location_mod.f90" >> Makefile
  echo "" >> Makefile
  echo "obs_sequence_mod.o: $MOD_PATH/observations/obs_sequence_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o location_mod.o random_seq_mod.o mpi_utilities_mod.o" >> Makefile
  echo "	\$(FC) \$(FFLAGS) -c $MOD_PATH/observations/obs_sequence_mod.f90" >> Makefile
  echo "" >> Makefile
  echo "obs_utilities_mod.o: $OBS_CONV_UTIL_PATH/obs_utilities_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o location_mod.o obs_sequence_mod.o random_seq_mod.o mpi_utilities_mod.o" >> Makefile
  echo "	\$(FC) \$(FFLAGS) -c $OBS_CONV_UTIL_PATH/obs_utilities_mod.f90" >> Makefile
  echo "" >> Makefile
  echo "obs_kind_mod.o: $MOD_PATH/observations/obs_kind_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o location_mod.o obs_sequence_mod.o obs_utilities_mod.o random_seq_mod.o mpi_utilities_mod.o" >> Makefile
  echo "	\$(FC) \$(FFLAGS) -c $MOD_PATH/observations/obs_kind_mod.f90" >> Makefile
  echo "" >> Makefile
  echo "netcdf_utilities_mod.o: $UTIL_PATH/netcdf_utilities_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o location_mod.o obs_sequence_mod.o obs_utilities_mod.o obs_kind_mod.o random_seq_mod.o mpi_utilities_mod.o" >> Makefile
  echo "	\$(FC) \$(FFLAGS) -c $UTIL_PATH/netcdf_utilities_mod.f90" >> Makefile
  echo "" >> Makefile
  echo "../${CONVERTER}_to_obs.o: ../${CONVERTER}_to_obs.f90 utilities_mod.o types_mod.o time_manager_mod.o location_mod.o obs_sequence_mod.o obs_utilities_mod.o obs_kind_mod.o netcdf_utilities_mod.o random_seq_mod.o mpi_utilities_mod.o" >> Makefile
  echo "	\$(FC) \$(FFLAGS) -c ../${CONVERTER}_to_obs.f90" >> Makefile
  echo "" >> Makefile
  echo "clean:" >> Makefile
  echo "	rm -f *.o *.mod $CONVERTER" >> Makefile
  echo "" >> Makefile
  echo ".PHONY: all clean" >> Makefile

  echo "Generated Makefile:"
  cat Makefile

  # Build converter
  make

  # Run converter
  if [[ -f "$CONVERTER" ]]; then
    ./$CONVERTER
  else
    echo "Build failed. Converter executable not found."
    exit 1
  fi

  # Clean up
  \rm -f -- *.o *.mod

  echo -e "Converter built successfully!\nTo get started, check the output files in the current directory."
}

main "$@"
