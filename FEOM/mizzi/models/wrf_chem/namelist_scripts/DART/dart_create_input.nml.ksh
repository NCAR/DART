#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART input,nml 
#
#########################################################################
#
echo off
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_assim_model_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_assim_tools_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_cov_cutoff_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_dart_to_wrf_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_ensemble_manager_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_filter_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_location_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_model_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_def_radar_mod_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_diag_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_kind_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_selection_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_seq_coverage_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_seq_to_netcdf_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_sequence_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_sequence_tool_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_perfect_model_obs_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_preprocess_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_reg_factor_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_replace_wrf_fields_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_restart_file_tool_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_restart_file_utility_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_schedule_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_smoother_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_utilities_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_wrf_obs_preproc_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_wrf_to_dart_nml.ksh
echo on
