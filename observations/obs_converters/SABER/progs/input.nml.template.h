
&preprocess_nml
    input_obs_kind_mod_file = '../../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90',
   output_obs_kind_mod_file = '../../../../assimilation_code/modules/observations/obs_kind_mod.f90',
     input_obs_def_mod_file = '../../../../observations/forward_operators/DEFAULT_obs_def_mod.F90',
    output_obs_def_mod_file = '../../../../observations/forward_operators/obs_def_mod.f90',
   input_files              = '../../../../observations/forward_operators/obs_def_ocean_mod.f90' /

&obs_kind_nml
 /

&obs_def_gps_nml
 /

&location_nml
 /

&utilities_nml
 module_details = .false.,
 write_nml = 'none',
 /

&obs_sequence_nml
   write_binary_obs_sequence = .false.  /

&obs_sequence_tool_nml
   filename_seq      = '',
   filename_seq_list = 'olist',
   filename_out      = 'OUTDIR/obs_seqYYYYMMDDHH',
   print_only        = .false.,
   gregorian_cal     = .true.,
   first_obs_days    = GREG0,
   first_obs_seconds = SECS0,
   last_obs_days     = GREG1,
   last_obs_seconds  = SECS1,
/

