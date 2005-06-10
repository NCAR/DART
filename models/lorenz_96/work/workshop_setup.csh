csh mkmf_preprocess
make
rm -f ../../../obs_def/obs_def_mod.f90
./preprocess
csh mkmf_create_fixed_network_seq
make
csh mkmf_create_obs_sequence
make
csh mkmf_obs_diag
make
csh mkmf_perfect_model_obs
make
csh mkmf_filter
make
./perfect_model_obs
./filter
