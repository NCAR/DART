## To convert coamps hdf5 files to netCDF
   - run *trans_coamps_to_dart*
   - requires a `convert.nml`  that specifies the cdtg for the variable names.
   - requires a `state.vars` file (that can be created from *populate_state_vars.sh* and a template `state.dat` file. This requires knowledge of the number of layers.
   - input hdf5 filename is hardcoded to 'coamps.hdf5'
   - This creates a `dart_vector.nc` file that has the Exner functions, the grid information, and the required model variables.

## To create synthetic observations
   - run *perfect_model_obs*
   - edit the `input.nml:model_nml:cdtg` to match file

## To convert innovation vectors to observation sequence files
   - run *innov_to_obs_seq*
