## To convert coamps hdf5 files to netCDF
   - run *trans_coamps_to_dart* (which is unfathomably slow ...)
   - requires a `convert.nml`  that specifies the cdtg for the variable names.
   - requires a `state.vars` file (that can be created from *populate_state_vars.sh* and a template `state.dat` file. This requires knowledge of the number of layers.
   - input hdf5 filename is hardcoded to 'coamps.hdf5'
   - This creates a `dart_vector.nc` file that has the Exner functions, the grid information, and the required model variables.

## To create the initial state-space inflation files
   - edit the `input.nml:fill_inflation_restart_nml` to desired values
   - run *fill_inflation_restart*
   - This will create inflation files for each domain
     * `input_priorinf_mean_d0[1,2].nc` and
     * `input_priorinf_sd_d0[1,2].nc`
   - given the logic of the scripting, these should be renamed from 'input' to 'output' 
   - TJH: inflation files to not appear to have the THBM_g0x variables ... why?

## To create synthetic observations
   - run *perfect_model_obs*
   - edit the `input.nml:model_nml:cdtg` to match file

## To convert innovation vectors to observation sequence files
   - link the `innov_a<CDTG>` to `innov.out`, 
   - link the `ngt<CDTG>` to `ngt.out`, 
   - run *innov_to_obs_seq*
   - make sure the `obs_seq.out` has a filename like `obs_seq_<CDTG>.out`

## To run an assimilation on an existing ensemble without advancing the model
   - configure and run *shell_scripts/stage_experiment.csh*
   - cd to the EXPERIMENTDIR ... check things out
   - make sure inflation files are `output_[prior,poste]inf_[mean,sd].nc`
   - eddy:  sbatch run_filter.csh

## problems so far:
   - the inflation files do not have all the variables
   - no innovation observations for the time and domain of the ensemble
