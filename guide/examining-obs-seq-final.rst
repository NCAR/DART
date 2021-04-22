Examining the obs_seq.final file
================================

1. If you are testing with a single observation, just look in the file. If this
   file is in binary format, edit the ``&obs_sequence_nml`` namelist in 
   *input.nml* so the output observation sequence file will be written in
   ASCII:

   .. code-block:: fortran

      &obs_sequence_nml
        write_binary_obs_sequence = .false.
      /
    
   Then rerun *filter* to regenerate an *obs_seq.final* file in ASCII. For 
   an explanation of the contents of your *obs_seq.final* file, see
   :doc:`detailed-structure-obs-seq`.

2. If you are using many observations, run the 
   :doc:`obs_diag <../assimilation_code/programs/readme>` program appropriate
   for your model. The :doc:`matlab-observation-space` will help to summarize
   your output and to explore what is going on.

If there are no changes in the model state after assimilation and a visual examination
of *obs_seq.final* was not informative, convert the *obs_seq.final* file to netCDF with 
:doc:`obs_seq_to_netcdf <../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf>`
and either use the Matlab tools distributed with DART or something of your own. Actually,
*obs_seq_to_netcdf* works on all observation sequence files, not just *obs_seq.final* files.

