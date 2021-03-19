Examining the obs_seq.final file
================================

If there are no changes in the model state after assimilation, then examine the
``obs_seq.final`` file. There are two ways to do this.

1. If you are testing with a single observation, just look in the file. If this
   file is in binary format, edit the ``&obs_sequence_nml`` namelist in 
   ``input.nml`` so the output observation sequence file will be written in
   ASCII:

   .. code-block:: fortran

      &obs_sequence_nml
        write_binary_obs_sequence = .false.
      /
    
   Then rerun ``filter`` to regenerate an ``obs_seq.final`` file in ASCII. For 
   an explanation of the contents of your ``obs_seq.final`` file, see
   :doc:`obs-seq-file`.

2. If you are using many observations, run the ``obs_diag`` program appropriate
   for your model. The :doc:`matlab-observation-space` will help to summarize
   your output and to explore what is going on.
