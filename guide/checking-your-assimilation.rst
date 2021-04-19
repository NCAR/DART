Checking your initial assimilation
==================================

You may require several attempts to get your assimilation configured correctly.
The next section, :doc:`computing-filter-increments`, describes how to take the
difference between two assimilation stages to determine whether your initial
assimilation worked as intented.

If your assimilation does not change anything in the model state, you may need
to rerun ``filter`` multiple times to understand what is wrong.

Thus you should make ``filter`` very fast to run. You can do this by:

1. Making an observation sequence file containing a single observation.
2. Configuring your run so that filter does a single assimilation and exits
   without having to advance the ensemble of models or do other work.
  
Making an observation sequence file containing a single observation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can use one of these methods to make an ``obs_seq`` with just a single
observation:

1. Run ``create_obs_sequence`` to make a new, short, observation sequence file.
2. Use the ``obs_sequence_tool`` to cut an existing ``obs_seq.out`` file down
   to just a few obs by selecting only a subset of the types and setting a very
   short time window, such as a second or two when you know there are
   observations available.

These programs are described in the 
:doc:`Programs directory <../assimilation_code/programs/readme>`.

Configuring your run so that filter does a single assimilation and exits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To configure ``filter`` to only do a single assimilation:

1. Edit the ``&filter_nml`` namelist in ``input.nml`` to set the
   ``init_time_days`` and ``init_time_seconds`` to match the observation time
   in your truncated observation sequence file. This overrides any times in the
   input files and ensures that ``filter`` will only assimilate and not try to
   advance the model.
2. Make sure the truncated observation sequence file contains only a single 
   observation or observations close enough together in time to fit into a
   single assimilation window.
