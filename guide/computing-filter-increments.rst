Computing filter increments
===========================

.. note::

   This document is written as if your experiment was run with
   ``single_file_out = .true.``. The potential permutations of filenames output
   by filter is enormous, so it isn't feasible to write documentation for all
   possible cases.

After *filter* executes without error and produces an ``obs_seq.final`` file, a
``preassim.nc`` file, and an ``analysis.nc`` file, the first questions to ask
are:

1. Is the model state output from ``filter`` different from the input?
2. Were any observations successfully assimilated?

You can check check if the output model state data was changed by the
assimilation by using the ``ncdiff`` tool to create a file containing the
difference of the ``preassim.nc`` and ``analysis.nc`` files. If you are running
with ``single_file_in = .true.`` and ``single_file_out = .true.`` use
``ncdiff`` on the files output for the analysis and preassim stages:

.. code-block::

   $ ncdiff analysis.nc preassim.nc increments.nc
  
Otherwise, if you are running with ``single_file_in = .false.`` and
``single_file_out = .false.``, use ``ncdiff`` on the ensemble mean files for
the analysis and preassim stages:

.. code-block::

   $ ncdiff analysis_mean.nc preassim_mean.nc increments.nc

``ncdiff`` generates a file, ``increments.nc``, that contains the increments,
or innovations, created by ``filter``. You can view the increments using
``ncview``:

.. code-block::

   $ ncview increments.nc

to examine the ensemble mean variables. If all values are 0, then the
assimilation changed nothing in the state.
