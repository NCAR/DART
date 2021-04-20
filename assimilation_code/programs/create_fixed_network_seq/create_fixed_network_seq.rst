program ``create_fixed_network_seq``
====================================

Overview
--------

Reads in an observation sequence file and creates a second observation sequence file. Any time information in the input
file is ignored entirely. All of the observations in the input file define a set of observations. The output sequence
replicates this set multiple times, either with a fixed period in time or at arbitrarily selected times. The program is
driven by input from standard input, either the terminal or a text file.

First, one must select either a regularly repeating time sequence of observations (option 1) or an arbitrarily repeating
sequence (option 2). For the fixed period, the total number of observation times, the first observation time and the
period of the observations is input and an output observation sequence is generated. For the arbitrary period, the user
is queried for the number of observing times and then a set of monotonically increasing times. Finally, the user selects
a file name (traditionally obs_seq.in) to which the output file is written. The format of the output file is controlled
by the namelist options in `obs_sequence_mod <../../modules/observations/obs_sequence_mod.html#Namelist>`__.

Any data values or quality control flags associated with the input set are replicated to the output, but this program is
typically used with perfect model experiments to create observations without data, which are then filled in by running
:doc:`../perfect_model_obs/perfect_model_obs`.

Modules used
------------

::

   types_mod
   utilities_mod
   obs_def_mod
   obs_sequence_mod
   time_manager_mod
   model_mod

Files
-----

-  Input observation sequence (set_def.out is standard).
-  Output observation sequence (obs_seq.in is standard).

References
----------

-  none
