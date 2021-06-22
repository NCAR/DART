program ``create_obs_sequence``
===============================

Overview
--------

This program creates an observation sequence file using values read from standard input. It is typically used to create
synthetic observations, or shorter sequences of observations (although there is no limit on the number of observations).
For creating observation sequence files directly from large, real-world observation datasets, see the
`observations <../../../observations/obs_converters/README.rst>`__ directory.

This program can be run interactively (input from a terminal), or input files can be created with a text editor, perl or
matlab script, or any other convenient method, and then run with standard input redirected from this file. The latter
method is most commonly used to create larger observation sequence files for perfect model applications.

The program can create complete observation sequences ready to be assimilated, or it can create observations with only
partial data which is later filled in by another program. Each observation needs to have a type, location, time,
expected error, and optionally a data value and/or a quality control indicator. For perfect model applications, it is
usually convenient to define 0 quality control fields and 0 copies of the data for each observation. The output of
create_obs_sequence can be read by :doc:`../../../assimilation_code/programs/perfect_model_obs/perfect_model_obs` which
will then create a synthetic (perfect_model) observation sequence complete with two copies of the data for each
observation: the observed value and the 'true' value.

Another common approach for perfect model applications is to use create_obs_sequence to define a set of observation
locations and types, and where observations will be repeatedly sampled in time. When running create_obs_sequence,
specify a single observation for each different location and type, with 0 copies of data and giving all the observations
the same time. Then the program :doc:`../create_fixed_network_seq/create_fixed_network_seq` can read the output of
create_obs_sequence and create an observation sequence file that will contain the set of input observations at a number
of different times. This models a fixed observation station, observing the system at some frequency in time.

This program can also create what are called "identity observations". These are observations located directly at one of
the state variables, so that computing the value requires no model interpolation but simply returns the actual state
variable value. To specify these types of observations, the convention is to put in the negative index number for the
offset of that state variable in the state vector. By specifying the index both the observation kind and location are
defined by the kind and location of that state variable.

The types of observations which can be created by this program is controlled by the observation types built into the
source files created by the :doc:`../../../assimilation_code/programs/preprocess/preprocess` program. The preprocess
namelist sets the available observation types, and must be run each time it is changed, and then the create_obs_sequence
program must be recompiled to incorporate the updated source files.

Other modules used
------------------

::

   utilities_mod
   obs_sequence_mod
   assim_model_mod

Namelist
--------

This program does not use a namelist. All user input is prompted for at the command line.

Files
-----

-  A file containing the output sequence is created.
   (``set_def.out`` is the recommended name)

References
----------

-  none
