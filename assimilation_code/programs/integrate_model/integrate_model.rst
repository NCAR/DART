PROGRAM ``integrate_model``
===========================

Overview
--------

Generic main program which can be compiled with a model-specific ``model_mod.f90`` file. The model must provide an
``adv_1step()`` subroutine which advances one copy of the model forward in time.

The executable built by this program can be used by the serial program ``perfect_model_obs``, or either the serial or
parallel version of the ``filter`` program. This program is called by the default script in the template directory
called ``advance_model.csh``, and is selected by setting the corresponding ``"async = "`` namelist setting to 2.

This program only advances a single ensemble member per execution and is expected to be run as a serial program. It can
be compiled with the MPI wrappers and called with mpirun with more than 1 task, however, it will only call the model
advance subroutine from a single task (task 0). This can be useful in testing various scripting options using simpler
and smaller models in preparation for running a larger parallel model.

Namelist
--------

There is no namelist for this program.

Modules used
------------

::

   types_mod
   time_manager_mod
   utilities_mod
   assim_model_mod
   obs_model_mod
   ensemble_manager_mod
   mpi_utilities_mod

Files
-----

-  inputfile (temp_ic)
-  outputfile (temp_ud)

References
----------

-  none
