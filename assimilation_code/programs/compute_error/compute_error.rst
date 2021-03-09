PROGRAM ``compute_error``
=========================

Overview
--------

Utility program to compute the time-mean ensemble error and spread in the same manner that the DART MATLAB diagnostic
routine 'plot_total_err' does. It runs from the command line, opens no windows, and outputs several types of numerical
results on standard output. Grep for 'Total' to get the 2 lines with total error and total spread. Intended for scripts
where only the numeric results are wanted instead of a time-series plot. This routine does not do any weighted
computations.

The default is to compare a True_State.nc file output from perfect_model_obs to a Prior_Diag.nc file output from filter.
Other filenames can be specified in the namelist. These files must have at least one overlapping value in the 'time'
array. The statistics will be done on the overlapping time region only.

The output includes the min and max error and spread values, and the time index and time value where that occurs. There
is also an option to recompute the time mean ensemble error and spread after skipping the first N times. This can be
useful to skip an initial error spike while the model is spinning up which can result in a larger than expected total
error.

Namelist interface ``&compute_error_nml`` is read from file ``input.nml``.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &compute_error_nml
      truth_file_name   = 'true_state.nc'
      diag_file_name    = 'preassim.nc'
      skip_first_ntimes = 0
     /

| 

.. container::

   +-------------------+--------------------+---------------------------------------------------------------------------+
   | Item              | Type               | Description                                                               |
   +===================+====================+===========================================================================+
   | truth_file_name   | character(len=256) | State-space diagnostic file from the 'perfect_model_obs' program.         |
   +-------------------+--------------------+---------------------------------------------------------------------------+
   | diag_file_name    | character(len=256) | State space diagnostic file output from the 'filter' program.             |
   +-------------------+--------------------+---------------------------------------------------------------------------+
   | skip_first_ntimes | integer            | If set to a value greater than 0, the error values will be recomputed a   |
   |                   |                    | second time, skipping the first N times. This can be useful when running  |
   |                   |                    | an experiment that has an initial error spike as the model spins up and   |
   |                   |                    | then decays down to a more steady state.                                  |
   +-------------------+--------------------+---------------------------------------------------------------------------+

| 

Modules used
------------

::

   types_mod
   utilities_mod

Files
-----

-  DART diagnosic files (True_State.nc, Prior_Diag.nc)
-  compute_error.nml

References
----------

-  none
