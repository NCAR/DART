system simulation programs
==========================

Overview
--------

A collection of standalone programs for simulating various properties of ensembles.

-  ``gen_sampling_err_table.f90``
-  ``full_error.f90``
-  ``obs_sampling_err.f90``
-  ``sampling_error.f90``
-  ``system_simulation.f90``
-  ``test_sampling_err_table.f90``
-  ``correl_error.f90``

**The program of most interest here is ``gen_sampling_err_table.f90`` which generates the lookup table needed when using
sampling error correction in ``filter``.** Talk to Jeff Anderson about the other programs in this directory.

To enable the sampling error correction algorithm in ``filter``, set the namelist item `&assim_tools_nml :
sampling_error_correction <../../modules/assimilation/assim_tools_mod.html#Namelist>`__ to *.true.*, and copy the netCDF
file system_simulation/sampling_error_correction_table.nc into the run directory.
The supported set of precomputed ensemble sizes can be found by exploring the ``ens_sizes`` variable in
sampling_error_correction_table.nc. To add support for another ensemble size, build the executables in the
`work <../system_simulation/work>`__ directory, (usually by running ``quickbuild.csh``) set the ``ens_sizes`` (it takes
a list, but keep it short) namelist item in ``work/input.nml``, and run ``gen_sampling_err_table``. It generates a LARGE
number of samples *per ensemble size* for statistical rigor. Larger ensemble sizes take longer to generate, and compiler
optimizations matter - perhaps significantly. For example, the numbers below come from calculating one ensemble size at
a time on my desktop machine with gfortran and basic optimization:

============= ==================
ensemble size run-time (seconds)
============= ==================
10            57
50            273
100           548
============= ==================

The basic structure of sampling_error_correction_table.nc is shown below.

.. container::

   ::

      0[1095] desktop:system_simulation/work % ncdump -v ens_sizes *nc
      netcdf sampling_error_correction_table {
      dimensions:
              bins = 200 ;
              ens_sizes = UNLIMITED ; // (40 currently)
      variables:
              int count(ens_sizes, bins) ;
                      count:description = "number of samples in each bin" ;
              double true_corr_mean(ens_sizes, bins) ;
              double alpha(ens_sizes, bins) ;
                      alpha:description = "sampling error correction factors" ;
              int ens_sizes(ens_sizes) ;
                      ens_sizes:description = "ensemble size used for calculation" ;

      // global attributes:
                      :num_samples = 100000000 ;
                      :title = "Sampling Error Corrections for fixed ensemble sizes." ;
                      :reference = "Anderson, J., 2012: Localization and Sampling Error 
                                    Correction in Ensemble Kalman Filter Data Assimilation.
                                    Mon. Wea. Rev., 140, 2359-2371, doi: 10.1175/MWR-D-11-00013.1." ;
                      :version = "" ;
      data:

      These ensemble sizes are already supported!
       ens_sizes = 5,  6,  7,  8,  9, 10, 12, 14, 15, 16, 18, 20, 22, 24, 28, 30, 32, 36, 40, 44,
                  48, 49, 50, 52, 56, 60, 64, 70, 72, 80, 84, 88, 90, 96, 100, 120, 140, 160, 180, 200
      }

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &gen_sampling_error_table_nml
      ens_sizes = 5,  6,  7,  8,  9, 10, 12, 14, 15, 16, 18, 20, 22, 24, 28, 30, 32, 36, 40, 44,
                 48, 49, 50, 52, 56, 60, 64, 70, 72, 80, 84, 88, 90, 96, 100, 120, 140, 160, 180, 200
      debug = .false.
      /

| 

+-----------+--------------+-----------------------------------------------------------------------------------------+
| Item      | Type         | Description                                                                             |
+===========+==============+=========================================================================================+
| ens_sizes | integer(200) | An array of ensemble sizes to compute. Any new size gets appended to the variables in   |
|           |              | the netCDF file. Any order is fine, the array does not have to be monotonic. The        |
|           |              | numbers listed in the example exist in the file distributed with DART. *Do not get      |
|           |              | carried away by generating a lot of new ensemble sizes in one execution.* The table of  |
|           |              | run-time above should give you some indication of how long it takes to create a new     |
|           |              | entry.                                                                                  |
+-----------+--------------+-----------------------------------------------------------------------------------------+
| debug     | logical      | A switch to add some run-time output. Generally not needed.                             |
+-----------+--------------+-----------------------------------------------------------------------------------------+

Modules used
------------

::

   types_mod
   utilities_mod
   random_seq_mod

-  ``input.nml`` for the run-time input
-  ``sampling_error_correction_table.nc`` is both read and written. Any new ensemble sizes are simply appended to the
   file.
-  ``dart_log.out`` has the run-time output.

.. _section-1:

-  ``input.nml`` for the run-time input
-  ``final_full.N`` are created - N is the ensemble size.
-  ``dart_log.out`` has the run-time output.

References
----------

-  **Anderson, J. L.**, 2012: Localization and Sampling Error Correction in Ensemble Kalman Filter Data Assimilation.
   *Mon. Wea. Rev.*, **140**, 2359-2371 `doi:
   10.1175/MWR-D-11-00013.1 <http://dx.doi.org/doi:10.1175/MWR-D-11-00013.1>`__
