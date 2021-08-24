PROGRAM ``dart_to_clm``
=======================

Overview
--------

``dart_to_clm`` replaces the contents of a CLM restart file with the posterior values 
from DART. Only variables with posterior values are updated. The *_FillValue* in the 
DART posterior is used as a mask when replacing the value. The original CLM restart file values 
are maintained when the DART posterior value is *_FillValue*.

Another intended use for ``dart_to_clm`` is to use the posterior snow water equivalent 
(‘SWE’ - the CLM variable *H2OSNO*) to update the prognostic snow variables in CLM. 
**This is not implemented yet.** The method will be to simply use the existing ratios 
of snow in the layers and change them proportionally to achieve the desired 
posterior SWE.

Usage
-----
``dart_to_clm`` overwrites the output file. See the companion :doc:`clm_to_dart` for a 
discussion on replacing indeterminate values in the CLM restart snow variables with 
*_FillValue* so that the DART netCDF read routines handle the special values correctly. 
Where the DART posterior file has *_FillValue* values, the original CLM restart file 
is left unchanged, preserving the original values.

Namelist
--------

Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent
them from prematurely terminating the namelist. These are the defaults:

::

   &dart_to_clm_nml
      dart_to_clm_input_file  = 'dart_posterior.nc'
      dart_to_clm_output_file = 'clm_restart.nc'
      repartition_swe         = .false.
      verbose                 = 0
      /


.. container::


   ======================= =================== ================================================================= 
   Item                    Type                Description                                                     
   ======================= =================== ================================================================= 
   dart_to_clm_input_file  character(len=256)  Name of the DART posterior (the output of the filter program)
   dart_to_clm_output_file character(len=256)  Name of the CLM restart file to modify. 
                                               | This will be used for the next CLM model advance.
   repartition_swe         logical             Update the prognostic snow variables in CLM such that 
                                               the posterior snow water equivalent is correct.
                                               THIS IS NOT IMPLEMENTED YET.
   verbose                 integer             | Flag to control how much run-time output is created.
                                               | 0   is very little output.
                                               | 1   reports which variables are being updated.
   ======================= =================== ================================================================= 


Modules used
------------

::

   assimilation_code/location/threed_sphere/location_mod.f90
   assimilation_code/location/utilities/location_io_mod.f90
   assimilation_code/modules/assimilation/adaptive_inflate_mod.f90
   assimilation_code/modules/assimilation/assim_model_mod.f90
   assimilation_code/modules/io/dart_time_io_mod.f90
   assimilation_code/modules/io/direct_netcdf_mod.f90
   assimilation_code/modules/io/io_filenames_mod.f90
   assimilation_code/modules/io/state_structure_mod.f90
   assimilation_code/modules/io/state_vector_io_mod.f90
   assimilation_code/modules/observations/obs_kind_mod.f90
   assimilation_code/modules/observations/obs_sequence_mod.f90
   assimilation_code/modules/utilities/distributed_state_mod.f90
   assimilation_code/modules/utilities/ensemble_manager_mod.f90
   assimilation_code/modules/utilities/netcdf_utilities_mod.f90
   assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
   assimilation_code/modules/utilities/null_win_mod.f90
   assimilation_code/modules/utilities/options_mod.f90
   assimilation_code/modules/utilities/random_seq_mod.f90
   assimilation_code/modules/utilities/sort_mod.f90
   assimilation_code/modules/utilities/time_manager_mod.f90
   assimilation_code/modules/utilities/types_mod.f90
   assimilation_code/modules/utilities/utilities_mod.f90
   models/clm/model_mod.f90
   models/utilities/default_model_mod.f90
   observations/forward_operators/obs_def_mod.f90
   observations/forward_operators/obs_def_utilities_mod.f90


Files
-----

- ``input.nml`` is used for ``dart_to_clm``

- ``dart_posterior.nc`` is one of the netCDF files output from the *filter* program.

- ``clm_restart.nc`` is the netCDF file that is modified.

- ``dart_log.out`` list directed output from the ``dart_to_clm``.


References
----------

none
