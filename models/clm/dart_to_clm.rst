PROGRAM ``dart_to_clm``
=======================

Overview
--------

``dart_to_clm`` replaces the contents of a CLM restart file with the posterior values 
from DART. Only variables with posterior values are updated. The *_FillValue* in the 
DART posterior is used as a mask when replacing the value. The original CLM restart file values 
are maintained when the DART posterior value is *_FillValue*.

Usage
-----
``dart_to_clm`` overwrites the output file. See the companion :doc:`clm_to_dart` for a 
discussion on replacing indeterminate values in the CLM restart snow variables with 
*_FillValue* so that the DART netCDF read routines handle the special values correctly. 
Where the DART posterior file has *_FillValue* values, the original CLM restart file 
is left unchanged, preserving the original values.

``dart_to_clm`` also **includes an option that repartitions the snow water equivalent (SWE)**
within the CLM snow layer variables (``repartition_swe``).  This algorithm functions as
follows: if during the DART *filter* step an  increment is applied to the CLM SWE variable 
(``H2OSNO``), this algorithm partitions the SWE increment amongst the active snow
layers for each CLM column to maintain the same relative ice and liquid distribution 
amongst the CLM ice and liquid variables (``H2OSOI_ICE`` and ``H2OSOI_LIQ``).  Then, 
based upon the change in snow layer ice mass, the algorithm updates the snow layer
dimensions (``SNOW_DEPTH``, ``DZSNO``, ``ZSNO`` and ``ZISNO``) assuming the prior snow
density remains the same.  

This repartitioning approach has two main advantages over the traditional *filter* update
applied by DART:

  1. It is not required that a given snow layer exists amongst all CLM ensemble members
  in order for the snow layer to be adjusted. Because the CLM snow algorithm splits and
  aggregates snow layers based upon the vertical thickness in each layer, the uncertainty
  in the CAM reanalysis (ie. snowfall, melting of snow)  may lead to CLM members with 
  different number of snow layers. The default version of DART must ignore these snow layers
  and no update is applied. When the ``repartition_swe`` algorithm is turned on a snow layer
  which is missing in a subset of ensemble members can still be updated.  The only requirement
  is that each ensemble member includes at least one active snow layer for a CLM column.

  2. The repartitioning algorithm guarantees the snow layer adjustments (both mass and dimensions)
  are consistent with the update to the SWE clm variable (``H2OSNO``). In other words,
  the sum of the ice/liquid mass for all the vertical snow layers is identical to the increment
  of the total column ice/liquid mass.  In the default implementation of DART there is no
  consistency between the total column update of snow and the update of the snow layer variables,
  because the snow layer adjustment is based upon the regression between the expected
  observation and each snow layer.      

Whether to repartition the snow characteristics (mass, dimensions) across all active snow layers 
(repartition_swe = 1) or to the bottom snow layer only (repartition_swe = 2) depends upon
the application.  Repartitioning to all snow layers maintains the prior relative
mass distribution of ice and liquid which could be important for certain water-limited 
regimes with influence upon the carbon cycle.  On the other hand, repartitioning to all snow
layers indirectly influences the albedo characteristics of the snowcover as snow variables
such as black carbon and dust are concentrated/diluted depending upon the snow layer mass
adjustment. Repartitioning the snow/ice to only the bottom snow layer avoids 
indirect changes to the snow albedo characteristics because the black carbon and dust of
the topmost snow layer (exposed to sun) is conserved. This could be important for applications
where snow albedo, suface energy balance, and surface temperature are critical. However, 
applying the mass adjustment to the bottom snow layer only, does not maintain the same 
relative mass distribution and given the proximity to the soil, may provide a more immediate
impact on soil moisture conditions.    


Namelist
--------

Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent
them from prematurely terminating the namelist. These are the defaults:

::

   &dart_to_clm_nml
      dart_to_clm_input_file    = 'dart_posterior.nc'
      dart_to_clm_output_file   = 'clm_restart.nc'
      repartition_swe           = 0
      repartition_vhist_file    = 'clm_vector_history.nc'
      repartition_analysis_file = 'dart_posterior_vector.nc'
      verbose                   = 0
      /




.. container::


   ========================= =================== =============================================================== 
   Item                      Type                Description                                                     
   ========================= =================== =============================================================== 
   dart_to_clm_input_file    character(len=256)  Name of the DART posterior (the output of the filter program)
   dart_to_clm_output_file   character(len=256)  Name of the CLM restart file to modify. 
   repartition_swe           integer             | Flag to choose snow repartitioning algorithm  
                                                 | 0   Snow repartitioning off. Uses default DART update
                                                 | 1   Snow repartitioned to all active snow layers
                                                 | 2   Snow repartitioned to bottom snow layer only
   repartition_vhist_file    character(len=256)  Name of the CLM vector history file (prior SWE value)
                                                 Only required when repartition_swe = 1,2
   repartition_analysis_file character(len=256)  Name of the DART analysis file (posterior SWE value)
                                                 Only required when repartition_swe = 1,2
   verbose                   integer             | Flag to control how much run-time output is created.
                                                 | 0   is very little output.
                                                 | 1   reports which CLM variables are being updated.
                                                 | 3   if repartition_swe = 1,2 reports on snow layer diagnostic 
                                                 |     variables during repartitioning process
   ========================= =================== ===============================================================



Variable requirements when snow repartitioning is enabled
---------------------------------------------------------

When swe_repartition = 1,2 the DART code requires CLM mass and dimensional snow layer variables output in
a specific format described in the table below. A list of the  CLM variables (with dimensions) are listed
in the left most column of the table.  The CLM dimensions important for snow variables are as follows. A
'column' is the land unit at which snow properties are defined. The 'levtot' is a vertical dimension
of the total number of snow layers (1-12) and soil layers (13-37) within a column. The 'levsno' dimension
is the total number of snow layers.   


.. container::

   ========================== =================== ================================================================= 
   Variable (dimension)       File Type           Description                                                     
   ========================== =================== ================================================================= 
   H2OSNO (column)            clm vector history  The prior snow water equivalent. Needs to be output in vector 
                                                  format within the CLM history file (h2). Must be included in 
                                                  DART state ``clm_variables`` (``&model_nml``)
   H2OSNO (column)            DART analysis       The posterior snow water equivalent output from the *filter*
                              stage               analysis stage.  This requires that 'analysis' stage be output
                                                  within the ``stages_to_write`` (``&filter_nml``)
   H2OSOI_ICE (column,levtot) clm restart         The ice mass in each snow/soil layer.  Must be included in
                                                  DART state ``clm_variables`` (``&model_nml``)
   H2OSOI_LIQ (column,levtot) clm restart         The liquid mass in each snow/soil layer.  Must be included in
                                                  DART state ``clm_variables`` (``&model_nml``)
   DZSNO (column,levsno)      clm restart         The snow layer thickness. Must be included in
                                                  DART state ``clm_variables`` (``&model_nml``)
   ZISNO (column,levsno)      clm restart         The top interface depth of each snow layer. Must be included in
                                                  DART state ``clm_variables`` (``&model_nml``)
   ZSNO (column,levsno)       clm restart         The middle depth of each snow layer. Must be in included in
                                                  DART state ``clm_variables`` (``&model_nml``)
   SNOW_DEPTH (column)        clm restart         The total snow depth (sum of all layers). Must be included in
                                                  DART state ``clm_variables`` (``&model_nml``)       
   ========================== =================== =================================================================


When adjusting snow layer variables within an assimilation, at a minimum, the ``clm_variables``
within ``model_nml`` must include the following: 

::

  clm_variables  = 'H2OSNO',      'QTY_SNOW_WATER',             '0.0', 'NA', 'vector'  , 'NO_COPY_BACK', 
                 'SNOW_DEPTH',  'QTY_SNOW_THICKNESS',         '0.0', 'NA', 'restart' , 'UPDATE',
                 'H2OSOI_LIQ',  'QTY_SOIL_LIQUID_WATER',      '0.0', 'NA', 'restart' , 'UPDATE',
                 'H2OSOI_ICE',  'QTY_SOIL_ICE',               '0.0', 'NA', 'restart' , 'UPDATE',
                 'DZSNO',       'QTY_SNOW_THICKNESS',         '0.0', 'NA', 'restart' , 'UPDATE',
                 'ZSNO',        'QTY_SNOW_THICKNESS',         'NA',  'NA', 'restart' , 'UPDATE',
                 'ZISNO',       'QTY_SNOW_THICKNESS',         'NA',  'NA', 'restart' , 'UPDATE',
                 /


.. note::

     The H2OSOI_ICE and H2OSOI_LIQ variables include both snow layer and subsurface 
     soil layers. **Only the snow and liquid mass within the snow layers are repartitioned,**
     whereas the subsurface layers are updated through the default DART approach
     based upon the regression relationship between the subsurface layer property and an
     expected observation.    



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

- ``clm_restart.nc`` is the CLM generated netCDF file that is modified.

- ``clm_vector_history.nc`` is the CLM generated netCDF file that provides the prior estimate of snow water equivalent (SWE)
- ``dart_posterior_vector.nc`` is the DART generated analysis stage file that provide posterior estimate of SWE

- ``dart_log.out`` list directed output from the ``dart_to_clm``.


References
----------

none
