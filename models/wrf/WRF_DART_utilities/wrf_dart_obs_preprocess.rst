PROGRAM ``wrf_dart_obs_preprocess``
===================================

Overview
--------

Program to preprocess observations, with specific knowledge of the WRF domain.

This program will exclude all observations outside of the given WRF domain. There are options to exclude or increase the
error values of obs close to the domain boundaries. The program can superob (average) aircraft and satellite wind obs if
they are too dense.

This program can read up to 9 additional obs_seq files and merge their data in with the basic obs_sequence file which is
the main input.

This program can reject surface observations if the elevation encoded in the observation is too different from the wrf
surface elevation.

This program can exclude observations above a specified height or pressure.

This program can overwrite the incoming Data QC value with another.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &wrf_obs_preproc_nml

     file_name_input          = 'obs_seq.old'
     file_name_output         = 'obs_seq.new'
     
     sonde_extra              = 'obs_seq.rawin'
     land_sfc_extra           = 'obs_seq.land_sfc'
     metar_extra              = 'obs_seq.metar'
     marine_sfc_extra         = 'obs_seq.marine'
     sat_wind_extra           = 'obs_seq.satwnd'
     profiler_extra           = 'obs_seq.profiler'
     gpsro_extra              = 'obs_seq.gpsro'
     acars_extra              = 'obs_seq.acars'
     trop_cyclone_extra       = 'obs_seq.tc'
     
     overwrite_obs_time       = .false.  
     
     obs_boundary             = 0.0
     increase_bdy_error       = .false.  
     maxobsfac                = 2.5   
     obsdistbdy               = 15.0  
     
     sfc_elevation_check      = .false.  
     sfc_elevation_tol        = 300.0  
     obs_pressure_top         = 0.0  
     obs_height_top           = 2.0e10  
     
     include_sig_data         = .true.   
     tc_sonde_radii           = -1.0  
     
     superob_aircraft         = .false.  
     aircraft_horiz_int       = 36.0  
     aircraft_pres_int        = 2500.0  
     
     superob_sat_winds        = .false.    
     sat_wind_horiz_int       = 100.0   
     sat_wind_pres_int        = 2500.0  
     
     overwrite_ncep_satwnd_qc = .false.    
     overwrite_ncep_sfc_qc    = .false.  
   /

+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Item**                                                      | **Type**           | **Description**                                                                                                                                                                                                                                          |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Generic parameters:                                           |                    |                                                                                                                                                                                                                                                          |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| file_name_input                                               | character(len=129) | The input obs_seq file.                                                                                                                                                                                                                                  |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| file_name_output                                              | character(len=129) | The output obs_seq file.                                                                                                                                                                                                                                 |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| | sonde_extra, land_sfc_extra,                                | character(len=129) | The names of additional input obs_seq files, which if they exist, will be merged in with the obs from the file_name_input obs_seq file. If the files do not exist, they are silently ignored without error.                                              |
| | metar_extra,                                                |                    |                                                                                                                                                                                                                                                          |
| | marine_sfc_extra, sat_wind_extra,                           |                    |                                                                                                                                                                                                                                                          |
| | profiler_extra, gpsro_extra,                                |                    |                                                                                                                                                                                                                                                          |
| | acars_extra, trop_cyclone_extra                             |                    |                                                                                                                                                                                                                                                          |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| overwrite_obs_time                                            | logical            | If true, replace the incoming observation time with the analysis time. Not recommended.                                                                                                                                                                  |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Boundary-specific parameters:**                             |                    |                                                                                                                                                                                                                                                          |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| obs_boundary                                                  | real(r8)           | Number of grid points around domain boundary which will be considered the new extent of the domain. Observations outside this smaller area will be excluded.                                                                                             |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| increase_bdy_error                                            | logical            | If true, observations near the domain boundary will have their observation error increased by maxobsfac.                                                                                                                                                 |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| maxobsfac                                                     | real(r8)           | If increase_bdy_error is true, multiply the error by a ramped factor. This item sets the maximum error.                                                                                                                                                  |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| obsdistbdy                                                    | real(r8)           | If increase_bdy_error is true, this defines the region around the boundary (in number of grid points) where the observation error values will be altered. This is ramped, so when you reach the innermost points the change in observation error is 0.0. |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Parameters to reduce observation count :**                  |                    |                                                                                                                                                                                                                                                          |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| sfc_elevation_check                                           | logical            | If true, check the height of surface observations against the surface height in the model.                                                                                                                                                               |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| sfc_elevation_tol                                             | real(r8)           | If sfc_elevation_check is true, the maximum difference between the elevation of a surface observation and the model surface height, in meters. If the difference is larger than this value, the observation is excluded.                                 |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| obs_pressure_top                                              | real(r8)           | Observations with a vertical coordinate in pressure which are located above this pressure level (i.e. the obs vertical value is smaller than the given pressure) will be excluded.                                                                       |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| obs_height_top                                                | real(r8)           | Observations with a vertical coordinate in height which are located above this height value (i.e. the obs vertical value is larger than the given height) will be excluded.                                                                              |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Radio/Rawinsonde-specific parameters :**                    |                    |                                                                                                                                                                                                                                                          |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| include_sig_data                                              | logical            | If true, include significant level data from radiosondes.                                                                                                                                                                                                |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| tc_sonde_radii                                                | real(r8)           | If greater than 0.0 remove any sonde observations closer than this distance in Kilometers to the center of a Tropical Cyclone.                                                                                                                           |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Aircraft-specific parameters :                                |                    |                                                                                                                                                                                                                                                          |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| superob_aircraft                                              | logical            | If true, average all aircraft observations within the given radius and output only a single observation. Any observation that is used in computing a superob observation is removed from the list and is not used in any other superob computation.      |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| aircraft_horiz_int                                            | real(r8)           | If superob_aircraft is true, the horizontal distance in Kilometers which defines the superob area. All other unused aircraft observations within this radius will be averaged with the current observation.                                              |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| aircraft_vert_int                                             | real(r8)           | If superob_aircraft is true, the vertical distance in Pascals which defines the maximum separation for including an observation in the superob computation.                                                                                              |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Satellite Wind-specific parameters :**                      |                    |                                                                                                                                                                                                                                                          |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| superob_sat_winds                                             | logical            | If true, average all sat_wind observations within the given radius and output only a single observation. Any observation that is used in computing a superob observation is removed from the list and is not used in any other superob computation.      |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| sat_wind_horiz_int                                            | real(r8)           | If superob_sat_winds is true, the horizontal distance in Kilometers which defines the superob area. All other unused sat_wind observations within this radius will be averaged with the current observation.                                             |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| sat_wind_vert_int                                             | real(r8)           | If superob_sat_winds is true, the vertical distance in Pascals which defines the maximum separation for including an observation in the superob computation.                                                                                             |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| overwrite_ncep_satwnd_qc                                      | logical            | If true, replace the incoming Data QC value in satellite wind observations with 2.0.                                                                                                                                                                     |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Surface Observation-specific parameters :**                 |                    |                                                                                                                                                                                                                                                          |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| overwrite_ncep_sfc_qc                                         | logical            | If true, replace the incoming Data QC value in surface observations with 2.0.                                                                                                                                                                            |
+---------------------------------------------------------------+--------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+ 

 

Modules used
------------

::

   types_mod
   obs_sequence_mod
   utilities_mod
   obs_kind_mod
   time_manager_mod
   model_mod
   netcdf

Files
-----

-  Input namelist ; ``input.nml``
-  Input WRF state netCDF files; ``wrfinput_d01, wrfinput_d02, ...``
-  Input obs_seq files (as specified in namelist)
-  Output obs_seq file (as specified in namelist)

File formats
~~~~~~~~~~~~

This utility can read one or more obs_seq files and combine them while doing the rest of the processing. It uses the
standard DART observation sequence file format.

References
----------

-  Generously contributed by Ryan Torn.
