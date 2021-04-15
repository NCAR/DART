ROMS observations to DART observation sequences
===============================================

Overview
--------

The relationship between ROMS and DART is slightly different than most other models. ROMS has the ability to apply its
own forward operator as the model is advancing (a capability needed for variational assimilation) which produces
something the ROMS community calls '*verification*' observations. The observation file that is input to ROMS is
specified by the ``s4dvar.in``:``OBSname`` variable. The verification obs are written out to a netcdf file whose name is
specified by the ``s4dvar.in``:``MODname`` variable. Since each ROMS model is advancing independently, a set of
verification observation files are created during a DART/ROMS assimilation cycle. This set of files can be converted
using ``convert_roms_obs`` to produce a DART observation sequence file that has precomputed forward operators (FOs).
``convert_roms_obs`` can also convert ``s4dvar.in``:``OBSname``,\ ``MODname`` files to a DART observation sequence file
that does not have the precomputed FOs.

The ROMS verification observation files **must** contain the **obs_provenance as a global attribute** and the following
variables:

-  *obs_lat, obs_lon, obs_depth*
-  *obs_value*
-  *obs_error*
-  *obs_time*
-  *NLmodel_value*
-  *obs_scale*
-  *obs_provenance*

Note that the *obs_provenance:flag_values*, and *obs_provenance:flag_meanings* attributes are totally ignored - those
relationships are specified by the global attribute **obs_provenance**.

Locations only specified by *obs_Xgrid, obs_Ygrid, obs_depth* are **not** supported.

| The conversion of a (set of) ROMS verification observations requires metadata to coordinate the relationship of the
  ROMS observation provenance to a DART observation TYPE. ROMS provides significant flexibility when specifying the
  observation provenance and it is simply impractical for DART to try to support all of them. An example of the current
  practice is described in the PROGRAMS section below.
| **Important:** ``filter`` and ``perfect_model_obs`` must also be informed which DART observation types use precomputed
  forward operators. This is done by setting the ``input.nml``\ ``&obs_kind_nml`` namelist. An example is shown at the
  end of the PROGRAMS section below.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &convert_roms_obs_nml
      ens_size               = 1
      roms_mod_obs_files     = ''
      roms_mod_obs_filelist  = 'filelist.txt'
      dart_output_obs_file   = 'obs_seq.out'
      append_to_existing     = .false.
      use_precomputed_values = .true.
      add_random_noise       = .false.
      pert_amplitude         = 0.01
      verbose                = 0
      type_translations      = 'NULL'
     /

| 

.. container::

   +------------------------+------------------------------------+------------------------------------------------------+
   | Item                   | Type                               | Description                                          |
   +========================+====================================+======================================================+
   | ens_size               | integer                            | Number of ensemble members which are expected to be  |
   |                        |                                    | found when creating the expected obs values. This    |
   |                        |                                    | must match the number of ROMS "mod" files listed in  |
   |                        |                                    | either the 'roms_mod_obs_files' or                   |
   |                        |                                    | 'roms_mod_obs_filelist' namelist items. It is an     |
   |                        |                                    | error if they are not the same length.               |
   +------------------------+------------------------------------+------------------------------------------------------+
   | roms_mods_obs_files    | character(len=256), dimension(100) | List of filenames, one per ensemble member, that     |
   |                        |                                    | contain the observation values for each ensemble     |
   |                        |                                    | member. These are output from the ROMS program. If   |
   |                        |                                    | listing the files explicitly in this list,           |
   |                        |                                    | 'roms_mod_obs_filelist' must be ' ' (null).          |
   +------------------------+------------------------------------+------------------------------------------------------+
   | roms_mods_obs_filelist | character(len=256)                 | The name of an ASCII file which contains, one per    |
   |                        |                                    | line, a list of filenames, one per ensemble member,  |
   |                        |                                    | that contain the expected obs values for each        |
   |                        |                                    | ensemble member. The filenames should NOT be quoted. |
   |                        |                                    | These are output from the ROMS program. If using a   |
   |                        |                                    | filelist, then 'roms_mod_obs_files' must be ' '      |
   |                        |                                    | (null).                                              |
   +------------------------+------------------------------------+------------------------------------------------------+
   | dart_output_obs_file   | character(len=256)                 | The name of the DART obs_seq file to create. If a    |
   |                        |                                    | file already exists with this name, it is either     |
   |                        |                                    | appended to or overwritten depending on the          |
   |                        |                                    | 'append_to_existing' setting below.                  |
   +------------------------+------------------------------------+------------------------------------------------------+
   | append_to_existing     | logical                            | If an existing 'dart_output_obs_file' is found, this |
   |                        |                                    | namelist item controls how it is handled. If .true.  |
   |                        |                                    | the new observations are appended to the existing    |
   |                        |                                    | file. If .false. the new observations overwrite the  |
   |                        |                                    | existing file.                                       |
   +------------------------+------------------------------------+------------------------------------------------------+
   | use_precomputed_values | logical                            | flag to indicate that the output DART observation    |
   |                        |                                    | sequence file should include the verification        |
   |                        |                                    | observation values from all of the ROMS observation  |
   |                        |                                    | files. If ``.true.`` this will result in the DART    |
   |                        |                                    | file having the precomputed FOs to be used in the    |
   |                        |                                    | DART assimilation. If ``.false.`` this will result   |
   |                        |                                    | in DART files having the instrument values only.     |
   +------------------------+------------------------------------+------------------------------------------------------+
   | add_random_noise       | logical                            | Almost always should be .false. . The exception is   |
   |                        |                                    | the first cycle of an assimilation if all the ROMS   |
   |                        |                                    | input files are identical (no ensemble currently     |
   |                        |                                    | exists). To create differences in the forward        |
   |                        |                                    | operator values (since they are computed by ROMS),   |
   |                        |                                    | we can add gaussian noise here to give them          |
   |                        |                                    | perturbed values. This should be set as well as the  |
   |                        |                                    | "perturb_from_single_instance = .true." namelist in  |
   |                        |                                    | the ``&filter_nml`` namelist. After the first cycle, |
   |                        |                                    | both these should be set back to .false. .           |
   +------------------------+------------------------------------+------------------------------------------------------+
   | pert_amplitude         | real(r8)                           | Ignored unless 'add_random_noise' is .true. .        |
   |                        |                                    | Controls the range of random values added to the     |
   |                        |                                    | expected obs values. Sets the width of a gaussian.   |
   +------------------------+------------------------------------+------------------------------------------------------+
   | verbose                | integer                            | If greater than 0, prints more information during    |
   |                        |                                    | the conversion.                                      |
   +------------------------+------------------------------------+------------------------------------------------------+
   | type_translations      | character(256), dimension(2, 100)  | A set of strings which control the mapping of ROMS   |
   |                        |                                    | observation types to DART observation types. These   |
   |                        |                                    | should be specified in pairs. The first column       |
   |                        |                                    | should be a string that occurs in the global         |
   |                        |                                    | attribute '``obs_provenance``'. Note that the        |
   |                        |                                    | ``obs_provenance:flag_values`` and                   |
   |                        |                                    | ``obs_provenance:flag_meanings`` attributes are      |
   |                        |                                    | ignored. The second column should be a DART specific |
   |                        |                                    | obs type that is found in                            |
   |                        |                                    | ``DART/assimi                                        |
   |                        |                                    | lation_code/modules/observations/obs_kind_mod.f90``, |
   |                        |                                    | which is created by the DART ``preprocess`` program. |
   +------------------------+------------------------------------+------------------------------------------------------+

| 

Data sources
------------

The origin of the input observation files used by ROMS are completely unknown to me.

Programs
--------

-  ``convert_roms_obs``
-  :doc:`../../../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf`
-  :doc:`../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool`
-  :doc:`../../../assimilation_code/programs/preprocess/preprocess`
-  :doc:`../../../assimilation_code/programs/advance_time/advance_time`

Only ``convert_roms_obs`` will be discussed here.

The **global attribute** ``obs_provenance`` is used to relate the observation provenance to DART observation TYPES. The
ROMS 'MODname' netCDF file(s) must have both the ``obs_provenance`` variable and a ``obs_provenance`` **global
attribute**. The **exact** strings must be repeated in the DART ``convert_roms_obs_nml``:``type_translations`` variable
to be able to convert from the integer value of the obs_provenance to th DART type in the following example:

``ncdump -h roms_mod_obs.nc`` (the output has been pruned for clarity)

::

   netcdf roms_mod_obs {
   dimensions:
           record = 2 ;
           survey = 5376 ;
           state_var = 8 ;
           datum = 2407217 ;
   variables:
           {snip}
           int obs_provenance(datum) ;
                   obs_provenance:long_name = "observation origin" ;
                   obs_provenance:flag_values = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ;
           double obs_time(datum) ;
                   obs_time:long_name = "time of observation" ;
                   obs_time:units = "days since 1900-01-01 00:00:00 GMT" ;
                   obs_time:calendar = "gregorian" ;
           double obs_lon(datum) ;
                   obs_lon:long_name = "observation longitude" ;
                   obs_lon:units = "degrees_east" ;
           double obs_lat(datum) ;
                   obs_lat:long_name = "observation latitude" ;
                   obs_lat:units = "degrees_north" ;
           double obs_depth(datum) ;
                   obs_depth:long_name = "ROMS internal depth of observation variable" ;
                   obs_depth:units = "meters or fractional z-levels" ;
                   obs_depth:negative_value = "downwards" ;
                   obs_depth:missing_value = 1.e+37 ;
           double obs_error(datum) ;
                   obs_error:long_name = "observation error covariance" ;
           double obs_value(datum) ;
                   obs_value:long_name = "observation value" ;
           double obs_scale(datum) ;
                   obs_scale:long_name = "observation screening/normalization scale" ;
                   obs_scale:_FillValue = 0. ;
           double NLmodel_value(datum) ;
                   NLmodel_value:long_name = "nonlinear model at observation locations" ;
                   NLmodel_value:_FillValue = 1.e+37 ;
           {snip}
        :obs_provenance = "\n",
                "1: gridded AVISO sea level anomaly (zeta)\n",
                "2: gridded Aquarius SSS (salinity)\n",
                "3: XBT from Met Office (temperature)\n",
                "4: CTD from Met Office (temperature)\n",
                "5: CTD from Met Office (salinity)\n",
                "6: ARGO floats (temperature)\n",
                "7: ARGO floats (salinity)\n",
                "8: glider UCSD (temperature)\n",
                "9: glider UCSD (salinity)\n",
                "10: blended satellite SST (temperature)" ;
           {snip}

| Note the integer values that start the obs_provenance strings are used to interpret the integer contents of the
  obs_provenance variable. They need not be consecutive, nor in any particular order, but they must not appear more than
  once.
| The following is the relevent section of the DART ``input.nml``:

::

   &convert_roms_obs_nml
      ens_size               = 32
      roms_mod_obs_filelist  = 'precomputed_files.txt'
      dart_output_obs_file   = 'obs_seq.out'
      append_to_existing     = .false.
      use_precomputed_values = .true.
      add_random_noise       = .false.
      verbose                = 1
      type_translations = "gridded AVISO sea level anomaly (zeta)", "SATELLITE_SSH",
                          "gridded Aquarius SSS (salinity)",        "SATELLITE_SSS",
                          "XBT from Met Office (temperature)",      "XBT_TEMPERATURE",
                          "CTD from Met Office (temperature)",      "CTD_TEMPERATURE",
                          "CTD from Met Office (salinity)",         "CTD_SALINITY",
                          "ARGO floats (temperature)",              "ARGO_TEMPERATURE",
                          "ARGO floats (salinity)",                 "ARGO_SALINITY",
                          "glider UCSD (temperature)",              "GLIDER_TEMPERATURE",
                          "glider UCSD (salinity)",                 "GLIDER_SALINITY",
                          "blended satellite SST (temperature)",    "SATELLITE_BLENDED_SST"
     /

A complete list of DART observation TYPES for oceans is described in
:doc:`../../forward_operators/obs_def_ocean_mod`

Any or all of the DART observation types that appear in the second column of ``type_translations`` must also be
designated as observations that have precomputed forward operators. This is done by setting the
``input.nml``\ ``&obs_kind_nml`` namelist as follows:

::

   &obs_kind_nml
      assimilate_these_obs_types =          'SATELLITE_SSH',
                                            'SATELLITE_SSS',
                                            'XBT_TEMPERATURE',
                                            'CTD_TEMPERATURE',
                                            'CTD_SALINITY',
                                            'ARGO_TEMPERATURE',
                                            'ARGO_SALINITY',
                                            'GLIDER_TEMPERATURE',
                                            'GLIDER_SALINITY',
                                            'SATELLITE_BLENDED_SST'
      use_precomputed_FOs_these_obs_types = 'SATELLITE_SSH',
                                            'SATELLITE_SSS',
                                            'XBT_TEMPERATURE',
                                            'CTD_TEMPERATURE',
                                            'CTD_SALINITY',
                                            'ARGO_TEMPERATURE',
                                            'ARGO_SALINITY',
                                            'GLIDER_TEMPERATURE',
                                            'GLIDER_SALINITY',
                                            'SATELLITE_BLENDED_SST'
     /
