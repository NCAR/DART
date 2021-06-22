PROGRAM preprocess
==================

Overview
--------

Preprocess is a DART-supplied preprocessor program. Preprocess is used to insert observation specific code into DART at
compile time.

In DART, forward operators are not specific to any one model. To achieve this separation between models and forward
operators DART makes a distinction between an observation *type* and a physical *quantity*. For example, a radiosonde
used to measure windspeed would be a *type* of observation. Zonal wind and meridional wind are *quantities* used to
calculate windspeed. Specifying many observation types allows DART to be able to evaluate some observations and
assimilate others even if the instruments measure the same quantity.

Preprocess takes user supplied observation and quantity files and combines them with template files to produce code for
DART. Use the namelist option 'obs_type_files' to specify the input observation files and the namelist option
'quantity_files' to specify the input quantity files.

-  If no quantity files are given, a default list of quantities is used.
-  If no obs_type_files are given, only identity observations can be used in the filter (i.e. the state variable values
   are directly observed; forward operator is an identity)

The template files ``DEFAULT_obs_def_mod.F90`` and ``DEFAULT_obs_kind_mod.F90`` contain specially formatted comment
lines. These comment lines are used as markers to insert observation specific information. Prepreocess relies these
comment lines being used *verbatim*.

There is no need to to alter ``DEFAULT_obs_def_mod.F90`` or ``DEFAULT_obs_kind_mod.F90``. Detailed instructions for
adding new observation types can be found in :doc:`../../../observations/forward_operators/obs_def_mod`. New quantities
should be added to a quantity file, for example a new atmosphere quantity should be added to
``atmosphere_quantities_mod.f90``.

Every line in a quantity file between the start and end markers must be a comment or a quantity definition (QTY_string).
Multiple name-value pairs can be specified for a quantity but are not required. For example, temperature may be defined:
``! QTY_TEMPERATURE units="K" minval=0.0``. Comments are allowed between quantity definitions or on the same line as the
definition. The code snippet below shows acceptable formats for quantity definitions

.. container :: routine

   ::

      ! BEGIN DART PREPROCESS QUANTITY DEFINITIONS
      !
      ! Formats accepted:
      !
      ! QTY_string 
      ! QTY_string name=value 
      ! QTY_string name=value name2=value2 
      ! 
      ! QTY_string ! comments 
      ! 
      ! ! comment
      ! 
      ! END DART PREPROCESS QUANTITY DEFINITIONS  


| The output files produced by preprocess are named ``assimilation_code/modules/observations/obs_kind_mod.f90`` and
  ``observations/forward_operators/obs_def_mod.f90``, but can be renamed by namelist control if needed. Be aware that if
  you change the name of these output files, you will need to change the path_names files for DART executables.

Namelist
--------

When you run preprocess, the namelist is read from the file ``input.nml`` in the directory where preprocess is run.

Namelists start with an ampersand '&' and terminate with a slash '/'. Character strings that contain a '/' must be
enclosed in quotes to prevent them from prematurely terminating the namelist. These are the defaults:

::

   &preprocess_nml
     overwrite_output        = .true.,
     input_obs_def_mod_file  = '../../../observations/forward_operators/DEFAULT_obs_def_mod.F90',
     output_obs_def_mod_file = '../../../observations/forward_operators/obs_def_mod.f90',
     input_obs_qty_mod_file  = '../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90',
     output_obs_qty_mod_file = '../../../assimilation_code/modules/observations/obs_kind_mod.f90',
     quantity_files          = '../../../assimilation_code/modules/observations/atmosphere_quantities_mod.f90',
     obs_type_files          = '../../../observations/forward_operators/obs_def_reanalysis_bufr_mod.f90',
                               '../../../observations/forward_operators/obs_def_rel_humidity_mod.f90',
                               '../../../observations/forward_operators/obs_def_altimeter_mod.f90'
    /

| 

.. container::

   +-------------------------+-----------------------+-------------------------------------------------------------------------------+
   | Item                    | Type                  | Description                                                                   |
   +=========================+=======================+===============================================================================+
   | input_obs_def_mod_file  | character(len=256)    | Path name of the template observation definition module                       |
   |                         |                       | to be preprocessed. The default is                                            |
   |                         |                       | ``../../../observations/forward_operators/DEFAULT_obs_def_mod.F90``.          |
   |                         |                       | This file must have the appropriate                                           |
   |                         |                       | commented lines indicating where the                                          |
   |                         |                       | different parts of the input special                                          |
   |                         |                       | obs definition modules are to be                                              |
   |                         |                       | inserted.                                                                     |
   +-------------------------+-----------------------+-------------------------------------------------------------------------------+
   | output_obs_def_mod_file | character(len=256)    | Path name of output observation definition module                             |
   |                         |                       | to be created by preprocess. The default is                                   |
   |                         |                       | ``../../../observations/forward_operators/obs_def_mod.f90``.                  |
   +-------------------------+-----------------------+-------------------------------------------------------------------------------+
   | input_obs_qty_mod_file  | character(len=256)    | Path name of input obs quantity file                                          |
   |                         |                       | to be preprocessed. The default path name is                                  |
   |                         |                       | ``../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90``. |
   |                         |                       | This file must have the appropriate                                           |
   |                         |                       | commented lines indicating where the                                          |
   |                         |                       | different quantity modules are to be                                          |
   |                         |                       | inserted.                                                                     |
   +-------------------------+-----------------------+-------------------------------------------------------------------------------+
   | output_obs_qty_mod_file | character(len=256)    | Path name of output obs quantity                                              |
   |                         |                       | module to be created by preprocess.                                           |
   |                         |                       | The default is                                                                |
   |                         |                       | ``../../../assimilation_code/modules/observations/obs_kind_mod.f90``.         |
   +-------------------------+-----------------------+-------------------------------------------------------------------------------+
   | obs_type_files          | character(len=256)(:) | A list of files containing                                                    |
   |                         |                       | observation definitions for the type                                          |
   |                         |                       | of observations you want to use with                                          |
   |                         |                       | DART. The maximum number of files is                                          |
   |                         |                       | limited to MAX_OBS_TYPE_FILES = 1000.                                         |
   |                         |                       | The DART obs_def files are in                                                 |
   |                         |                       | ``observations/forward_operators/obs_def_*.mod.f90``.                         |
   +-------------------------+-----------------------+-------------------------------------------------------------------------------+
   | overwrite_output        | logical               | By default, preprocess will overwrite                                         |
   |                         |                       | the existing obs_kind_mod.f90 and                                             |
   |                         |                       | obs_def_mod.f90 files. Set                                                    |
   |                         |                       | ``overwrite_output = .false.`` if you                                         |
   |                         |                       | want to preprocess to not overwrite                                           |
   |                         |                       | existing files.                                                               |
   +-------------------------+-----------------------+-------------------------------------------------------------------------------+

| 

Modules used
------------

::

   parse_arges_mod
   types_mod
   utilities_mod

Namelist interface ``&preprocess_nml`` must be read from file ``input.nml``.

Files
-----

-  input_obs_def_mod_file, specified by namelist; usually ``DEFAULT_obs_def_mod.F90``.
-  output_obs_def_mod_file, specified by namelist; usually ``obs_def_mod.f90``.
-  input_obs_qty_mod_file, specified by namelist; usually ``DEFAULT_obs_kind_mod.F90``.
-  output_obs_qty_mod_file, specified by namelist; usually ``obs_kind_mod.f90``.
-  obs_type_files, specified by namelist; usually files like ``obs_def_reanalysis_bufr_mod.f90``.
-  quantity_files, specified by namelist; usually files like ``atmosphere_quantities_mod.f90``.
-  namelistfile; ``input.nml``

References
----------

-  none
