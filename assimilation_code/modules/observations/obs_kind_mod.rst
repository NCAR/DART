MODULE ``obs_kind_mod``
=======================

Overview
--------

Introduction
^^^^^^^^^^^^

This module provides definitions of specific observation types and generic variable quantities, routines for mapping
between integer identifiers and string names, routines for reading and writing this information, and routines for
determining whether and how to process observations from an observation sequence file.

The distinction between quantities and types is this: ``Quantities`` apply both to observations and to state vector
variables. Knowing the type of an observation must be sufficient to compute the correct forward operator. The quantity
associated with an observation must be sufficient to identify which variable in the state vector should be used to
compute the expected value. ``Types`` only apply to observations, and are usually observation-platform dependent. Making
distinctions between different observation sources by using different types allows users to selectively assimilate,
evaluate, or ignore them.

Examples and use
^^^^^^^^^^^^^^^^

Generic quantities are associated with an observation type or with a model state variable. An example quantity is
``QTY_U_WIND_COMPONENT``. Multiple different specific observation types can be associated with this generic quantity,
for instance ``RADIOSONDE_U_WIND_COMPONENT``, ``ACARS_U_WIND_COMPONENT``, and ``SAT_U_WIND_COMPONENT``. Generic
quantities are defined via an integer parameter statement at the start of this module. As new generic quantities are
needed they are added to this list. Generic quantity integer parameters are required to start with ``QTY_`` and
observation types are NOT allowed to start with ``QTY_``.

Typically quantities are used by model-interface files ``models/xx/model_mod.f90``, observation forward operator files
``observations/forward_operators/obs_def_xx_mod.f90``, and observation converter programs
``observations/obs_converters/xx/xx.f90``.

The obs_kind module being described here is created by the program ``preprocess`` from two categories of input files.
First, a DEFAULT obs_kind module (normally called ``DEFAULT_obs_kind_mod.F90`` and documented in this directory) is used
as a template into which the preprocessor incorporates information from zero or more special obs_def modules (such as
``obs_def_1d_state_mod.f90`` or ``obs_def_reanalysis_bufr_mod.f90``) which are documented in the obs_def directory. If
no special obs_def files are included in the preprocessor namelist, a minimal ``obs_kind_mod.f90`` is created which can
only support identity forward observation operators.

All of the build scripts in DART remove the existing ``obs_kind_mod.f90`` file and regenerate it using the
``preprocess`` program. Do not add new quantities to ``obs_kind_mod.f90``, because these changes will not be kept when
you run *quickbuild.sh*.

Adding additional quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

New quantities should be added to a quantity file, for example a new ocean quantity should be added to
``ocean_quantities_mod.f90``. The quantity files are in ``assimilation_code/modules/observations/``.

Every line in a quantity file between the start and end markers must be a comment or a quantity definition (QTY_string).
Multiple name-value pairs can be specified for a quantity but are not required. For example, temperature may be defined:
``! QTY_TEMPERATURE units="K" minval=0.0``. Comments are allowed between quantity definitions or on the same line as the
definition. The code snippet below shows acceptable formats for quantity definitions

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
  ! 
  ! comment 
  ! 
  ! END DART PREPROCESS QUANTITY DEFINITIONS

Implementation details
^^^^^^^^^^^^^^^^^^^^^^

The obs_kind module contains an automatically-generated list of integer parameters, derived from the obs_def files, an
integer parameter ``max_defined_types_of_obs``, and an automatically-generated list of initializers for the
``obs_type_type`` derived type that defines the details of each observation type that has been created by the preprocess
program. Each entry contains the integer index of the observation type, the string name of the observation type (which
is identical to the F90 identifier), the integer index of the associated generic quantities, and three logicals
indicating whether this observation type is to be assimilated, evaluated only (forward operator is computed but not
assimilated), assimilated but has externally computed forward operator values in the input observation sequence file, or
ignored entirely. The logicals initially default to .false. and are set to .true. via the ``&obs_kind_nml`` namelist. A
second derived type ``obs_qty_type`` maps the integer parameter for a quantity to the quantity name (a string), and
stores any additional pair-value metadata for that quantity.

| 

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &obs_kind_nml
      assimilate_these_obs_types          = 'null',
      evaluate_these_obs_types            = 'null'
      use_precomputed_FOs_these_obs_types = 'null'
    /

| 

Controls what observation types are to be assimilated, evaluated, or ignored. For each entry, a list of observation type
names can be specified. Any name in the obs_type_type table is eligible. Specifying a name that is not in the table
results in an error. Specifying the same name for both namelist entries also results in an error. Observation types
specified in the list for assimilate_these_obs_types are assimilated. Those in the evaluate_these_obs_types list have
their forward operators computed and included in diagnostic files but are not assimilated. An observation type that is
specified in neither list is ignored. Identity observations, however, are always assimilated if present in the
obs_seq.out file.

.. container::

   +-------------------------------------+---------------------------------+--------------------------------------------+
   | Item                                | Type                            | Description                                |
   +=====================================+=================================+============================================+
   | assimilate_these_obs_types          | character(len=31), dimension(:) | Names of observation types to be           |
   |                                     |                                 | assimilated.                               |
   +-------------------------------------+---------------------------------+--------------------------------------------+
   | evaluate_these_obs_types            | character(len=31), dimension(:) | Names of observation types to be evaluated |
   |                                     |                                 | only.                                      |
   +-------------------------------------+---------------------------------+--------------------------------------------+
   | use_precomputed_FOs_these_obs_types | character(len=31), dimension(:) | If the forward operator values have been   |
   |                                     |                                 | precomputed outside of filter, for example |
   |                                     |                                 | for radiances or other compute intensive   |
   |                                     |                                 | computations, the ensemble of forward      |
   |                                     |                                 | operator values can be stored in the       |
   |                                     |                                 | observation sequence file. For any type    |
   |                                     |                                 | listed here, the forward operator          |
   |                                     |                                 | interpolation code will not be called and  |
   |                                     |                                 | the values in the file will be used        |
   |                                     |                                 | instead.                                   |
   +-------------------------------------+---------------------------------+--------------------------------------------+

For example:

::

   &obs_kind_nml
      assimilate_these_obs_types = 'RADIOSONDE_TEMPERATURE',
                                   'RADIOSONDE_U_WIND_COMPONENT',
                                   'RADIOSONDE_V_WIND_COMPONENT',
      evaluate_these_obs_types   = 'RADIOSONDE_SURFACE_PRESSURE',
     use_precomputed_FOs_these_obs_types = 'RADIANCE'
   /

| would assimilate temperature and wind observations, but only compute the forward operators for surface pressure obs.
  Radiance observations have precomputed values for each ensemble member in the input observation sequence file which
  would be used instead of calling the forward operator code.

Modules used
------------

::

   utilities_mod

| 

Public interfaces
-----------------

========================= ============================
*use obs_def_mod, only :* max_defined_types_of_obs
\                         get_num_types_of_obs
\                         get_num_quantities
\                         get_name_for_type_of_obs
\                         get_name_for_quantity
\                         get_index_for_type_of_obs
\                         get_index_for_quantity
\                         assimilate_this_type_of_obs
\                         evaluate_this_type_of_obs
\                         get_quantity_for_type_of_obs
\                         write_type_of_obs_table
\                         read_type_of_obs_table
\                         get_type_of_obs_from_menu
\                         map_type_of_obs_table
\                         GENERIC_QTY_DEFINITIONS
\                         OBSERVATION_TYPES
========================= ============================

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *integer, parameter :: max_defined_types_of_obs*

.. container:: indent1

   The total number of available observation types in the obs_type_type table. This value is added by the preprocess
   program and depends on which ``obs_def_xxx_mod.f90`` files are listed in the
   `&preprocess_nml <../../programs/preprocess/preprocess.html#Namelist>`__ namelist.

   There is also a function interface which is an alternate method to get this value. In some cases the code requires a
   parameter value known at compile time (for declaring a fixed length array, for example). For an array allocated at
   run time the size can be returned by the function interface.

| 

.. container:: routine

   *var = get_num_types_of_obs()*
   ::

      integer :: get_num_types_of_obs

.. container:: indent1

   Returns the number of different specific observation types (e.g. RADIOSONDE_TEMPERATURE, AIRCRAFT_SPECIFIC_HUMIDITY)
   defined in the obs_kind_mod.f90 file. This file is generated by the preprocess program. This is the same value as the
   public 'max_defined_types_of_obs' above.

   ======= =========================================================================================
   ``var`` Integer count of the total number of specific types defined in the obs_kind_mod.f90 file.
   ======= =========================================================================================

| 

.. container:: routine

   *var = get_num_quantities()*
   ::

      integer :: get_num_quantities

.. container:: indent1

   Returns the number of different generic quantities (e.g. QTY_TEMPERATURE, QTY_SPECIFIC_HUMIDITY) defined in the
   obs_kind_mod.f90 file. This file is generated by the preprocess program.

   ======= =============================================================================================
   ``var`` Integer count of the total number of generic quantities defined in the obs_kind_mod.f90 file.
   ======= =============================================================================================

| 

.. container:: routine

   *var = get_name_for_type_of_obs(obs_type_ind)*
   ::

      character(len=32)              :: get_name_for_type_of_obs
      integer, intent(in)            :: obs_type_ind

.. container:: indent1

   Given an integer index return the string name of the corresponding specific observation type (e.g.
   "RADIOSONDE_TEMPERATURE", "AIRCRAFT_SPECIFIC_HUMIDITY"). This string is the same as the F90 identifier associated
   with the integer index.

   ================ ==================================================================
   ``var``          Name string associated with this entry in the obs_type_type table.
   ``obs_type_ind`` An integer index into the obs_type_type table.
   ================ ==================================================================

| 

.. container:: routine

   *var = get_name_for_quantity(obs_qty_ind)*
   ::

      character(len=32)              :: get_name_for_quantity
      integer, intent(in)            :: obs_qty_ind

.. container:: indent1

   Given an integer index return the string name of the corresponding generic quantity (e.g. "QTY_TEMPERATURE",
   "QTY_SPECIFIC_HUMIDITY"). This string is the same as the F90 identifier associated with the integer index.

   =============== =================================================================
   ``var``         Name string associated with this entry in the obs_qty_type table.
   ``obs_qty_ind`` An integer index into the obs_qty_type table.
   =============== =================================================================

| 

.. container:: routine

   *var = get_index_for_type_of_obs(obs_type_name)*
   ::

      integer                       :: get_index_for_type_of_obs
      character(len=*), intent(in)  :: obs_type_name

.. container:: indent1

   Given the name of a specific observation type (e.g. "RADIOSONDE_TEMPERATURE", "AIRCRAFT_SPECIFIC_HUMIDITY"), returns
   the index of the entry in the obs_type_type table with this name. If the name is not found in the table, a -1 is
   returned. The integer returned for a successful search is the value of the integer parameter with the same identifier
   as the name string.

   +-------------------------------+-------------------------------------------------------------------------------------+
   | ``get_index_for_type_of_obs`` | Integer index into the obs_type_type table entry with name string corresponding to  |
   |                               | obs_type_name.                                                                      |
   +-------------------------------+-------------------------------------------------------------------------------------+
   | ``obs_type_name``             | Name of specific observation type found in obs_type_type table.                     |
   +-------------------------------+-------------------------------------------------------------------------------------+

| 

.. container:: routine

   *var = get_index_for_quantity(obs_qty_name)*
   ::

      integer                       :: get_index_for_quantity
      character(len=32), intent(in) :: obs_qty_name

.. container:: indent1

   Given the name of a generic quantity (e.g. "QTY_TEMPERATURE", "QTY_SPECIFIC_HUMIDITY"), returns the index of the
   entry in the obs_qty_type table with this name. If the name is not found in the table, a -1 is returned. The integer
   returned for a successful search is the value of the integer parameter with the same identifier as the name string.

   +----------------------------+----------------------------------------------------------------------------------------+
   | ``get_index_for_quantity`` | Integer index into the obs_qty_type table entry with name string corresponding to      |
   |                            | obs_qty_name.                                                                          |
   +----------------------------+----------------------------------------------------------------------------------------+
   | ``obs_qty_name``           | Name of generic kind found in obs_qty_type table.                                      |
   +----------------------------+----------------------------------------------------------------------------------------+

| 

.. container:: routine

   *var = assimilate_this_type_of_obs(obs_type_ind)*
   ::

      logical              :: assimilate_this_type_of_obs
      integer, intent(in)  :: obs_type_ind

.. container:: indent1

   Given the integer index associated with a specific observation type (e.g. RADIOSONDE_TEMPERATURE,
   AIRCRAFT_SPECIFIC_HUMIDITY), return true if this observation type is to be assimilated, otherwise false. The
   parameter defined by this name is used as an integer index into the obs_type_type table to return the status of this
   type.

   ================ ===========================================================================
   ``var``          Returns true if this entry in the obs_type_type table is to be assimilated.
   ``obs_type_ind`` An integer index into the obs_type_type table.
   ================ ===========================================================================

| 

.. container:: routine

   *var = evaluate_this_type_of_obs(obs_type_ind)*
   ::

      logical              :: evaluate_this_type_of_obs
      integer, intent(in)  :: obs_type_ind

.. container:: indent1

   Given the integer index associated with a specific observation type (e.g. RADIOSONDE_TEMPERATURE,
   AIRCRAFT_SPECIFIC_HUMIDITY), return true if this observation type is to be evaluated only, otherwise false. The
   parameter defined by this name is used as an integer index into the obs_type_type table to return the status of this
   type.

   ================ =========================================================================
   ``var``          Returns true if this entry in the obs_type_type table is to be evaluated.
   ``obs_type_ind`` An integer index into the obs_type_type table.
   ================ =========================================================================

| 

.. container:: routine

   *var = get_quantity_for_type_of_obs(obs_type_ind)*
   ::

      integer              :: get_quantity_for_type_of_obs
      integer, intent(in)  :: obs_type_ind

.. container:: indent1

   Given the integer index associated with a specific observation type (e.g. RADIOSONDE_TEMPERATURE,
   AIRCRAFT_SPECIFIC_HUMIDITY), return the generic quantity associated with this type (e.g. QTY_TEMPERATURE,
   QTY_SPECIFIC_HUMIDITY). The parameter defined by this name is used as an integer index into the obs_type_type table
   to return the generic quantity associated with this type.

   ================ =========================================================================
   ``var``          Returns the integer GENERIC quantity index associated with this obs type.
   ``obs_type_ind`` An integer index into the obs_type_type table.
   ================ =========================================================================

| 

.. container:: routine

   *call write_type_of_obs_table(ifile [, fform, use_list])*
   ::

      integer,                    intent(in) :: ifile
      character(len=*), optional, intent(in) :: fform
      integer,          optional, intent(in) :: use_list(:)

.. container:: indent1

   Writes out information about all defined observation types from the obs_type_type table. For each entry in the table,
   the integer index of the observation type and the associated string are written. These appear in the header of an
   obs_sequence file. If given, the *use_list(:)* must be the same length as the max_obs_specific count. If greater than
   0, the corresponding index will be written out; if 0 this entry is skipped. This allows a table of contents to be
   written which only includes those types actually being used.

   +---------------+-----------------------------------------------------------------------------------------------------+
   | ``ifile``     | Unit number of output observation sequence file being written.                                      |
   +---------------+-----------------------------------------------------------------------------------------------------+
   | *fform*       | Optional format for file. Default is FORMATTED.                                                     |
   +---------------+-----------------------------------------------------------------------------------------------------+
   | *use_list(:)* | Optional integer array the same length as the number of specific types (from get_num_types_of_obs() |
   |               | or the public max_defined_types_of_obs). If value is larger than 0, the corresponding type          |
   |               | information will be written out. If 0, it will be skipped. If this argument is not specified, all   |
   |               | values will be written.                                                                             |
   +---------------+-----------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call read_type_of_obs_table(ifile, pre_I_format [, fform])*
   ::

      integer,                    intent(in) :: ifile
      logical,                    intent(in) :: pre_I_format !(deprecated)
      character(len=*), optional, intent(in) :: fform

.. container:: indent1

   Reads the mapping between integer indices and observation type names from the header of an observation sequence file
   and prepares mapping to convert these to values defined in the obs_type_type table. If pre_I_format is true, there is
   no header in the observation sequence file and it is assumed that the integer indices for observation types in the
   file correspond to the storage order of the obs_type_type table (integer index 1 in the file corresponds to the first
   table entry, etc.) Support for pre_I_format is deprecated and may be dropped in future releases of DART.

   ================ ===========================================================================
   ``ifile``        Unit number of output observation sequence file being written.
   ``pre_I_format`` True if the file being read has no obs type definition header (deprecated).
   *fform*          Optional format for file. Default is FORMATTED.
   ================ ===========================================================================

| 

.. container:: routine

   *var = get_type_of_obs_from_menu()*
   ::

      integer              :: get_type_of_obs_from_menu

.. container:: indent1

   Interactive input of observation type. Prompts user with list of available types and validates entry before
   returning.

   ======= ==================================
   ``var`` Integer index of observation type.
   ======= ==================================

| 

.. container:: routine

   *var = map_type_of_obs_table(obs_def_index)*
   ::

      integer              :: map_type_of_obs_table
      integer, intent(in)  :: obs_def_index

.. container:: indent1

   Maps from the integer observation type index in the header block of an input observation sequence file into the
   corresponding entry in the obs_type_type table. This allows observation sequences that were created with different
   obs_kind_mod.f90 versions to be used with the current obs_kind_mod.

   ================= ===============================================================
   ``var``           Index of this observation type in obs_type_type table.
   ``obs_def_index`` Index of observation type from input observation sequence file.
   ================= ===============================================================

| 

.. container:: routine

   ``integer, parameter :: QTY_.....``

.. container:: indent1

   All generic quantities available are public parameters that begin with ``QTY_``.

| 

.. container:: routine

   *integer, parameter :: SAMPLE_OBS_TYPE*

.. container:: indent1

   A list of all observation types that are available is provided as a set of integer parameter statements. The F90
   identifiers are the same as the string names that are associated with this identifier in the obs_type_type table.

| 

Files
-----

-  &obs_kind_nml in input.nml
-  Files containing input or output observation sequences.

| 

References
----------

-  none

| 
