MODULE model_mod
================

Overview
--------

| Every model that is DART compliant must provide an set of interfaces that will be called by DART code. For models
  which have no special code for some of these routines, they can pass through the call to this default module, which
  satisfies the call but does no work. To use these routines in a ``model_mod.f90``, add at the top:

::

   use default_model_mod, only : xxx, yyy

and then leave them in the public list.

Namelist
--------

The default routines have no namelist.

Other modules used
------------------

::

   types_mod
   time_manager_mod
   location_mod
   utilities_mod
   netcdf_utilities_mod
   ensemble_manager_mod
   dart_time_io_mod

Public interfaces
-----------------

======================= ===================================
*use model_mod, only :* get_model_size
\                       adv_1step
\                       get_state_meta_data
\                       model_interpolate
\                       shortest_time_between_assimilations
\                       static_init_model
\                       init_time
\                       fail_init_time
\                       init_conditions
\                       fail_init_conditions
\                       nc_write_model_atts
\                       nc_write_model_vars
\                       pert_model_copies
\                       get_close_obs
\                       get_close_state
\                       convert_vertical_obs
\                       convert_vertical_state
\                       read_model_time
\                       write_model_time
\                       end_model
======================= ===================================

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *model_size = get_model_size( )*
   ::

      integer(i8) :: get_model_size

.. container:: indent1

   Returns the length of the model state vector as 1. Probably not what you want. The model_mod should set this to the
   right size and not use this routine.

   ============== =====================================
   ``model_size`` The length of the model state vector.
   ============== =====================================

| 

.. container:: routine

   *call adv_1step(x, time)*
   ::

      real(r8), dimension(:), intent(inout) :: x
      type(time_type),        intent(in)    :: time

.. container:: indent1

   Throws a fatal error. If the model_mod can advance the model it should provide a real routine. This default routine
   is intended for use by models which cannot advance themselves from inside filter.

   ======== ==================================
   ``x``    State vector of length model_size.
   ``time`` Current time of model state.
   ======== ==================================

| 

.. container:: routine

   *call get_state_meta_data (index_in, location, [, var_type] )*
   ::

      integer,             intent(in)  :: index_in
      type(location_type), intent(out) :: location
      integer, optional,   intent(out) ::  var_type 

.. container:: indent1

   Sets the location to missing and the variable type to 0. The model_mod should provide a routine that sets a real
   location and a state vector type for the requested item in the state vector.

   ============ ===================================================================
   ``index_in`` Index of state vector element about which information is requested.
   ``location`` The location of state variable element.
   *var_type*   The generic quantity of the state variable element.
   ============ ===================================================================

| 

.. container:: routine

   *call model_interpolate(state_handle, ens_size, location, obs_quantity, expected_obs, istatus)*
   ::

      type(ensemble_type),    intent(in)  :: state_handle
      integer,                intent(in)  :: ens_size
      type(location_type),    intent(in)  :: location
      integer,                intent(in)  :: obs_quantity
      real(r8),               intent(out) :: expected_obs(ens_size)
      integer,                intent(out) :: istatus(ens_size)

.. container:: indent1

   Sets the expected obs to missing and returns an error code for all obs. This routine should be supplied by the
   model_mod.

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``state_handle`` | The handle to the state structure containing information about the state vector about which      |
   |                  | information is requested.                                                                        |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``ens_size``     | The ensemble size.                                                                               |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``location``     | Location to which to interpolate.                                                                |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``obs_quantity`` | Quantity of state field to be interpolated.                                                      |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``expected_obs`` | The interpolated values from the model.                                                          |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``istatus``      | Integer values return 0 for success. Other positive values can be defined for various failures.  |
   +------------------+--------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *var = shortest_time_between_assimilations()*
   ::

      type(time_type) :: shortest_time_between_assimilations

.. container:: indent1

   Returns 1 day.

   ======= ===================================
   ``var`` Smallest advance time of the model.
   ======= ===================================

| 

.. container:: routine

   *call static_init_model()*

.. container:: indent1

   Does nothing.

| 

.. container:: routine

   *call init_time(time)*
   ::

      type(time_type), intent(out) :: time

.. container:: indent1

   Returns a time of 0.

   ======== ===================
   ``time`` Initial model time.
   ======== ===================

| 

.. container:: routine

   *call fail_init_time(time)*
   ::

      type(time_type), intent(out) :: time

.. container:: indent1

   Throws a fatal error. This is appropriate for models that cannot start from arbitrary initial conditions.

   ======== ============================
   ``time`` NOT SET. Initial model time.
   ======== ============================

| 

.. container:: routine

   *call init_conditions(x)*
   ::

      real(r8), dimension(:), intent(out) :: x

.. container:: indent1

   Returns x(:) = 0.0

   ===== ====================================
   ``x`` Initial conditions for state vector.
   ===== ====================================

| 

.. container:: routine

   *call fail_init_conditions(x)*
   ::

      real(r8), dimension(:), intent(out) :: x

.. container:: indent1

   Throws a fatal error. This is appropriate for models that cannot start from arbitrary initial conditions.

   ===== =============================================
   ``x`` NOT SET: Initial conditions for state vector.
   ===== =============================================

| 

.. container:: routine

   *call nc_write_model_atts(ncFileID, domain_id)*
   ::

      integer, intent(in) :: ncFileID
      integer, intent(in) :: domain_id

.. container:: indent1

   Does nothing.

   +---------------+-----------------------------------------------------------------------------------------------------+
   | ``ncFileID``  | Integer file descriptor to previously-opened netCDF file.                                           |
   +---------------+-----------------------------------------------------------------------------------------------------+
   | ``domain_id`` | integer describing the domain (which can be a nesting level, a component model ...) Models with     |
   |               | nested grids are decomposed into 'domains' in DART. The concept is extended to refer to 'coupled'   |
   |               | models where one model component may be the atmosphere, another component may be the ocean, or      |
   |               | land, or ionosphere ... these would be referenced as different domains.                             |
   +---------------+-----------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call nc_write_model_vars(ncFileID, domain_id, state_ens_handle [, memberindex] [, timeindex])*
   ::

      integer,             intent(in) :: ncFileID
      integer,             intent(in) :: domain_id
      type(ensemble_type), intent(in) :: state_ens_handle
      integer, optional,   intent(in) :: memberindex
      integer, optional,   intent(in) :: timeindex

.. container:: indent1

   Does nothing

   +----------------------+----------------------------------------------------------------------------------------------+
   | ``ncFileID``         | file descriptor to previously-opened netCDF file.                                            |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``domain_id``        | integer describing the domain (which can be a nesting level, a component model ...)          |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``state_ens_handle`` | The handle to the state structure containing information about the state vector about which  |
   |                      | information is requested.                                                                    |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``memberindex``      | Integer index of ensemble member to be written.                                              |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``timeindex``        | The timestep counter for the given state.                                                    |
   +----------------------+----------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)*
   ::

      type(ensemble_type), intent(inout) :: state_ens_handle
      integer,             intent(in)    :: ens_size
      real(r8),            intent(in)    :: pert_amp
      logical,             intent(out)   :: interf_provided

.. container:: indent1

   Returns 'interface provided' flag as false, so the default perturb routine in DART will add small amounts of gaussian
   noise to all parts of the state vector.

   +----------------------+----------------------------------------------------------------------------------------------+
   | ``state_ens_handle`` | The handle containing an ensemble of state vectors to be perturbed.                          |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``ens_size``         | The number of ensemble members to perturb.                                                   |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``pert_amp``         | the amplitude of the perturbations. The interpretation is based on the model-specific        |
   |                      | implementation.                                                                              |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``interf_provided``  | Returns false if model_mod cannot do this, else true.                                        |
   +----------------------+----------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, num_close, close_ind [, dist] [,
   state_handle)*
   ::

      type(get_close_type),          intent(in)  :: gc
      type(location_type),           intent(in)  :: base_loc
      integer,                       intent(in)  :: base_type
      type(location_type),           intent(in)  :: locs(:)
      integer,                       intent(in)  :: loc_qtys(:)
      integer,                       intent(in)  :: loc_types(:)
      integer,                       intent(out) :: num_close
      integer,                       intent(out) :: close_ind(:)
      real(r8),            optional, intent(out) :: dist(:)
      type(ensemble_type), optional, intent(in)  :: state_handle

.. container:: indent1

   Passes the call through to the location module code.

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``gc``           | The get_close_type which stores precomputed information about the locations to speed up          |
   |                  | searching                                                                                        |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``base_loc``     | Reference location. The distances will be computed between this location and every other         |
   |                  | location in the obs list                                                                         |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``base_type``    | The DART quantity at the ``base_loc``                                                            |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``locs(:)``      | Compute the distance between the ``base_loc`` and each of the locations in this list             |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``loc_qtys(:)``  | The corresponding quantity of each item in the ``locs`` list                                     |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``loc_types(:)`` | The corresponding type of each item in the ``locs`` list. This is not available in the default   |
   |                  | implementation but may be used in custom implementations.                                        |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``num_close``    | The number of items from the ``locs`` list which are within maxdist of the base location         |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``close_ind(:)`` | The list of index numbers from the ``locs`` list which are within maxdist of the base location   |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``dist(:)``      | If present, return the distance between each entry in the close_ind list and the base location.  |
   |                  | If not present, all items in the obs list which are closer than maxdist will be added to the     |
   |                  | list but the overhead of computing the exact distances will be skipped.                          |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``state_handle`` | The handle to the state structure containing information about the state vector about which      |
   |                  | information is requested.                                                                        |
   +------------------+--------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call get_close_state(gc, base_loc, base_type, state_loc, state_qtys, state_indx, num_close, close_ind, dist,
   state_handle*)
   ::

      type(get_close_type), intent(in)    :: gc
      type(location_type),  intent(inout) :: base_loc
      integer,              intent(in)    :: base_type
      type(location_type),  intent(inout) :: state_loc(:)
      integer,              intent(in)    :: state_qtys(:)
      integer(i8),          intent(in)    :: state_indx(:)
      integer,              intent(out)   :: num_close
      integer,              intent(out)   :: close_ind(:)
      real(r8),             intent(out)   :: dist(:)
      type(ensemble_type),  intent(in)    :: state_handle

.. container:: indent1

   Passes the call through to the location module code.

   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``gc``            | The get_close_type which stores precomputed information about the locations to speed up         |
   |                   | searching                                                                                       |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``base_loc``      | Reference location. The distances will be computed between this location and every other        |
   |                   | location in the obs list                                                                        |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``base_type``     | The DART quantity at the ``base_loc``                                                           |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``state_loc(:)``  | Compute the distance between the ``base_loc`` and each of the locations in this list            |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``state_qtys(:)`` | The corresponding quantity of each item in the ``state_loc`` list                               |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``state_indx(:)`` | The corresponding DART index of each item in the ``state_loc`` list. This is not available in   |
   |                   | the default implementation but may be used in custom implementations.                           |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``num_close``     | The number of items from the ``state_loc`` list which are within maxdist of the base location   |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``close_ind(:)``  | The list of index numbers from the ``state_loc`` list which are within maxdist of the base      |
   |                   | location                                                                                        |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``dist(:)``       | If present, return the distance between each entry in the ``close_ind`` list and the base       |
   |                   | location. If not present, all items in the ``state_loc`` list which are closer than maxdist     |
   |                   | will be added to the list but the overhead of computing the exact distances will be skipped.    |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``state_handle``  | The handle to the state structure containing information about the state vector about which     |
   |                   | information is requested.                                                                       |
   +-------------------+-------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call convert_vertical_obs(state_handle, num, locs, loc_qtys, loc_types, which_vert, status)*
   ::

      type(ensemble_type), intent(in)  :: state_handle
      integer,             intent(in)  :: num
      type(location_type), intent(in)  :: locs(:)
      integer,             intent(in)  :: loc_qtys(:)
      integer,             intent(in)  :: loc_types(:)
      integer,             intent(in)  :: which_vert
      integer,             intent(out) :: status(:)

.. container:: indent1

   Passes the call through to the location module code.

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``state_handle`` | The handle to the state.                                                                         |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``num``          | the number of observation locations                                                              |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``locs``         | the array of observation locations                                                               |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``loc_qtys``     | the array of observation quantities.                                                             |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``loc_types``    | the array of observation types.                                                                  |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``which_vert``   | the desired vertical coordinate system. There is a table in the ``location_mod.f90`` that        |
   |                  | relates integers to vertical coordinate systems.                                                 |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``status``       | Success or failure of the vertical conversion. If ``istatus = 0``, the conversion was a success. |
   |                  | Any other value is a failure.                                                                    |
   +------------------+--------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call convert_vertical_state(state_handle, num, locs, loc_qtys, loc_indx, which_vert, status)*
   ::

      type(ensemble_type), intent(in)  :: state_handle
      integer,             intent(in)  :: num
      type(location_type), intent(in)  :: locs(:)
      integer,             intent(in)  :: loc_qtys(:)
      integer(i8),         intent(in)  :: loc_indx(:)
      integer,             intent(in)  :: which_vert
      integer,             intent(out) :: status(:)

.. container:: indent1

   Passes the call through to the location module code.

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``state_handle`` | The handle to the state.                                                                         |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``num``          | the number of state locations                                                                    |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``locs``         | the array of state locations                                                                     |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``loc_qtys``     | the array of state quantities.                                                                   |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``loc_indx``     | the array of state vector indices.                                                               |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``which_vert``   | the desired vertical coordinate system. There is a table in the ``location_mod.f90`` that        |
   |                  | relates integers to vertical coordinate systems.                                                 |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``status``       | Success or failure of the vertical conversion. If ``istatus = 0``, the conversion was a success. |
   |                  | Any other value is a failure.                                                                    |
   +------------------+--------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *model_time = read_model_time(filename)*
   ::

      character(len=*), intent(in) :: filename
      type(time_type)              :: model_time

.. container:: indent1

   Passes the call through to the dart_time_io module code.

   ============== ====================================
   ``filename``   netCDF file name
   ``model_time`` The current time of the model state.
   ============== ====================================

| 

.. container:: routine

   *call write_model_time(ncid, dart_time)*
   ::

      integer,          intent(in) :: ncid
      type(time_type),  intent(in) :: dart_time

.. container:: indent1

   Passes the call through to the dart_time_io module code.

   ============= ====================================
   ``ncid``      handle to an open netCDF file
   ``dart_time`` The current time of the model state.
   ============= ====================================

| 

.. container:: routine

   *call end_model()*

.. container:: indent1

   Does nothing.

Files
-----

none

References
----------

#. none

Private components
------------------

N/A
