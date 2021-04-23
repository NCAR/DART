pe2lyr
======

.. attention::

   ``pe2lyr`` works with versions of DART *before* Manhattan (9.x.x) and has yet to be updated. If you are interested in
   using ``pe2lyr`` with more recent versions of DART, contact DAReS staff to assess the feasibility of an update.
   Until that time, you should consider this documentation as out-of-date.


Overview
--------

DART standard interfaces for a two-layer isentropic primitive equation model.

The 16 public interfaces are standardized for all DART compliant models. These interfaces allow DART to advance the
model, get the model state and metadata describing this state, find state variables that are close to a given
location, and do spatial interpolation for model state variables.

This model is a 2-layer, isentropic, primitive equation model on a sphere. TODO: add more detail here, including
equations, etc.

Contact: Jeffrey.S.Whitaker@noaa.gov

Other modules used
------------------

::

   types_mod
   time_manager_mod
   utilities_mod
   random_seq_mod
   threed_sphere/location_mod

Public interfaces
-----------------

======================= ======================
*use model_mod, only :* get_model_size
\                       adv_1step
\                       get_state_meta_data
\                       model_interpolate
\                       get_model_time_step
\                       static_init_model
\                       end_model
\                       init_time
\                       init_conditions
\                       nc_write_model_atts
\                       nc_write_model_vars
\                       pert_model_state
\                       get_close_maxdist_init
\                       get_close_obs_init
\                       get_close_obs
\                       ens_mean_for_model
======================= ======================

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *model_size = get_model_size( )*
   ::

      integer :: get_model_size

.. container:: indent1

   Returns the size of the model as an integer. For this model the default grid size is 96 (lon) by 48 (lat) by 2
   levels, and 3 variables (U, V, Z) at each grid location, for a total size of 27,648. There are alternative include
   files which, if included at compile time instead of the default file, defines a grid at twice and 4 times this
   resolution. They have corresponding truncation values of T63 and T127 (the default grid uses T31).

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

   Advances the model for a single time step. The time associated with the initial model state is also input although it
   is not used for the computation.

   ======== ==========================================
   ``x``    State vector of length model_size.
   ``time`` Specifies time of the initial model state.
   ======== ==========================================

| 

.. container:: routine

   *call get_state_meta_data (index_in, location, [, var_type] )*
   ::

      integer,             intent(in)  :: index_in
      type(location_type), intent(out) :: location
      integer, optional,   intent(out) ::  var_type 

.. container:: indent1

   | Returns metadata about a given element, indexed by ``index_in``, in the model state vector. The ``location``
     defines where the state variable is located.
   | For this model, the default grid is a global lat/lon grid, 96 (lon) by 48 (lat) by 2 levels. The variable types are
     U, V, and Z:

   -  1 = TYPE_u
   -  2 = TYPE_v
   -  901 = TYPE_z

   Grids at twice and 4 times the resolution can be compiled in instead by using one of the alternative header files
   (see ``resolt31.h`` (the default), ``resolt63.h``, and ``resolt127.h``).

   ============ ===================================================================
   ``index_in`` Index of state vector element about which information is requested.
   ``location`` The location of state variable element.
   *var_type*   The type of the state variable element.
   ============ ===================================================================

| 

.. container:: routine

   *call model_interpolate(x, location, itype, obs_val, istatus)*
   ::

      real(r8), dimension(:), intent(in)  :: x
      type(location_type),    intent(in)  :: location
      integer,                intent(in)  :: itype
      real(r8),               intent(out) :: obs_val
      integer,                intent(out) :: istatus

.. container:: indent1

   Given a state vector, a location, and a model state variable type, interpolates the state variable field to that
   location and returns the value in obs_val. The istatus variable is always returned as 0 (OK).

   ============ ===========================================================================================
   ``x``        A model state vector.
   ``location`` Location to which to interpolate.
   ``itype``    Type of state field to be interpolated.
   ``obs_val``  The interpolated value from the model.
   ``istatus``  Integer value returning 0 for successful, other values can be defined for various failures.
   ============ ===========================================================================================

| 

.. container:: routine

   *var = get_model_time_step()*
   ::

      type(time_type) :: get_model_time_step

.. container:: indent1

   Returns the the time step of the model; the smallest increment in time that the model is capable of advancing the
   state in a given implementation. For this model the default value is 20 minutes (1200 seconds), but also comes with
   header files with times steps of 10 and 5 minutes (for higher grid resolution and truncation constants).

   ======= ============================
   ``var`` Smallest time step of model.
   ======= ============================

| 

.. container:: routine

   *call static_init_model()*

.. container:: indent1

   | Used for runtime initialization of a model, for instance calculating storage requirements, initializing model
     parameters, etc. This is the first call made to a model by any DART compliant assimilation routines.
   | In this model, it allocates space for the grid, and initializes the grid locations, data values, and various
     parameters, including spherical harmonic weights.

| 

.. container:: routine

   *call end_model()*

.. container:: indent1

   A stub since the pe2lyr model does no cleanup.

| 

.. container:: routine

   *call init_time(time)*
   ::

      type(time_type), intent(out) :: time

.. container:: indent1

   Returns the time at which the model will start if no input initial conditions are to be used. This model sets the
   time to 0.

   ======== ===================
   ``time`` Initial model time.
   ======== ===================

| 

.. container:: routine

   *call init_conditions(x)*
   ::

      real(r8), dimension(:), intent(out) :: x

.. container:: indent1

   Returns default initial conditions for model; generally used for spinning up initial model states. This model sets
   the default state vector based on the initialized fields in the model. (TODO: which are what?)

   ===== ====================================
   ``x`` Initial conditions for state vector.
   ===== ====================================

| 

.. container:: routine

   *ierr = nc_write_model_atts(ncFileID)*
   ::

      integer             :: nc_write_model_atts
      integer, intent(in) :: ncFileID

.. container:: indent1

   This routine writes the model-specific attributes to a netCDF file. This includes coordinate variables and any
   metadata, but NOT the model state vector. This model writes out the data as U, V, and Z arrays on a lat/lon/height
   grid, so the attributes are organized in the same way.

   ============ =========================================================
   ``ncFileID`` Integer file descriptor to previously-opened netCDF file.
   ``ierr``     Returns a 0 for successful completion.
   ============ =========================================================

| 

.. container:: routine

   *ierr = nc_write_model_vars(ncFileID, statevec, copyindex, timeindex)*
   ::

      integer                            :: nc_write_model_vars
      integer,                intent(in) :: ncFileID
      real(r8), dimension(:), intent(in) :: statevec
      integer,                intent(in) :: copyindex
      integer,                intent(in) :: timeindex

.. container:: indent1

   This routine writes the model-specific state vector (data) to a netCDF file. This model writes out the data as U, V,
   and Z arrays on a lat/lon/height grid.

   ============= =================================================
   ``ncFileID``  file descriptor to previously-opened netCDF file.
   ``statevec``  A model state vector.
   ``copyindex`` Integer index of copy to be written.
   ``timeindex`` The timestep counter for the given state.
   ``ierr``      Returns 0 for normal completion.
   ============= =================================================

| 

.. container:: routine

   *call pert_model_state(state, pert_state, interf_provided)*
   ::

      real(r8), dimension(:), intent(in)  :: state
      real(r8), dimension(:), intent(out) :: pert_state
      logical,                intent(out) :: interf_provided

.. container:: indent1

   Given a model state vector, perturbs this vector. Used to generate initial conditions for spinning up ensembles. This
   model has no code to generate these values, so it returns ``interf_provided`` as .false. and the default algorithms
   in filter are then used by the calling code.

   =================== =============================================
   ``state``           State vector to be perturbed.
   ``pert_state``      Perturbed state vector
   ``interf_provided`` Returned false; interface is not implemented.
   =================== =============================================

| 

.. container:: routine

   *call get_close_maxdist_init(gc, maxdist)*
   ::

      type(get_close_type), intent(inout) :: gc
      real(r8),             intent(in)    :: maxdist

.. container:: indent1

   In distance computations any two locations closer than the given ``maxdist`` will be considered close by the
   ``get_close_obs()`` routine. Pass-through to the 3-D sphere locations module. See
   `get_close_maxdist_init() <../../location/threed_sphere/location_mod.html#get_close_maxdist_init>`__ for the
   documentation of this subroutine.

   =========== =================================================================================================
   ``gc``      The get_close_type which stores precomputed information about the locations to speed up searching
   ``maxdist`` Anything closer than this will be considered close.
   =========== =================================================================================================

| 

.. container:: routine

   *call get_close_obs_init(gc, num, obs)*
   ::

      type(get_close_type), intent(inout) :: gc
      integer,              intent(in)    :: num
      type(location_type),  intent(in)    :: obs(num)

.. container:: indent1

   Pass-through to the 3-D sphere locations module. See
   `get_close_obs_init() <../../location/threed_sphere/location_mod.html#get_close_obs_init>`__ for the documentation of
   this subroutine.

| 

.. container:: routine

   *call get_close_obs(gc, base_obs_loc, base_obs_kind, obs, obs_kind, num_close, close_ind [, dist])*
   ::

      type(get_close_type), intent(in)  :: gc
      type(location_type),  intent(in)  :: base_obs_loc
      integer,              intent(in)  :: base_obs_kind
      type(location_type),  intent(in)  :: obs(:)
      integer,              intent(in)  :: obs_kind(:)
      integer,              intent(out) :: num_close
      integer,              intent(out) :: close_ind(:)
      real(r8), optional,   intent(out) :: dist(:)

.. container:: indent1

   | Given a location and kind, compute the distances to all other locations in the ``obs`` list. The return values are
     the number of items which are within maxdist of the base, the index numbers in the original obs list, and
     optionally the distances. The ``gc`` contains precomputed information to speed the computations.
   | Pass-through to the 3-D sphere locations module. See
     `get_close_obs() <../../location/threed_sphere/location_mod.html#get_close_obs>`__ for the documentation of this
     subroutine.

| 

.. container:: routine

   *call ens_mean_for_model(ens_mean)*
   ::

      real(r8), dimension(:), intent(in) :: ens_mean

.. container:: indent1

   Stub only. Not needed by this model.

   ============ ==========================================
   ``ens_mean`` State vector containing the ensemble mean.
   ============ ==========================================

| 

This model currently has no values settable by namelist.

Files
-----

-  The model source is in pe2lyr_mod.f90, and the spherical harmonic code is in spharmt_mod.f90. The various resolution
   settings are in resolt31.h, resolt63.h, and resolt127.h.

References
----------

Zou, X., Barcilon, A., Navon, I.M., Whitaker, J., Cacuci, D.G.. 1993: An Adjoint Sensitivity Study of Blocking in a
Two-Layer Isentropic Model. Monthly Weather Review: Vol. 121, No. 10, pp. 2833-2857.

Private components
------------------

N/A
