COAMPS
======

.. attention::

   ``COAMPS`` works with versions of DART *before* Manhattan (9.x.x) and has yet to be updated. If you are interested in
   using ``COAMPS`` with more recent versions of DART, contact DAReS staff to assess the feasibility of an update.
   Until that time, you should consider this documentation as out-of-date.


Overview
--------

| DART interface module for the Coupled Ocean / Atmosphere Mesoscale Prediction (COAMPS ®) model. The 16 public
  interfaces listed here are standardized for all DART compliant models. These interfaces allow DART to advance the
  model, get the model state and metadata describing this state, find state variables that are close to a given
  location, and do spatial interpolation for a variety of variables required in observational operators.
| The following model description is taken from the `COAMPS overview web
  page: <http://www.nrlmry.navy.mil/coamps-web/web/view>`__

   "The Coupled Ocean/Atmosphere Mesoscale Prediction System (COAMPS) has been developed by the Marine Meteorology
   Division (MMD) of the Naval Research Laboratory (NRL). The atmospheric components of COAMPS, described below, are
   used operationally by the U.S. Navy for short-term numerical weather prediction for various regions around the world.

   The atmospheric portion of COAMPS represents a complete three-dimensional data assimilation system comprised of data
   quality control, analysis, initialization, and forecast model components. Features include a globally relocatable
   grid, user-defined grid resolutions and dimensions, nested grids, an option for idealized or real-time simulations,
   and code that allows for portability between mainframes and workstations. The nonhydrostatic atmospheric model
   includes predictive equations for the momentum, the non-dimensional pressure perturbation, the potential temperature,
   the turbulent kinetic energy, and the mixing ratios of water vapor, clouds, rain, ice, grauple, and snow, and
   contains advanced parameterizations for boundary layer processes, precipitation, and radiation.

   The distributed version of the COAMPS code that can be downloaded from the web site has been designed to use the
   message-passing interface (MPI), OpenMP directives, and horizontal domain decomposition to achieve parallelism. The
   code is capable of executing efficiently across vector, parallel, or symmetric muti-processor (SMP) machines by
   simply changing run-time options."

Other modules used
------------------

::

   types_mod
   time_manager_mod
   threed_sphere/location_mod
   utilities_mod
   obs_kind_mod
   random_seq_mod
   netcdf
   typesizes
   coamps_grid_mod
   coamps_interp_mod
   coamps_restart_mod
   coamps_util_mod

Public interfaces
-----------------

======================= ======================
*use model_mod, only :* get_model_size
\                       get_state_meta_data
\                       model_interpolate
\                       get_model_time_step
\                       static_init_model
\                       nc_write_model_atts
\                       nc_write_model_vars
\                       pert_model_state
\                       get_close_maxdist_init
\                       get_close_obs_init
\                       get_close_obs
\                       ens_mean_for_model
\                       adv_1step
\                       end_model
\                       init_time
\                       init_conditions
======================= ======================

The last 4 interfaces are only required for low-order models where advancing the model can be done by a call to a
subroutine. The COAMPS model only advances by executing the coamps program. Thus the last 4 interfaces only appear as
stubs in this module.

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *model_size = get_model_size( )*
   ::

      integer :: get_model_size

.. container:: indent1

   Returns the length of the model state vector as an integer. This includes all nested domains.

   ============== =====================================
   ``model_size`` The length of the model state vector.
   ============== =====================================

| 

.. container:: routine

   *call get_state_meta_data (index_in, location, [, var_type] )*
   ::

      integer,             intent(in)  :: index_in
      type(location_type), intent(out) :: location
      integer, optional,   intent(out) ::  var_type 

.. container:: indent1

   Returns metadata about a given element, indexed by index_in, in the model state vector. The location defines where
   the state variable is located while the type of the variable (for instance temperature, or u wind component) is
   returned by var_type. The integer values used to indicate different variable types in var_type are themselves defined
   as public interfaces to model_mod if required.

   +--------------+------------------------------------------------------------------------------------------------------+
   | ``index_in`` | Index of state vector element about which information is requested.                                  |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``location`` | Returns location of indexed state variable. The location should use a location_mod that is           |
   |              | appropriate for the model domain. For realistic atmospheric models, for instance, a                  |
   |              | three-dimensional spherical location module that can represent height in a variety of ways is        |
   |              | provided.                                                                                            |
   +--------------+------------------------------------------------------------------------------------------------------+
   | *var_type*   | Returns the type of the indexed state variable as an optional argument.                              |
   +--------------+------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call model_interpolate(x, location, obs_kind, obs_val, istatus)*
   ::

      real(r8), dimension(:), intent(in)  :: x
      type(location_type),    intent(in)  :: location
      integer,                  intent(in)  ::  obs_kind 
      real(r8),               intent(out) :: obs_val
      integer,                intent(out) :: istatus

.. container:: indent1

   Given model state, returns the value of observation type interpolated to a given location by a method of the model's
   choosing. All observation kinds defined in obs_kind_mod are supported. In the case where the observational operator
   is not defined at the given location (e.g. the observation is below the model surface or outside the domain), obs_val
   is returned as -888888.0 and istatus = 1. Otherwise, istatus = 0. The interpolation is performed in the domain with
   the highest resolution containing the observation.

   ============ =================================================================
   ``x``        A model state vector.
   ``location`` Location to which to interpolate.
   ``obs_kind`` Integer indexing which type of observation is to be interpolated.
   ``obs_val``  The interpolated value from the model.
   ``istatus``  Integer flag indicating the result of the interpolation.
   ============ =================================================================

| 

.. container:: routine

   *var = get_model_time_step()*
   ::

      type(time_type) :: get_model_time_step

.. container:: indent1

   Returns the model base time step as a time_type. For now this is set to 1 minute.

   ======= ============================
   ``var`` Smallest time step of model.
   ======= ============================

| 

.. container:: routine

   *call static_init_model()*

.. container:: indent1

   Used for runtime initialization of the model. This is the first call made to the model by any DART compliant
   assimilation routine. It reads the model namelist parameters, initializes the pressure levels for the state vector,
   and generates the location data for each member of the state.

| 

.. container:: routine

   *ierr = nc_write_model_atts(ncFileId)*
   ::

      integer             ::  nc_write_model_atts
      integer, intent(in) ::  ncFileId 

.. container:: indent1

   Function to write model specific attributes to a netCDF file. At present, DART is using the NetCDF format to output
   diagnostic information. This is not a requirement, and models could choose to provide output in other formats. This
   function writes the metadata associated with the model to a NetCDF file opened to a file identified by ncFileID.

   ============ ==============================================
   ``ncFileId`` Integer file descriptor opened to NetCDF file.
   ``ierr``     Returned error code.
   ============ ==============================================

| 

.. container:: routine

   *ierr = nc_write_model_vars(ncFileID, statevec, copyindex, timeindex)*
   ::

      integer                            ::  nc_write_model_vars
      integer,                intent(in) ::  ncFileID 
      real(r8), dimension(:), intent(in) ::  statevec 
      integer,                intent(in) ::  copyindex
      integer,                intent(in) ::  timeindex 

.. container:: indent1

   Writes a copy of the state variables to a NetCDF file. Multiple copies of the state for a given time are supported,
   allowing, for instance, a single file to include multiple ensemble estimates of the state.

   ============= =========================================================
   ``ncFileID``  Integer file descriptor opened to NetCDF file.
   ``statevec``  State vector.
   ``copyindex`` Integer index to which copy is to be written.
   ``timeindex`` Integer index of which time in the file is being written.
   ``ierr``      Returned error code.
   ============= =========================================================

| 

.. container:: routine

   *call pert_model_state(state, pert_state, interf_provided)*
   ::

      real(r8), dimension(:),   intent(in)    ::  state 
      real(r8), dimension(:),   intent(out)   ::  pert_state 
      logical,                  intent(out)   ::  interf_provided

.. container:: indent1

   Given a model state, produces a perturbed model state. This is used to generate initial ensemble conditions perturbed
   around some control trajectory state when one is preparing to spin-up ensembles. In the COAMPS interface, this can be
   done three different ways:

   -  No perturbation
   -  Uniform perturbation - each element of the field has the same additive perturbation
   -  Individual perturbation - each element of the field has a different additive perturbation The perturbation
      magnitude and option are supplied out of the dynamic restart vector definition - this allows us to supply a
      variance appropriate for each type of variable at each level.

   =================== ===================================
   ``state``           State vector to be perturbed.
   ``pert_state``      Perturbed state vector is returned.
   ``interf_provided`` Returns .true. for this model.
   =================== ===================================

| 

.. container:: routine

   *call get_close_maxdist_init(gc, maxdist)*
   ::

      type(get_close_type), intent(inout) :: gc
      real(r8),             intent(in)    :: maxdist

.. container:: indent1

   Pass-through to the 3-D sphere locations module. See
   `get_close_maxdist_init() <../../location/threed-sphere/location_mod.html#get_close_maxdist_init>`__ for the
   documentation of this subroutine.

| 

.. container:: routine

   *call get_close_obs_init(gc, num, obs)*
   ::

      type(get_close_type), intent(inout) :: gc
      integer,              intent(in)    :: num
      type(location_type),  intent(in)    :: obs(num)

.. container:: indent1

   Pass-through to the 3-D sphere locations module. See
   `get_close_obs_init() <../../location/threed-sphere/location_mod.html#get_close_obs_init>`__ for the documentation of
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

   Pass-through to the 3-D sphere locations module. See
   `get_close_obs() <../../location/threed-sphere/location_mod.html#get_close_obs>`__ for the documentation of this
   subroutine.

| 

.. container:: routine

   *call ens_mean_for_model(ens_mean)*
   ::

      real(r8), dimension(:), intent(in)  :: ens_mean

.. container:: indent1

   A local copy is available here for use during other computations in the model_mod code.

   ============ ==========================
   ``ens_mean`` Ensemble mean state vector
   ============ ==========================

| 

.. container:: routine

   *call adv_1step(x, time)*
   ::

      real(r8), dimension(:),   intent(inout) ::  x 
      type(time_type),          intent(in)    ::  time 

.. container:: indent1

   This operation is not defined for the COAMPS model. This interface is only required if \`synchronous' model state
   advance is supported (the model is called directly as a Fortran90 subroutine from the assimilation programs). This is
   generally not the preferred method for large models and a stub for this interface is provided for the COAMPS model.

   +----------+----------------------------------------------------------------------------------------------------------+
   | ``x``    | State vector of length model_size.                                                                       |
   +----------+----------------------------------------------------------------------------------------------------------+
   | ``time`` | Gives time of the initial model state. Needed for models that have real time state requirements, for     |
   |          | instance the computation of radiational parameters. Note that DART provides a time_manager_mod module    |
   |          | that is used to support time computations throughout the facility.                                       |
   +----------+----------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call end_model( )*

.. container:: indent1

   Called when use of a model is completed to clean up storage, etc. A stub is provided for the COAMPS model.

| 

.. container:: routine

   *call init_time(i_time)*
   ::

      type(time_type),        intent(in)  ::  i_time 

.. container:: indent1

   Returns the time at which the model will start if no input initial conditions are to be used. This is frequently used
   to spin-up models from rest, but is not meaningfully supported for the COAMPS model.

| 

.. container:: routine

   *call init_conditions( x )*
   ::

      real(r8), dimension(:), intent(out) ::  x 

.. container:: indent1

   Returns default initial conditions for model; generally used for spinning up initial model states. For the COAMPS
   model just return 0's since initial state is always to be provided from input files.

   ===== ===================
   ``x`` Model state vector.
   ===== ===================

| 

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &model_nml
     cdtg = '2006072500',
     y_bound_skip = 3,
     x_bound_skip = 3,
     need_mean = .false.,
   /

| 

.. container::

   ========================== ================= ==========================================================================
   Item                       Type              Description
   ========================== ================= ==========================================================================
   cdtg                       character(len=10) Date/time group.
   x_bound_skip, y_bound_skip integer           Number of x and y boundary points to skip when perturbing the model state.
   need_mean                  logical           Does the forward operator computation need the ensemble mean?
   ========================== ================= ==========================================================================

| 

Files
-----

=========================== ===========================================================================
filename                    purpose
=========================== ===========================================================================
input.nml                   to read the model_mod namelist
preassim.nc                 the time-history of the model state before assimilation
analysis.nc                 the time-history of the model state after assimilation
dart_log.out [default name] the run-time diagnostic output
dart_log.nml [default name] the record of all the namelists actually USED - contains the default values
=========================== ===========================================================================

References
----------

The COAMPS registration web site is http://www.nrlmry.navy.mil/coamps-web/web/home and COAMPS is a registered trademark
of the Naval Research Laboratory.

Private components
------------------

N/A
