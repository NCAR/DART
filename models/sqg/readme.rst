SQG
===


.. attention::

   ``sqg`` works with versions of DART *before* Manhattan (9.x.x) and has yet to be updated. If you are interested in
   using ``sqg`` with more recent versions of DART, contact DAReS staff to assess the feasibility of an update.
   Until that time, you should consider this documentation as out-of-date.


Overview
--------

This is a uniform PV two-surface QG+1 spectral model contributed by Rahul Majahan.

The underlying model is described in: Hakim, Gregory J., 2000: Role of Nonmodal Growth and Nonlinearity in Cyclogenesis
Initial-Value Problems. J. Atmos. Sci., 57, 2951-2967. doi: 10.1175/1520-0469(2000)057<2951:RONGAN>2.0.CO;2

Other modules used
------------------

::

   types_mod
   time_manager_mod
   threed_sphere/location_mod
   utilities_mod

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

Optional namelist interface ``&model_nml`` may be read from file ``input.nml``.

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *model_size = get_model_size( )*
   ::

      integer :: get_model_size

.. container:: indent1

   Returns the length of the model state vector.

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

   Returns metadata about a given element, indexed by index_in, in the model state vector. The location defines where
   the state variable is located.

   ============ ==================================================================================
   ``index_in`` Index of state vector element about which information is requested.
   ``location`` The location of state variable element.
   *var_type*   Returns the type (always 1) of the indexed state variable as an optional argument.
   ============ ==================================================================================

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

   Given model state, returns the value interpolated to a given location.

   ============ ===============================================
   ``x``        A model state vector.
   ``location`` Location to which to interpolate.
   ``itype``    Not used.
   ``obs_val``  The interpolated value from the model.
   ``istatus``  Quality control information, always returned 0.
   ============ ===============================================

| 

.. container:: routine

   *var = get_model_time_step()*
   ::

      type(time_type) :: get_model_time_step

.. container:: indent1

   Returns the time step (forecast length) of the model;

   ======= ============================
   ``var`` Smallest time step of model.
   ======= ============================

| 

.. container:: routine

   *call static_init_model()*

.. container:: indent1

   Used for runtime initialization of model; reads namelist, initializes model parameters, etc. This is the first call
   made to the model by any DART-compliant assimilation routine.

| 

.. container:: routine

   *call end_model()*

.. container:: indent1

   A stub.

| 

.. container:: routine

   *call init_time(time)*
   ::

      type(time_type), intent(out) :: time

.. container:: indent1

   Returns the time at which the model will start if no input initial conditions are to be used. This is used to spin-up
   the model from rest.

   ======== ===================
   ``time`` Initial model time.
   ======== ===================

| 

.. container:: routine

   *call init_conditions(x)*
   ::

      real(r8), dimension(:), intent(out) :: x

.. container:: indent1

   Returns default initial conditions for the model; generally used for spinning up initial model states.

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

   Function to write model specific attributes to a netCDF file. At present, DART is using the NetCDF format to output
   diagnostic information. This is not a requirement, and models could choose to provide output in other formats. This
   function writes the metadata associated with the model to a NetCDF file opened to a file identified by ncFileID.

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

   Writes a copy of the state variables to a netCDF file. Multiple copies of the state for a given time are supported,
   allowing, for instance, a single file to include multiple ensemble estimates of the state.

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

   Given a model state, produces a perturbed model state.

   =================== =============================================
   ``state``           State vector to be perturbed.
   ``pert_state``      Perturbed state vector: NOT returned.
   ``interf_provided`` Returned false; interface is not implemented.
   =================== =============================================

| 

.. container:: routine

   *call get_close_maxdist_init(gc, maxdist)*
   ::

      type(get_close_type), intent(inout) :: gc
      real(r8),             intent(in)    :: maxdist

.. container:: indent1

   Pass-through to the 3D Sphere locations module. See
   `get_close_maxdist_init() <../../location/threed_sphere/location_mod.html#get_close_maxdist_init>`__ for the
   documentation of this subroutine.

| 

.. container:: routine

   *call get_close_obs_init(gc, num, obs)*
   ::

      type(get_close_type), intent(inout) :: gc
      integer,              intent(in)    :: num
      type(location_type),  intent(in)    :: obs(num)

.. container:: indent1

   Pass-through to the 3D Sphere locations module. See
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

   Pass-through to the 3D Sphere locations module. See
   `get_close_obs() <../../location/threed_sphere/location_mod.html#get_close_obs>`__ for the documentation of this
   subroutine.

| 

.. container:: routine

   *call ens_mean_for_model(ens_mean)*
   ::

      real(r8), dimension(:), intent(in) :: ens_mean

.. container:: indent1

   A NULL INTERFACE in this model.

   ============ ==========================================
   ``ens_mean`` State vector containing the ensemble mean.
   ============ ==========================================

| 

Namelist
--------

We adhere to the F90 standard of starting a namelist with an ampersand '&' and terminating with a slash '/' for all our
namelist input.

::

   &model_nml 
     output_state_vector = .false.
     channel_center = 45.0
     channel_width = 40.0
     assimilation_period_days = 0
     assimilation_period_seconds = 21600
     debug = .false.
   /

.. container:: indent1

   This namelist is read in a file called ``input.nml``

   +-----------------------------+----------+---------------------------------------------------------------------------+
   | Contents                    | Type     | Description                                                               |
   +=============================+==========+===========================================================================+
   | output_state_vector         | logical  | If .true. write state vector as a 1D array to the diagnostic output file. |
   |                             |          | If .false. break state vector up into fields before writing to the        |
   |                             |          | outputfile.                                                               |
   +-----------------------------+----------+---------------------------------------------------------------------------+
   | channel_center              | real(r8) | Channel center                                                            |
   +-----------------------------+----------+---------------------------------------------------------------------------+
   | channel_width               | real(r8) | Channel width                                                             |
   +-----------------------------+----------+---------------------------------------------------------------------------+
   | assimilation_period_days    | integer  | Number of days for timestep                                               |
   +-----------------------------+----------+---------------------------------------------------------------------------+
   | assimilation_period_seconds | integer  | Number of seconds for timestep                                            |
   +-----------------------------+----------+---------------------------------------------------------------------------+
   | debug                       | logical  | Set to .true. for more output                                             |
   +-----------------------------+----------+---------------------------------------------------------------------------+

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

| The underlying model is described in:
| Hakim, Gregory J., 2000: Role of Nonmodal Growth and Nonlinearity in Cyclogenesis Initial-Value Problems. J. Atmos.
  Sci., 57, 2951-2967. doi: 10.1175/1520-0469(2000)057<2951:RONGAN>2.0.CO;2

Private components
------------------

N/A
