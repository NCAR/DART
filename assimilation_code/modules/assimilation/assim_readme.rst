MODULE assim_model_mod
======================

Overview
--------

This module acts as an intermediary between DART compliant models and the filter. At one time the assim_model_type,
which combines a state vector and a time_type, was envisioned as being fundamental to how DART views model states. This
paradigm is gradually being abandoned so that model state vectors and times are handled as separate data types. It is
important to call static_init_assim_model before using routines in assim_model_mod. Interfaces to work with model time
stepping, restart files, and computations about the locations of model state variables and the distance between
observations and state variables. Many of the interfaces are passed through nearly directly to the model_mod.

Notes
~~~~~

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

Namelist
--------

This module does not have a namelist.

Other modules used
------------------

::

   types_mod
   location_mod (model dependent choice)
   time_manager_mod
   utilities_mod
   model_mod
   netcdf
   typeSizes (part of netcdf)

Public interfaces
-----------------

============================= =============================
*use assim_model_mod, only :* 
\                             adv_1step
\                             aoutput_diagnostics
\                             aread_state_restart
\                             assim_model_type
\                             awrite_state_restart
\                             close_restart
\                             copy_assim_model
\                             end_assim_model
\                             ens_mean_for_model
\                             finalize_diag_output
\                             get_close_maxdist_init
\                             get_close_obs
\                             get_close_obs_init
\                             get_closest_state_time_to
\                             get_diag_input_copy_meta_data
\                             get_initial_condition
\                             get_model_size
\                             get_model_state_vector
\                             get_model_time
\                             get_model_time_step
\                             get_state_meta_data
\                             init_assim_model
\                             init_diag_input
\                             init_diag_output
\                             input_diagnostics
\                             interpolate
\                             nc_append_time
\                             nc_get_tindex
\                             nc_write_calendar_atts
\                             netcdf_file_type
\                             open_restart_read
\                             open_restart_write
\                             output_diagnostics
\                             pert_model_state
\                             read_state_restart
\                             set_model_state_vector
\                             set_model_time
\                             static_init_assim_model
\                             write_state_restart
============================= =============================

| 

.. container:: routine

   ::

      type assim_model_type
         private
         real(r8), pointer   :: state_vector(:) 
         type(time_type)     :: time
         integer             :: model_size
         integer             :: copyID
      end type assim_model_type

.. container:: indent1

   This type is used to represent both the state and time of a state from a model.

   ============ ===========================================================
   Component    Description
   ============ ===========================================================
   state_vector A one dimensional representation of the model state vector.
   time         The time of the model state.
   model_s      Size of the model state vector.
   copyID       Not used in present implementation.
   ============ ===========================================================

| 

.. container:: routine

   ::

      type netcdf_file_type
         integer             :: ncid
         integer             :: Ntimes
         integer             :: NtimesMAX
         real(r8), pointer   :: rtimes(:)
         type(time_type), pointer :: times(:)
         character(len = 80)      :: fname
      end type netcdf_file_type

.. container:: indent1

   Basically, we want to keep a local mirror of the unlimited dimension coordinate variable (i.e. time) because
   dynamically querying it causes unacceptable performance degradation over "long" integrations.

   ========= ===========================
   Component Description
   ========= ===========================
   ncid      The netcdf file unit id.
   Ntimes    The current working length.
   NtimesMAX Allocated length.
   rtimes    Times as real (r8).
   times     Times as time_types.
   fname     Netcdf file name.
   ========= ===========================

| 

.. container:: routine

   *call static_init_assim_model()*

.. container:: indent1

   Initializes the assim_model class. Must be called before any other assim_model_mod interfaces are used. Also calls
   the static initialization for the underlying model. There are no arguments.

| 

.. container:: routine

   *ncFileID = init_diag_output(FileName, global_meta_data, copies_of_field_per_time, meta_data_per_copy [, lagID])*
   ::

      type(netcdf_file_type)          :: init_diag_output 
      character (len = *), intent(in) :: FileName 
      character (len = *), intent(in) :: global_meta_data 
      integer, intent(in)             :: copies_of_field_per_time 
      character (len = *), intent(in) :: meta_data_per_copy(copies_of_field_per_time) 
      integer, optional, intent(in)   :: lagID 

.. container:: indent1

   Initializes a netCDF file for output of state space diagnostics. A handle to the channel on which the file is opened
   is returned.

   +------------------------------+--------------------------------------------------------------------------------------+
   | ``ncFileID``                 | Identifier for the netcdf file is returned. This is not an integer unit number, but  |
   |                              | a derived type containing additional information about the opened file.              |
   +------------------------------+--------------------------------------------------------------------------------------+
   | ``FileName``                 | Name of file to open.                                                                |
   +------------------------------+--------------------------------------------------------------------------------------+
   | ``global_meta_data``         | Global metadata that describes the contents of this file.                            |
   +------------------------------+--------------------------------------------------------------------------------------+
   | ``copies_of_field_per_time`` | Number of copies of data to be written at each time. For instance, these could be    |
   |                              | the prior ensemble members, prior ensemble mean, prior ensemble spread, posterior    |
   |                              | ensemble members, posterior spread and mean, etc..                                   |
   +------------------------------+--------------------------------------------------------------------------------------+
   | ``meta_data_per_copy``       | Metadata describing each of the copies.                                              |
   +------------------------------+--------------------------------------------------------------------------------------+
   | *lagID*                      | If using the smoother, which lag number this output is for.                          |
   +------------------------------+--------------------------------------------------------------------------------------+

| 

.. container:: routine

   *var = get_model_size()*
   ::

      integer :: get_model_size 

.. container:: indent1

   Returns the size of the model state vector. This is a direct pass through to the model_mod.

| 

.. container:: routine

   *var = get_closest_state_time_to(model_time, time)*
   ::

      type(time_type)              ::  get_closest_state_time_to 
      type(time_type), intent(in)  ::  model_time 
      type(time_type), intent(in)  ::  time

.. container:: indent1

   Returns the closest time that a model is capable of advancing a given state to a specified time. For instance, what
   is the closest time to 12GMT 01 January, 2004 that a model state at 00GMT 01 January, 2004 can be advanced? If the
   model time is past the time, the model time is returned (new feature in releases after Hawaii).

   ============== ================================================================
   ``var``        The closest time to which the model can be advanced is returned.
   ``model_time`` The time of a model state vector.
   ``time``       A time that one would like to get close to with the model.
   ============== ================================================================

| 

.. container:: routine

   *call get_state_meta_data()*

.. container:: indent1

   Pass through to model_mod. See model_mod documentation for arguments and description.

| 

.. container:: routine

   *var = get_model_time(assim_model)*
   ::

      type(time_type)                    :: get_model_time
      type(assim_model_type), intent(in) :: assim_model

.. container:: indent1

   Returns time from an assim_model type.

   =============== ===========================================
   ``var``         Returned time from assim_model
   ``assim_model`` Assim_model type from which to extract time
   =============== ===========================================

| 

.. container:: routine

   *var = get_model_state_vector(assim_model)*
   ::

      real(r8)                           :: get_model_state_vector(model_size)
      type(assim_model_type), intent(in) :: assim_model

.. container:: indent1

   Returns the state vector component from an assim_model_type.

   =============== ======================
   ``var``         Returned state vector
   ``assim_model`` Input assim_model_type
   =============== ======================

| 

.. container:: routine

   *call copy_assim_model(model_out, model_in)*
   ::

      type(assim_model_type), intent(out) :: model_out
      type(assim_model_type), intent(in)  :: model_in

.. container:: indent1

   Copies one assim_model_type to another.

   ============= ==================
   ``model_out`` Copy.
   ``model_in``  Data to be copied.
   ============= ==================

| 

.. container:: routine

   *call interpolate(x, location, loctype, obs_vals, istatus)*
   ::

      real(r8),            intent(in)  :: x(:)
      type(location_type), intent(in)  :: location
      integer,             intent(in)  :: loctype
      real(r8),            intent(out) :: obs_vals
      integer,             intent(out) :: istatus

.. container:: indent1

   Interpolates a given model state variable type to a location given the model state vector. Nearly direct call to
   model_interpolate in model_mod. See model_mod for the error return values in istatus.

   ============ ==================================================
   ``x``        Model state vector.
   ``location`` Location to which to interpolate.
   ``loctype``  Type of variable to interpolate.
   ``obs_vals`` Returned interpolated value.
   ``istatus``  Returned as 0 if all is well, else various errors.
   ============ ==================================================

| 

.. container:: routine

   *call set_model_time(assim_model, time)*
   ::

      type(assim_model_type), intent(inout) :: assim_model
      type(time_type), intent(in)           :: time

.. container:: indent1

   Sets the time in an assim_model_type.

   =============== ======================================
   ``assim_model`` Set the time in this assim_model_type.
   ``time``        Set to this time
   =============== ======================================

| 

.. container:: routine

   *call set_model_state_vector(assim_model, state)*
   ::

      type(assim_model_type), intent(inout) :: assim_model
      real(r8), intent(in)                  :: state(:)

.. container:: indent1

   Set the state in an assim_model_type.

   =============== ==============================================
   ``assim_model`` Set the state vector in this assim_model_type.
   ``state``       The state vector to be inserted.
   =============== ==============================================

| 

.. container:: routine

   *call write_state_restart(assim_model, funit [, target_time])*
   ::

      type(assim_model_type),    intent(in) :: assim_model
      integer,                   intent(in) :: funit
      type(time_type), optional, intent(in) :: target_time

.. container:: indent1

   Writes a restart from an assim_model_type with an optional target_time.

   =============== ==================================================================
   ``assim_model`` Write a restart from this assim_model_type.
   ``funit``       Integer file unit id open for output of restart files.
   *target_time*   If present, put this target time at the front of the restart file.
   =============== ==================================================================

| 

.. container:: routine

   *call read_state_restart(assim_model, funit [, target_time])*
   ::

      type(assim_model_type),    intent(out) :: assim_model
      integer,                   intent(in)  :: funit
      type(time_type), optional, intent(out) :: target_time

.. container:: indent1

   Read a state restart file into assim_model_type. Optionally read a prepended target time.

   =============== ====================================================================
   ``assim_model`` Read the time and state vector from restart into this.
   ``funit``       File id that has been opened for reading restart files.
   *target_time*   If present, read a target time from the front of the file into this.
   =============== ====================================================================

| 

.. container:: routine

   *call output_diagnostics(ndFileID, state [, copy_index])*
   ::

      type(netcdf_file_type), intent(inout) :: ndFileID
      type(assim_model_type), intent(in)    :: state
      integer, optional,      intent(in)    :: copy_index

.. container:: indent1

   Writes one copy of the state time and vector to a netCDF file.

   ============ ===================================
   ``ndFileID`` An identifier for a netCDF file
   ``state``    State vector and time
   *copy_index* Which copy of state is to be output
   ============ ===================================

| 

.. container:: routine

   *call end_assim_model()*

.. container:: indent1

   Called to clean-up at end of assim_model use. For now just passes through to model_mod.

| 

.. container:: routine

   *call input_diagnostics(file_id, state, copy_index)*
   ::

      integer,                intent(in)    :: file_id
      type(assim_model_type), intent(inout) :: state
      integer,                intent(out)   :: copy_index

.. container:: indent1

   Used to read in a particular copy of the state vector from an open state diagnostics file.

   ============== ======================================================================
   ``file_id``    Integer descriptor (channel number) for a diagnostics file being read.
   ``state``      Assim_model_type to read in data.
   ``copy_index`` Which copy of state to be read.
   ============== ======================================================================

| 

.. container:: routine

   *var = init_diag_input(file_name, global_meta_data, model_size, copies_of_field_per_time)*
   ::

      integer                       :: init_diag_input
      character(len=*), intent(in)  :: file_name
      character(len=*), intent(out) :: global_meta_data
      integer,          intent(out) :: model_size
      integer,          intent(out) :: copies_of_field_per_time

.. container:: indent1

   Opens a state diagnostic file and reads the global meta data, model size, and number of data copies.

   ============================ ==================================================
   ``var``                      Returns the unit number on which the file is open.
   ``file_name``                File name of state diagnostic file.
   ``global_meta_data``         Global metadata string from file.
   ``model_size``               Size of model.
   ``copies_of_field_per_time`` Number of copies of the state vector at each time.
   ============================ ==================================================

| 

.. container:: routine

   *call init_assim_model(state)*
   ::

      type(assim_model_type), intent(inout) :: state

.. container:: indent1

   Creates storage for an assim_model_type.

   ========= ===============================================
   ``state`` An assim_model_type that needs storage created.
   ========= ===============================================

| 

.. container:: routine

   *call get_diag_input_copy_meta_data(file_id, model_size_out, num_copies, location, meta_data_per_copy)*
   ::

      integer,             intent(in)  :: file_id
      integer,             intent(in)  :: model_size_out
      integer,             intent(in)  :: num_copies
      type(location_type), intent(out) :: location(model_size_out)
      character(len = *)               :: meta_data_per_copy(num_copies)

.. container:: indent1

   Reads meta-data describing state vectors in a state diagnostics file. Given the file, the model_size, and the number
   of copies, returns the locations of each state variable and the text description of each copy.

   ====================== =========================================================
   ``file_id``            Integer channel open to state diagostic file being read
   ``Model_size_out``     model size
   ``num_copies``         Number of copies of state in file
   ``location``           Returned locations for state vector
   ``meta_data_per_copy`` Meta data describing what is in each copy of state vector
   ====================== =========================================================

| 

.. container:: routine

   *var = finalize_diag_output(ncFileID)*
   ::

      integer                               :: finalize_diag_output
      type(netcdf_file_type), intent(inout) :: ncFileID

.. container:: indent1

   Used to complete writing on and open netcdf file. An error return is provided for passing to the netcdf error
   handling routines.

   ============ ===============================
   ``var``      Returns an error value.
   ``ncFileID`` Netcdf file id of an open file.
   ============ ===============================

| 

.. container:: routine

   *call aread_state_restart(model_time, model_state, funit [, target_time])*
   ::

      type(time_type),           intent(out) :: model_time
      real(r8),                  intent(out) :: model_state(:)
      integer,                   intent(in)  :: funit
      type(time_type), optional, intent(out) :: target_time

.. container:: indent1

   Reads a model time and state, and optionally a prepended target time, from a state restart file.

   =============== =================================================================
   ``model_time``  Returned time of model state
   ``model_state`` Returned model state.
   ``funit``       Channel open for reading a state restart file.
   *target_time*   If present, this time is read from the front of the restart file.
   =============== =================================================================

| 

.. container:: routine

   *call aoutput_diagnostics(ncFileID, model_time, model_state [, copy_index])*
   ::

      type(netcdf_file_type), intent(inout) :: ncFileID
      type(time_type),        intent(in)    :: model_time
      real(r8),               intent(in)    :: model_state(:)
      integer, optional,      intent(in)    :: copy_index

.. container:: indent1

   Write a state vector to a state diagnostics netcdf file.

   =============== ==============================================================
   ``ncFileID``    Unit for a state vector netcdf file open for output.
   ``model_time``  The time of the state to be output
   ``model_state`` A model state vector to be output.
   *copy_index*    Which copy of state vector is to be written, default is copy 1
   =============== ==============================================================

| 

.. container:: routine

   *call awrite_state_restart(model_time, model_state, funit [, target_time])*
   ::

      type(time_type),           intent(in) :: model_time
      real(r8),                  intent(in) :: model_state(:)
      integer,                   intent(in) :: funit
      type(time_type), optional, intent(in) :: target_time

.. container:: indent1

   Writes a model time and state vector to a restart file and optionally prepends a target time.

   =============== ========================================================
   ``model_time``  Time of model state.
   ``model_state`` Model state vector.
   ``funit``       Channel of file open for restart output.
   *target_time*   If present, time to be prepended to state time / vector.
   =============== ========================================================

| 

.. container:: routine

   *call pert_model_state()*

.. container:: indent1

   Passes through to pert_model_state in model_mod. See model_mod documentation for arguments and details.

| 

.. container:: routine

   *var = nc_append_time(ncFileID, time)*
   ::

      integer                               :: nc_append_time
      type(netcdf_file_type), intent(inout) :: ncFileID
      type(time_type),        intent(in)    :: time

.. container:: indent1

   Appends the time to the time coordinate variable of the netcdf file. The new length of the time variable is returned.
   Requires that time is a coordinate variable AND it is the unlimited dimension.

   ============ ======================================
   ``var``      Returns new length of time variable.
   ``ncFileID`` Points to open netcdf file.
   ``time``     The next time to be added to the file.
   ============ ======================================

| 

.. container:: routine

   *var = nc_write_calendar_atts(ncFileID, TimeVarID)*
   ::

      integer                            :: nc_write_calendar_atts
      type(netcdf_file_type), intent(in) :: ncFileID
      integer,                intent(in) :: TimeVarID

.. container:: indent1

   Sets up the metadata for the appropriate calendar being used in the time manager an writes it to a netcdf file.

   ============= ===================================================
   ``var``       Returns a netcdf error code.
   ``ncFileID``  Netcdf file id pointing to a file open for writing.
   ``TimeVarID`` The index of the time variable in the netcdf file.
   ============= ===================================================

| 

.. container:: routine

   *var = nc_get_tindex(ncFileID, statetime)*
   ::

      integer                               :: nc_get_tindex
      type(netcdf_file_type), intent(inout) :: ncFileID
      type(time_type),        intent(in)    :: statetime

.. container:: indent1

   Returns the index of a time from the time variable in a netcdf file. This function has been replaced with more
   efficient approaches and may be deleted from future releases.

   ============= =========================================
   ``var``       The index of the time in the netcdf file.
   ``ncFileID``  File id for an open netcdf file.
   ``statetime`` The time to be found in the netcdf file.
   ============= =========================================

| 

.. container:: routine

   *var = get_model_time_step()*
   ::

      type(time_type) :: get_model_time_step

.. container:: indent1

   This passes through to model_mod. See model_mod documentation for arguments and details.

   ======= ===========================
   ``var`` Returns time step of model.
   ======= ===========================

| 

.. container:: routine

   *var = open_restart_read(file_name)*
   ::

      integer                      :: open_restart_read
      character(len=*), intent(in) :: file_name

.. container:: indent1

   Opens a restart file for readig.

   ============= ============================================
   ``var``       Returns a file descriptor (channel number).
   ``file_name`` Name of restart file to be open for reading.
   ============= ============================================

| 

.. container:: routine

   *var = open_restart_write(file_name)*
   ::

      integer                      :: open_restart_write
      character(len=*), intent(in) :: file_name

.. container:: indent1

   Open a restart file for writing.

   ============= =======================================================
   ``var``       Returns a file descriptor (channel) for a restart file.
   ``file_name`` File name of restart file to be opened.
   ============= =======================================================

| 

.. container:: routine

   *call close_restart(file_unit)*
   ::

      integer, intent(in) :: file_unit

.. container:: indent1

   Closes a restart file.

   ============= ======================================================
   ``file_unit`` File descriptor (channel number) of open restart file.
   ============= ======================================================

| 

.. container:: routine

   *call adv_1step()*

.. container:: indent1

   Advances a model by one step. Pass through to model_mod. See model_mod documentation for arguments and details.

| 

.. container:: routine

   *call get_initial_condition(time, x)*
   ::

      type(time_type), intent(out) :: time
      real(r8),        intent(out) :: x

.. container:: indent1

   Obtains an initial condition from models that support this option.

   ======== =================================
   ``time`` the valid time of the model state
   ``x``    the initial model state
   ======== =================================

| 

.. container:: routine

   *call ens_mean_for_model(ens_mean)*
   ::

      type(r8), intent(in) :: ens_mean(:)

.. container:: indent1

   An array of length model_size containing the ensemble means. This is a direct pass through to the model_mod.

   ============ ==================================================================================
   ``ens_mean`` Array of length model_size containing the mean for each entry in the state vector.
   ============ ==================================================================================

| 

.. container:: routine

   *call get_close_maxdist_init(gc, maxdist)*
   ::

      type(get_close_type), intent(inout) :: gc
      type(r8), intent(in)                :: maxdist

.. container:: indent1

   Sets the threshold distance. Anything closer than this is deemed to be close. This is a direct pass through to the
   model_mod, which in turn can pass through to the location_mod.

   =========== =======================================================
   ``gc``      Data for efficiently finding close locations.
   ``maxdist`` Anything closer than this distance is a close location.
   =========== =======================================================

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
      real(r8),  optional,  intent(out) :: dist(:)

.. container:: indent1

   Given a single location and a list of other locations, returns the indices of all the locations close to the single
   one along with the number of these and the distances for the close ones. The observation kinds are passed in to allow
   more sophisticated distance computations to be done if needed. This is a direct pass through to the model_mod, which
   in turn can pass through to the location_mod.

   ================= ===========================================================================
   ``gc``            Data for efficiently finding close locations.
   ``base_obs_loc``  Single given location.
   ``base_obs_kind`` Kind of the single location.
   ``obs``           List of observations from which close ones are to be found.
   ``obs_kind``      Kind associated with observations in obs list.
   ``num_close``     Number of observations close to the given location.
   ``close_ind``     Indices of those locations that are close.
   *dist*            Distance between given location and the close ones identified in close_ind.
   ================= ===========================================================================

| 

.. container:: routine

   *call get_close_obs_init(gc, num, obs)*
   ::

      type(get_close_type), intent(inout) :: gc
      integer,              intent(in)    :: num
      type(location_type),  intent(in)    :: obs(:)

.. container:: indent1

   Initialize storage for efficient identification of locations close to a given location. Allocates storage for keeping
   track of which 'box' each observation in the list is in. This is a direct pass through to the model_mod, which in
   turn can pass through to the location_mod.

   ======= ========================================================================
   ``gc``  Data for efficiently finding close locations.
   ``num`` The number of locations in the list.
   ``obs`` The location of each element in the list, not used in 1D implementation.
   ======= ========================================================================

| 

Files
-----

============== =============================================
filename       purpose/comment
============== =============================================
filter_restart specified in &filter_nml:restart_in_filename
filter_restart specified in &filter_nml:restart_out_filename
input.nml      to read namelists
============== =============================================

References
----------

-  none

Private components
------------------

N/A
