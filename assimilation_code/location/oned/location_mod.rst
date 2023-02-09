MODULE (1D) location_mod
========================

Overview
--------

The DART framework needs to be able to compute distances between locations, to pass location information to and from the
model interface code (in model_mod.f90), and to be able to read and write location information to files. DART isolates
all this location information into separate modules so that the main algorithms can operate with the same code
independent of whether the model uses latitude/longitude/height, one-d unit sphere coordinates, cylindrical coordinates,
etc. DART provides about half a dozen possible coordinate systems, and others can be added.

This locations module provides a representation of a physical location on a periodic 1D domain with location values
between 0 and 1. A type that abstracts the location is provided along with operators to set, get, read, write, and
compute distances between locations. This is a member of a class of similar location modules that provide the same
abstraction for different represenations of physical space.

All types of location modules define the same module name ``location_mod``. Therefore, the DART framework and any user
code should include a Fortran 90 ``use`` statement of ``location_mod``. The selection of which location module will be
compiled into the program is controlled by the LOCATION variable in ``quickbuild.sh``.

The model-specific ``model_mod.f90`` files need to define four ``get_close`` routines, but in most cases they can simply
put a ``use`` statement at the top which uses the routines in the locations module, and they do not have to provide any
additional code.

However, if the model interface code wants to intercept and alter the default behavior of the get_close routines, they
are able to. The correct usage of the ``get_close`` routines is as follows:

::


   call get_close_maxdist_init()  ! must be called before get_close_obs_init()
   call get_close_obs_init()
   ...
   call get_close_obs()           ! many, many times
   ...
   call get_close_obs_destroy()

Regardless of the fact that the names include the string 'obs', they are intended for use with any group of locations in
the system, frequently state vector items or observations, but any location is acceptable.

Namelist
--------

This version of the locations module does not have any namelist input.

Other modules used
------------------

::

   types_mod
   utilities_mod
   random_seq_mod

Public interfaces
-----------------

============================ ======================
``use location_mod, only :`` location_type
\                            get_close_type
\                            get_location
\                            set_location
\                            write_location
\                            read_location
\                            interactive_location
\                            set_location_missing
\                            query_location
\                            get_close_maxdist_init
\                            get_close_obs_init
\                            get_close_obs
\                            get_close_obs_destroy
\                            get_dist
\                            LocationDims
\                            LocationName
\                            LocationLName
\                            horiz_dist_only
\                            vert_is_undef
\                            vert_is_surface
\                            vert_is_pressure
\                            vert_is_level
\                            vert_is_height
\                            VERTISUNDEF
\                            VERTISSURFACE
\                            VERTISLEVEL
\                            VERTISPRESSURE
\                            VERTISHEIGHT
\                            operator(==)
\                            operator(/=)
============================ ======================

There is currently no namelist interface for the 1D location module.

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: type

   *type location_type*
   ::

         private
         real(r8) :: x
      end type location_type

.. container:: indent1

   Provides an abstract representation of physical location on a one-dimensional periodic domain.

   ========= =========================
   Component Description
   ========= =========================
   x         Location has range 0 to 1
   ========= =========================

| 

.. container:: type

   *type get_close_type*
   ::

         private
         integer  :: num
         real(r8) :: maxdist
      end type get_close_type

.. container:: indent1

   Provides a structure for doing efficient computation of close locations. Doesn't do anything in the 1D implementation
   except provide appropriate stubs.

   ========= ==============================================
   Component Description
   ========= ==============================================
   num       Number of locations in list
   maxdist   Threshhold distance. Anything closer is close.
   ========= ==============================================

| 

.. container:: routine

   *var = get_location(loc)*
   ::

      real(r8)                        :: get_location
      type(location_type), intent(in) :: loc

.. container:: indent1

   Extracts the real location value, range 0 to 1, from a location type.

   ================ =============================
   ``get_location`` The real value for a location
   ``loc``          A location derived type
   ================ =============================

| 

.. container:: routine

   *var = set_location(x)*
   ::

      type(location_type)   :: set_location
      real(r8), intent(in)  :: x

.. container:: indent1

   Returns a location type with the location x.

   ================ ====================================
   ``set_location`` A location derived type
   ``x``            Location value in the range 0. to 1.
   ================ ====================================

| 

.. container:: routine

   *call write_location(locfile, loc [, fform, charstring])*
   ::

      integer,               intent(in)       ::  locfile 
      type(location_type),   intent(in)       ::  loc 
      character(len=*), optional, intent(in)  ::  fform 
      character(len=*), optional, intent(out) ::  charstring 

.. container:: indent1

   Given an integer IO channel of an open file and a location, writes the location to this file. The *fform* argument
   controls whether write is "FORMATTED" or "UNFORMATTED" with default being formatted. If the final *charstring*
   argument is specified, the formatted location information is written to the character string only, and the
   ``locfile`` argument is ignored.

   +--------------+------------------------------------------------------------------------------------------------------+
   | ``locfile``  | the unit number of an open file.                                                                     |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``loc``      | location type to be written.                                                                         |
   +--------------+------------------------------------------------------------------------------------------------------+
   | *fform*      | Format specifier ("FORMATTED" or "UNFORMATTED"). Default is "FORMATTED" if not specified.            |
   +--------------+------------------------------------------------------------------------------------------------------+
   | *charstring* | Character buffer where formatted location string is written if present, and no output is written to  |
   |              | the file unit.                                                                                       |
   +--------------+------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *var = read_location(locfile [, fform])*
   ::

      type(location_type)                    :: read_location
      integer, intent(in)                    :: locfile
      character(len=*), optional, intent(in) :: fform

.. container:: indent1

   Reads a location_type from a file open on channel locfile using format *fform* (default is formatted).

   ================= ==============================================================================
   ``read_location`` Returned location type read from file
   ``locfile``       Integer channel opened to a file to be read
   *fform*           Optional format specifier ("FORMATTED" or "UNFORMATTED"). Default "FORMATTED".
   ================= ==============================================================================

| 

.. container:: routine

   *call interactive_location(location [, set_to_default])*
   ::

      type(location_type), intent(out) :: location
      logical, optional, intent(in)    :: set_to_default

.. container:: indent1

   Use standard input to define a location type. With set_to_default true get one with all elements set to 0.

   ================ ================================================
   ``location``     Location created from standard input
   *set_to_default* If true, sets all elements of location type to 0
   ================ ================================================

| 

.. container:: routine

   *var = query_location(loc [, attr])*
   ::

      real(r8)                               :: query_location
      type(location_type), intent(in)        :: loc
      character(len=*), optional, intent(in) :: attr

.. container:: indent1

   Returns the location value if attr = 'X' or if attr is not passed.

   ================== ===================
   ``query_location`` Returns value of x.
   ``loc``            A location type
   *attr*             Selects 'X'
   ================== ===================

| 

.. container:: routine

   *var = set_location_missing()*
   ::

      type(location_type) :: set_location_missing

.. container:: indent1

   Returns a location with location set to missing value from types_mod.

   ======================== ===============================
   ``set_location_missing`` A location set to missing value
   ======================== ===============================

| 

.. container:: routine

   *call get_close_maxdist_init(gc,maxdist , [maxdist_list])*
   ::

      type(get_close_type), intent(inout) :: gc
      real(r8), intent(in)                :: maxdist
      real(r8), intent(in), optional      :: maxdist_list(:)

.. container:: indent1

   Sets the threshhold distance. Anything closer than this is deemed to be close. This routine must be called first,
   before the other ``get_close`` routines. It allocates space so it is necessary to call ``get_close_obs_destroy`` when
   completely done with getting distances between locations.

   ============== =======================================================
   ``gc``         Data for efficiently finding close locations.
   ``maxdist``    Anything closer than this distance is a close location.
   *maxdist_list* Ignored for this location type.
   ============== =======================================================

| 

.. container:: routine

   *call get_close_obs_init(gc, num, obs)*
   ::

      type(get_close_type),             intent(inout) :: gc
      integer,                          intent(in)    :: num
      type(location_type), dimension(:) intent(in)    :: obs

.. container:: indent1

   Initialize storage for efficient identification of locations close to a given location. The oned implementation is
   minimal and just records the number of locations here. Must be called after ``get_close_maxdist_init``, and the list
   of locations here must be the same as the list of locations passed into ``get_close_obs()``. If the list changes,
   ``get_close_obs_destroy()`` must be called, and both the initialization routines must be called again. It allocates
   space so it is necessary to call ``get_close_obs_destroy`` when completely done with getting distances between
   locations.

   ======= =====================================================================================
   ``gc``  Structure that contains data to efficiently find locations close to a given location.
   ``num`` The number of locations in the list.
   ``obs`` The locations of each element in the list, not used in 1D implementation.
   ======= =====================================================================================

| 

.. container:: routine

   *call get_close_obs(gc, base_obs_loc, base_obs_kind, obs, obs_kind, num_close, close_ind, dist)*
   ::

      type(get_close_type),              intent(in)  :: gc
      type(location_type),               intent(in)  :: base_obs_loc
      integer,                           intent(in)  :: base_obs_kind
      type(location_type), dimension(:), intent(in)  :: obs
      integer, dimension(:),             intent(in)  :: obs_kind
      integer,                           intent(out) :: num_close
      integer, dimension(:),             intent(out) :: close_ind
      real(r8), dimension(:),            intent(out) :: dist

.. container:: indent1

   Given a single location and a list of other locations, returns the indices of all the locations close to the single
   one along with the number of these and the distances for the close ones. The list of locations passed in via the
   ``obs`` argument must be identical to the list of ``obs`` passed into the most recent call to
   ``get_close_obs_init()``. If the list of locations of interest changes ``get_close_obs_destroy()`` must be called and
   then the two initialization routines must be called before using ``get_close_obs()`` again.

   ================= ===================================================================================
   ``gc``            Structure to allow efficient identification of locations close to a given location.
   ``base_obs_loc``  Single given location.
   ``base_obs_kind`` Kind of the single location.
   ``obs``           List of locations from which close ones are to be found.
   ``obs_kind``      Kind associated with locations in obs list.
   ``num_close``     Number of locations close to the given location.
   ``close_ind``     Indices of those locations that are close.
   ``dist``          Distance between given location and the close ones identified in close_ind.
   ================= ===================================================================================

| 

.. container:: routine

   *call get_close_obs_destroy(gc)*
   ::

      type(get_close_type), intent(inout) :: gc

.. container:: indent1

   Releases memory associated with the ``gc`` derived type. Must be called whenever the list of locations changes, and
   then ``get_close_maxdist_init`` and ``get_close_obs_init`` must be called again with the new locations list.

   ====== =============================================
   ``gc`` Data for efficiently finding close locations.
   ====== =============================================

| 

.. container:: routine

   *var = get_dist(loc1, loc2, [, kind1, kind2])*
   ::

      real(r8)                        :: get_dist
      type(location_type), intent(in) :: loc1
      type(location_type), intent(in) :: loc2
      integer, optional,   intent(in) :: kind1
      integer, optional,   intent(in) :: kind2

.. container:: indent1

   Return the distance between 2 locations. Since this is a periodic domain, the shortest distance may wrap around.

   The kind arguments are not used by the default location code, but are available to any user-supplied distance
   routines which want to do specialized calculations based on the kinds associated with each of the two locations.

   ======== ====================================================
   ``loc1`` First of two locations to compute distance between.
   ``loc2`` Second of two locations to compute distance between.
   *kind1*  DART kind associated with location 1.
   *kind2*  DART kind associated with location 2.
   ``var``  distance between loc1 and loc2.
   ======== ====================================================

| 

.. container:: routine

   *var = vert_is_undef(loc)*
   ::

      logical                         :: vert_is_undef
      type(location_type), intent(in) :: loc

.. container:: indent1

   Always returns false; this locations module has no vertical coordinates. Provided only for compile-time compatibility
   with other location modules.

   ================= ======================
   ``vert_is_undef`` Always returns .FALSE.
   ``loc``           A location type
   ================= ======================

| 

.. container:: routine

   *var = vert_is_surface(loc)*
   ::

      logical                         :: vert_is_surface
      type(location_type), intent(in) :: loc

.. container:: indent1

   Always returns false; this locations module has no vertical coordinates. Provided only for compile-time compatibility
   with other location modules.

   =================== ======================
   ``vert_is_surface`` Always returns .FALSE.
   ``loc``             A location type
   =================== ======================

| 

.. container:: routine

   *var = vert_is_pressure(loc)*
   ::

      logical                         :: vert_is_pressure
      type(location_type), intent(in) :: loc

.. container:: indent1

   Always returns false; this locations module has no vertical coordinates. Provided only for compile-time compatibility
   with other location modules.

   ==================== ======================
   ``vert_is_pressure`` Always returns .FALSE.
   ``loc``              A location type
   ==================== ======================

| 

.. container:: routine

   *var = vert_is_level(loc)*
   ::

      logical                         :: vert_is_level
      type(location_type), intent(in) :: loc

.. container:: indent1

   Always returns false; this locations module has no vertical coordinates. Provided only for compile-time compatibility
   with other location modules.

   ================= ======================
   ``vert_is_level`` Always returns .FALSE.
   ``loc``           A location type
   ================= ======================

| 

.. container:: routine

   *var = vert_is_height(loc)*
   ::

      logical                         :: vert_is_height
      type(location_type), intent(in) :: loc

.. container:: indent1

   Always returns false; this locations module has no vertical coordinates. Provided only for compile-time compatibility
   with other location modules.

   ================== ======================
   ``vert_is_height`` Always returns .FALSE.
   ``loc``            A location type
   ================== ======================

| 

.. container:: routine

   *var = has_vertical_localization()*
   ::

      logical :: has_vertical_localization

.. container:: indent1

   Always returns false; this locations module has no vertical coordinates. Provided only for compile-time compatibility
   with other location modules.

   See note in threed_sphere locations module about the function name.

| 

.. container:: routine

   *loc1 == loc2*
   ::

      type(location_type), intent(in) :: loc1, loc2

.. container:: indent1

   Returns true if the two location types have identical values, else false.

| 

.. container:: routine

   *loc1 /= loc2*
   ::

      type(location_type), intent(in) :: loc1, loc2

.. container:: indent1

   Returns true if the two location types do NOT have identical values, else false.

| 

.. container:: routine

   ::

      integer, parameter :: VERTISUNDEF    = -2
      integer, parameter :: VERTISSURFACE  = -1
      integer, parameter :: VERTISLEVEL    =  1
      integer, parameter :: VERTISPRESSURE =  2
      integer, parameter :: VERTISHEIGHT   =  3

.. container:: indent1

   This locations module has no vertical coordinate, but for compatibility with other location modules, these are
   defined.

| 

.. container:: routine

   ::

      integer, parameter :: LocationDims = 1

.. container:: indent1

   This is a **constant**. Contains the number of real values in a location type. Useful for output routines that must
   deal transparently with many different location modules.

| 

.. container:: routine

   ::

      character(len=129), parameter :: LocationName = "loc1d"

.. container:: indent1

   This is a **constant**. A parameter to identify this location module in output metadata.

| 

.. container:: routine

   ::

      character(len=129), parameter :: LocationLName = "location on unit circle"

.. container:: indent1

   This is a **constant**. A parameter to identify this location module in output long name metadata.

| 

Files
-----

None.

References
----------

#. none

Private components
------------------

N/A
