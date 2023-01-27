MODULE location_mod (channel)
=============================

Overview
--------

THIS HAS NOT BEEN UPDATED YET - ONLY COPIED FROM 3D SPHERE VERSION

THIS HAS NOT BEEN UPDATED YET - ONLY COPIED FROM 3D SPHERE VERSION

THIS HAS NOT BEEN UPDATED YET - ONLY COPIED FROM 3D SPHERE VERSION

The DART framework needs to be able to compute distances between locations, to pass location information to and from the
model interface code (model_mod.f90), and to be able to read and write location information to files. DART isolates all
this location information into separate modules so that the main algorithms can operate with the same code independent
of whether the model uses latitude/longitude/height, 1D unit sphere coordinates, cylindrical coordinates, etc. DART
provides about half a dozen possible coordinate systems, and others can be added. The most common one for geophysical
models is this one: threed_sphere.

This location module provides a representation of a physical location on a 3-D spherical shell, using latitude and
longitude plus a vertical component with choices of vertical coordinate type such as pressure or height in meters. A
type that abstracts the location is provided along with operators to set, get, read, write, and compute distances
between locations. This is a member of a class of similar location modules that provide the same abstraction for
different represenations of physical space.

Location-independent code
^^^^^^^^^^^^^^^^^^^^^^^^^

All types of location modules define the same module name ``location_mod``. Therefore, the DART framework and any user
code should include a Fortran 90 ``use`` statement of ``location_mod``. The selection of which location module will be
compiled into the program is controlled by the LOCATION variable in ``quickbuild.sh``.

All types of location modules define the same Fortran 90 derived type ``location_type``. Programs that need to pass
location information to subroutines but do not need to interpret the contents can declare, receive, and pass this
derived type around in their code independent of which location module is specified at compile time. Model and
location-independent utilities should be written in this way. However, as soon as the contents of the location type
needs to be accessed by user code then it becomes dependent on the exact type of location module that it is compiled
with.

Usage of distance routines
^^^^^^^^^^^^^^^^^^^^^^^^^^

Regardless of the fact that the distance subroutine names include the string 'obs', there is nothing specific to
observations in these routines. They work to compute distances between any set of locations. The most frequent use of
these routines in the filter code is to compute the distance between a single observation and items in the state vector,
and also between a single observation and other nearby observations. However, any source for locations is supported.

In simpler location modules (like the ``oned`` version) there is no need for anything other than a brute force search
between the base location and all available state vector locations. However in the case of large geophysical models
which typically use the ``threed_sphere`` locations code, the brute-force search time is prohibitive. The location code
pre-processes all locations into a set of *bins* and then only needs to search the lists of locations in nearby bins
when looking for locations that are within a specified distance.

The expected calling sequence of the ``get_close`` routines is as follows:

::


   call get_close_maxdist_init()  ! is called before get_close_obs_init()
   call get_close_obs_init()

   call get_close_obs()           ! called many, many times

   call get_close_obs_destroy()

In the ``threed_sphere`` implementation the first routine initializes some data structures, the second one bins up the
list of locations, and then the third one is called multiple times to find all locations within a given radius of some
reference location, and to optionally compute the exact separation distance from the reference location. The last
routine deallocates the space. See the documentation below for the specific details for each routine.

All 4 of these routines must be present in every location module but in most other versions all but ``get_close_obs()``
are stubs. In this ``threed_sphere`` version of the locations module all are fully implemented.

Interaction with model_mod.f90 code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The filter and other DART programs could call the ``get_close`` routines directly, but typically do not. They declare
them (in a ``use`` statement) to be in the ``model_mod`` module, and all model interface modules are required to supply
them. However in many cases the model_mod only needs to contain another ``use`` statement declaring them to come from
the ``location_mod`` module. Thus they 'pass through' the model_mod but the user does not need to provide a subroutine
or any code for them.

However, if the model interface code wants to intercept and alter the default behavior of the get_close routines, it is
able to. Typically the model_mod still calls the location_mod routines and then adjusts the results before passing them
back to the calling code. To do that, the model_mod must be able to call the routines in the location_mod which have the
same names as the subroutines it is providing. To allow the compiler to distinguish which routine is to be called where,
we use the Fortran 90 feature which allows a module routine to be renamed in the use statement. For example, a common
case is for the model_mod to want to supply additions to the get_close_obs() routine only. At the top of the model_mod
code it would declare:

::


   use location_mod, only :: location_get_close_obs => get_close_obs,    &
                             get_close_maxdist_init, get_close_obs_init, &
                             get_close_obs_destroy

That makes calls to the maxdist_init, init, and destroy routines simply pass through to the code in the location_mod,
but the model_mod must supply a get_close_obs() subroutine. When it wants to call the code in the location_mod it calls
``location_get_close_obs()``.

One use pattern is for the model_mod to call the location get_close_obs() routine without the ``dist`` argument. This
returns a list of any potentially close locations without computing the exact distance from the base location. At this
point the list of locations is a copy and the model_mod routine is free to alter the list in any way it chooses: it can
change the locations to make certain types of locations appear closer or further away from the base location; it can
convert the vertical coordinates into a common coordinate type so that calls to the ``get_dist()`` routine can do full
3d distance computations and not just 2d (the vertical coordinates must match between the base location and the
locations in the list in order to compute a 3d distance). Then typically the model_mod code loops over the list calling
the ``get_dist()`` routine to get the actual distances to be returned to the calling code. To localize in the vertical
in a particular unit type, this is the place where the conversion to that vertical unit should be done.

Horizontal distance only
^^^^^^^^^^^^^^^^^^^^^^^^

If *horiz_distance_only* is .true. in the namelist, then the vertical coordinate is ignored and only the great-circle
distance between the two locations is computed, as if they were both on the surface of the sphere.

If *horiz_distance_only* is .false. in the namelist then the appropriate normalization constant determines the relative
impact of vertical and horizontal separation. Since only a single localization distance is specified, and the vertical
scales might have very different distance characteristics, the vert_normalization_xxx values can be used to scale the
vertical appropriately to control the desired influence of observations in the vertical.

Precomputation for run-time search efficiency
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For search efficiency all locations are pre-binned. The surface of the sphere is divided up into *nlon* by *nlat* boxes
and the index numbers of all items (both state vector entries and observations) are stored in the appropriate box. To
locate all points close to a given location, only the locations listed in the boxes within the search radius must be
checked. This speeds up the computations, for example, when localization controls which state vector items are impacted
by any given observation. The search radius is the localization distance and only those state vector items in boxes
closer than the radius to the observation location are processed.

The default values have given good performance on many of our existing model runs, but for tuning purposes the box
counts have been added to the namelist to allow adjustment. By default the code prints some summary information about
how full the average box is, how many are empty, and how many items were in the box with the largest count. The namelist
value *output_box_info* can be set to .true. to get even more information about the box statistics. The best performance
will be obtained somewhere between two extremes; the worst extreme is all the points are located in just a few boxes.
This degenerates into a (slow) linear search through the index list. The other extreme is a large number of empty or
sparsely filled boxes. The overhead of creating, managing, and searching a long list of boxes will impact performance.
The best performance lies somewhere in the middle, where each box contains a reasonable number of values, more or less
evenly distributed across boxes. The absolute numbers for best performance will certainly vary from case to case.

For latitude, the *nlat* boxes are distributed evenly across the actual extents of the data. (Locations are in radians,
so the maximum limits are the poles at :math:`-\pi/2` and :math:`+\pi/2`. For longitude, the code automatically determines if the data is
spread around more than half the sphere, and if so, the boxes are distributed evenly across the entire sphere (longitude
range :math:`0` to :math:`2\pi`). If the data spans less than half the sphere in longitude, the actual extent of the data is determined
(including correctly handling the cyclic boundary at :math:`0`) and the boxes are distributed only within the data extent. This
simplifies the actual distance calculations since the distance from the minimum longitude box to the maximum latitude
box cannot be shorter going the other way around the sphere. In practice, for a global model the boxes are evenly
distributed across the entire surface of the sphere. For local or regional models, the boxes are distributed only across
the the extent of the local grid.

For efficiency in the case where the boxes span less than half the globe, the 3D location module needs to be able to
determine the greatest longitude difference between a base point at latitude :math:`\phi_s` and all points that are separated from
that point by a central angle of :math:`\theta`. We might also want to know the latitude, :math:`\phi_f`, at which the largest separation
occurs. Note also that an intermediate form below allows the computation of the maximum longitude difference at a
particular latitude.

The central angle between a point at latitude :math:`\phi_s` and a second point at latitude :math:`\phi_f` that are separated in longitude
by :math:`\Delta\lambda` is:

.. math::

   \theta = cos^{-1}(sin\phi_s sin\phi_f + cos\phi_s cos\phi_f cos\Delta\lambda)

Taking the :math:`cos` of both sides gives:

.. math::

   cos\theta = (sin\phi_s sin\phi_f + cos\phi_s cos\phi_f cos\Delta\lambda)

Solving for :math:`cos\Delta\lambda` gives:

.. math::

   cos\Delta\lambda = \frac{a-bsin\phi_f}{c cos\phi_f}

   cos\Delta\lambda = \frac{a}{c sec\phi_f}-\frac{b}{c tan\phi_f}

where :math:`a = cos\theta`, :math:`b = sin\phi_s`, and :math:`c = cos\phi_s`. We want to maximize :math:`\Delta\lambda` which
implies minimizing :math:`cos\Delta\lambda` subject to constraints.

Taking the derivative with respect to :math:`\phi_f` gives:

.. math::

   \frac{d cos\Delta\lambda}{d\phi_f} = \frac{a}{c sec\phi_f tan\phi_f}-\frac{b}{c sec^2\phi_f}=0

Factoring out :math:`sec\phi_f` which can never be :math:`0` and using the definitions of :math:`sec` and :math:`tan` gives:

.. math::

   \frac{a sin\phi_f}{c cos\phi_f}-\frac{b}{c cos\phi_f}=0

Solving in the constrained range from :math:`0` to :math:`\pi/2` gives:

.. math::

   sin\phi_f = \frac{b}{a}=\frac{sin\phi_s}{cos\theta}

So knowing base point (:math:`\phi_s`, :math:`\lambda_s`), latitude :math:`\phi_f`, and distance :math:`\theta` we can
use the great circle equation to find the longitude difference at the greatest separation point:

.. math::

   \Delta\lambda = cos^{-1}\left(\frac{a- b sin\phi_f}{c cos\phi_f}\right)

Note that if the angle between the base point and a pole is less than or equal to the central angle, all longitude
differences will occur as the pole is approached.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand ``&`` and terminate with a slash ``/``.
Character strings that contain a ``/`` must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &location_nml
      use_octree      = .false.
      nboxes          = 1000
      maxdepth        = 4
      filled          = 10
      output_box_info = .false.
      print_box_level = 0
      compare_to_correct = .false.
   /

|

Items in this namelist either control the way in which distances are computed and/or influence the code performance.

 +--------------------+---------+-------------------------------------------------------------------------------------+
 | Item               | Type    | Description                                                                         |
 +====================+=========+=====================================================================================+
 | use_octree         | logical | There are two variations of search code. For now, this must be set to .false.       |
 +--------------------+---------+-------------------------------------------------------------------------------------+
 | nboxes             | integer | An optimization parameter which controls how many boxes the space is divided up     |
 |                    |         | into for precomputing nearby points. Larger numbers use more memory but may make    |
 |                    |         | searching faster if the model contains a large grid.                                |
 +--------------------+---------+-------------------------------------------------------------------------------------+
 | filled             | integer | An optimization parameter for the octree code. Set the lower item count limit where |
 |                    |         | a box no longer splits.                                                             |
 +--------------------+---------+-------------------------------------------------------------------------------------+
 | output_box_info    | logical | If true, print more details about the distribution of locations across the array of |
 |                    |         | boxes.                                                                              |
 +--------------------+---------+-------------------------------------------------------------------------------------+
 | print_box_level    | integer | If output_box_info is true, controls the amount of output.                          |
 +--------------------+---------+-------------------------------------------------------------------------------------+
 | compare_to_correct | logical | If true do an exhaustive (and slow) search to ensure the results are the same as    |
 |                    |         | using optimized search code. Should only be used for debugging.                     |
 +--------------------+---------+-------------------------------------------------------------------------------------+

|

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
\                            vert_is_scale_height
\                            vert_is_level
\                            vert_is_height
\                            VERTISUNDEF
\                            VERTISSURFACE
\                            VERTISLEVEL
\                            VERTISPRESSURE
\                            VERTISHEIGHT
\                            VERTISSCALEHEIGHT
\                            operator(==)
\                            operator(/=)
============================ ======================

Namelist interface ``&location_nml`` must be read from file ``input.nml``.

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

|

.. container:: type

   *type location_type*
   ::

         private
         real(r8) :: lon, lat, vloc
         integer  :: which_vert
      end type location_type

.. container:: indent1

   Provides an abstract representation of physical location on a three-d spherical shell.

   +------------+--------------------------------------------------------------------------------------------------------+
   | Component  | Description                                                                                            |
   +============+========================================================================================================+
   | lon        | longitude in radians                                                                                   |
   +------------+--------------------------------------------------------------------------------------------------------+
   | lat        | latitude in radians                                                                                    |
   +------------+--------------------------------------------------------------------------------------------------------+
   | vloc       | vertical location, units as selected by which_vert                                                     |
   +------------+--------------------------------------------------------------------------------------------------------+
   | which_vert | type of vertical location: -2=no specific vert location; -1=surface; 1=level; 2=pressure; 3=height,    |
   |            | 4=scale height                                                                                         |
   +------------+--------------------------------------------------------------------------------------------------------+

   The vertical types have parameters defined for them so they can be referenced by name instead of number.

|

.. container:: type

   *type get_close_type*
   ::

         private
         integer  :: num
         real(r8) :: maxdist
         integer, pointer :: lon_offset(:, :)
         integer, pointer :: obs_box(:)
         integer, pointer :: count(:, :)
         integer, pointer :: start(:, :)
      end type get_close_type

.. container:: indent1

   Provides a structure for doing efficient computation of close locations.

   +------------+--------------------------------------------------------------------------------------------------------+
   | Component  | Description                                                                                            |
   +============+========================================================================================================+
   | num        | Number of locations in list                                                                            |
   +------------+--------------------------------------------------------------------------------------------------------+
   | maxdist    | Threshhold distance. Anything closer is close.                                                         |
   +------------+--------------------------------------------------------------------------------------------------------+
   | lon_offset | Dimensioned nlon by nlat. For a given offset in longitude boxes and difference in latitudes, gives max |
   |            | distance from base box to a point in offset box.                                                       |
   +------------+--------------------------------------------------------------------------------------------------------+
   | obs_box    | Dimensioned num. Gives index of what box each location is in.                                          |
   +------------+--------------------------------------------------------------------------------------------------------+
   | count      | Dimensioned nlon by nlat. Number of obs in each box.                                                   |
   +------------+--------------------------------------------------------------------------------------------------------+
   | start      | Dimensioned nlon by nlat. Index in straight storage list where obs in each box start.                  |
   +------------+--------------------------------------------------------------------------------------------------------+

|

.. container:: routine

   *var = get_location(loc)*
   ::

      real(r8), dimension(3)          :: get_location
      type(location_type), intent(in) :: loc

.. container:: indent1

   Extracts the longitude and latitude (converted to degrees) and the vertical location from a location type and returns
   in a 3 element real array.

   ================ =============================================================
   ``get_location`` The longitude and latitude (in degrees) and vertical location
   ``loc``          A location type
   ================ =============================================================

|

.. container:: routine

   *var = set_location(lon, lat, vert_loc, which_vert)*
   ::

      type(location_type)   :: set_location
      real(r8), intent(in)    :: lon
      real(r8), intent(in)    :: lat
      real(r8), intent(in)    :: vert_loc
      integer,  intent(in)    :: which_vert

.. container:: indent1

   Returns a location type with the input longitude and latitude (input in degrees) and the vertical location of type
   specified by which_vert.

   ================ ============================================
   ``set_location`` A location type
   ``lon``          Longitude in degrees
   ``lat``          Latitude in degrees
   ``vert_loc``     Vertical location consistent with which_vert
   ``which_vert``   The vertical location type
   ================ ============================================

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

   Returns the value of which_vert, latitude, longitude, or vertical location from a location type as selected by the
   string argument attr. If attr is not present or if it is 'WHICH_VERT', the value of which_vert is converted to real
   and returned. Otherwise, attr='LON' returns longitude, attr='LAT' returns latitude and attr='VLOC' returns the
   vertical location.

   ================== =================================================================================
   ``query_location`` Returns longitude, latitude, vertical location, or which_vert (converted to real)
   ``loc``            A location type
   *attr*             Selects 'WHICH_VERT', 'LON', 'LAT' or 'VLOC'
   ================== =================================================================================

|

.. container:: routine

   *var = set_location_missing()*
   ::

      type(location_type) :: set_location_missing

.. container:: indent1

   Returns a location with all elements set to missing values defined in types module.

   ======================== ==================================================
   ``set_location_missing`` A location with all elements set to missing values
   ======================== ==================================================

|

.. container:: routine

   *call get_close_maxdist_init(gc,maxdist, [maxdist_list])*
   ::

      type(get_close_type), intent(inout) :: gc
      real(r8), intent(in)                :: maxdist
      real(r8), intent(in), optional      :: maxdist_list(:)

.. container:: indent1

   Sets the threshhold distance. ``maxdist`` is in units of radians. Anything closer than this is deemed to be close.
   This routine must be called first, before the other ``get_close`` routines. It allocates space so it is necessary to
   call ``get_close_obs_destroy`` when completely done with getting distances between locations.

   If the last optional argument is not specified, maxdist applies to all locations. If the last argument is specified,
   it must be a list of exactly the length of the number of specific types in the obs_kind_mod.f90 file. This length can
   be queried with the `get_num_types_of_obs() <../../modules/observations/obs_kind_mod.html#get_num_types_of_obs>`__
   function to get count of obs types. It allows a different maximum distance to be set per base type when get_close()
   is called.

   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``gc``      | Data for efficiently finding close locations.                                                         |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``maxdist`` | Anything closer than this number of radians is a close location.                                      |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | *maxdist*   | If specified, must be a list of real values. The length of the list must be exactly the same length   |
   |             | as the number of observation types defined in the obs_def_kind.f90 file. (See                         |
   |             | `get_num_types_of_obs() <../../modules/observations/obs_kind_mod.html#get_num_types_of_obs>`__ to get |
   |             | count of obs types.) The values in this list are used for the obs types as the close distance instead |
   |             | of the maxdist argument.                                                                              |
   +-------------+-------------------------------------------------------------------------------------------------------+

|

.. container:: routine

   *call get_close_obs_init(gc, num, obs)*
   ::

      type(get_close_type),             intent(inout) :: gc
      integer,                          intent(in)    :: num
      type(location_type), dimension(:) intent(in)    :: obs

.. container:: indent1

   Initialize storage for efficient identification of locations close to a given location. Allocates storage for keeping
   track of which 'box' each location in the list is in. Must be called after ``get_close_maxdist_init``, and the list
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
      integer,             dimension(:), intent(in)  :: obs_kind
      integer,                           intent(out) :: num_close
      integer,             dimension(:), intent(out) :: close_ind
      real(r8), optional,  dimension(:), intent(out) :: dist

.. container:: indent1

   Given a single location and a list of other locations, returns the indices of all the locations close to the single
   one along with the number of these and the distances for the close ones. The list of locations passed in via the
   ``obs`` argument must be identical to the list of ``obs`` passed into the most recent call to
   ``get_close_obs_init()``. If the list of locations of interest changes ``get_close_obs_destroy()`` must be called and
   then the two initialization routines must be called before using ``get_close_obs()`` again.

   If called without the optional *dist* argument, all locations that are potentially close are returned, which is
   likely a superset of the locations that are within the threshold distance specified in the
   ``get_close_maxdist_init()`` call. This can be useful to collect a list of potential locations, and then to convert
   all the vertical coordinates into one consistent unit (pressure, height in meters, etc), and then the list can be
   looped over, calling get_dist() directly to get the exact distance, either including vertical or not depending on the
   setting of ``horiz_dist_only``.

   ================= ===================================================================================
   ``gc``            Structure to allow efficient identification of locations close to a given location.
   ``base_obs_loc``  Single given location.
   ``base_obs_kind`` Kind of the single location.
   ``obs``           List of locations from which close ones are to be found.
   ``obs_kind``      Kind associated with locations in obs list.
   ``num_close``     Number of locations close to the given location.
   ``close_ind``     Indices of those locations that are close.
   *dist*            Distance between given location and the close ones identified in close_ind.
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

   *var = get_dist(loc1, loc2, [, kind1, kind2, no_vert])*
   ::

      real(r8)                        :: get_dist
      type(location_type), intent(in) :: loc1
      type(location_type), intent(in) :: loc2
      integer, optional,   intent(in) :: kind1
      integer, optional,   intent(in) :: kind2
      logical, optional,   intent(in) :: no_vert

.. container:: indent1

   Returns the distance between two locations in radians. If ``horiz_dist_only`` is set to .TRUE. in the locations
   namelist, it computes great circle distance on sphere. If ``horiz_dist_only`` is false, then it computes an
   ellipsoidal distance with the horizontal component as above and the vertical distance determined by the types of the
   locations and the normalization constants set by the namelist for the different vertical coordinate types. The
   vertical normalization gives the vertical distance that is equally weighted as a horizontal distance of 1 radian. If
   *no_vert* is present, it overrides the value in the namelist and controls whether vertical distance is included or
   not.

   The kind arguments are not used by the default location code, but are available to any user-supplied distance
   routines which want to do specialized calculations based on the kinds associated with each of the two locations.

   ========= =====================================================================================
   ``loc1``  First of two locations to compute distance between.
   ``loc2``  Second of two locations to compute distance between.
   *kind1*   DART kind associated with location 1.
   *kind2*   DART kind associated with location 2.
   *no_vert* If true, no vertical component to distance. If false, vertical component is included.
   ``var``   distance between loc1 and loc2.
   ========= =====================================================================================

|

.. container:: routine

   *var = vert_is_undef(loc)*
   ::

      logical                         :: vert_is_undef
      type(location_type), intent(in) :: loc

.. container:: indent1

   Returns true if which_vert is set to undefined, else false. The meaning of 'undefined' is specific; it means there is
   no particular vertical location associated with this type of measurement; for example a column-integrated value.

   ================= ========================================================
   ``vert_is_undef`` Returns true if vertical coordinate is set to undefined.
   ``loc``           A location type
   ================= ========================================================

|

.. container:: routine

   *var = vert_is_surface(loc)*
   ::

      logical                         :: vert_is_surface
      type(location_type), intent(in) :: loc

.. container:: indent1

   Returns true if which_vert is for surface, else false.

   =================== ===================================================
   ``vert_is_surface`` Returns true if vertical coordinate type is surface
   ``loc``             A location type
   =================== ===================================================

|

.. container:: routine

   *var = vert_is_pressure(loc)*
   ::

      logical                         :: vert_is_pressure
      type(location_type), intent(in) :: loc

.. container:: indent1

   Returns true if which_vert is for pressure, else false.

   ==================== ====================================================
   ``vert_is_pressure`` Returns true if vertical coordinate type is pressure
   ``loc``              A location type
   ==================== ====================================================

|

.. container:: routine

   *var = vert_is_scale_height(loc)*
   ::

      logical                         :: vert_is_scale_height
      type(location_type), intent(in) :: loc

.. container:: indent1

   Returns true if which_vert is for scale_height, else false.

   ======================== ========================================================
   ``vert_is_scale_height`` Returns true if vertical coordinate type is scale_height
   ``loc``                  A location type
   ======================== ========================================================

|

.. container:: routine

   *var = vert_is_level(loc)*
   ::

      logical                         :: vert_is_level
      type(location_type), intent(in) :: loc

.. container:: indent1

   Returns true if which_vert is for level, else false.

   ================= =================================================
   ``vert_is_level`` Returns true if vertical coordinate type is level
   ``loc``           A location type
   ================= =================================================

|

.. container:: routine

   *var = vert_is_height(loc)*
   ::

      logical                         :: vert_is_height
      type(location_type), intent(in) :: loc

.. container:: indent1

   Returns true if which_vert is for height, else false.

   ================== ==================================================
   ``vert_is_height`` Returns true if vertical coordinate type is height
   ``loc``            A location type
   ================== ==================================================

|

.. container:: routine

   *var = has_vertical_localization()*
   ::

      logical :: has_vertical_localization

.. container:: indent1

   Returns .TRUE. if the namelist variable ``horiz_dist_only`` is .FALSE. meaning that vertical separation between
   locations is going to be computed by ``get_dist()`` and by ``get_close_obs()``.

   This routine should perhaps be renamed to something like 'using_vertical_for_distance' or something similar. The
   current use for it is in the localization code inside filter, but that doesn't make this a representative function
   name. And at least in current usage, returning the opposite setting of the namelist item makes the code read more
   direct (fewer double negatives).

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

      integer, parameter :: VERTISUNDEF       = -2
      integer, parameter :: VERTISSURFACE     = -1
      integer, parameter :: VERTISLEVEL       =  1
      integer, parameter :: VERTISPRESSURE    =  2
      integer, parameter :: VERTISHEIGHT      =  3
      integer, parameter :: VERTISSCALEHEIGHT =  4

.. container:: indent1

   Constant parameters used to differentiate vertical types.

|

.. container:: routine

   ::

      integer, parameter :: LocationDims = 3

.. container:: indent1

   This is a **constant**. Contains the number of real values in a location type. Useful for output routines that must
   deal transparently with many different location modules.

|

.. container:: routine

   ::

      character(len=129), parameter :: LocationName = "loc3Dsphere"

.. container:: indent1

   This is a **constant**. A parameter to identify this location module in output metadata.

|

.. container:: routine

   ::

      character(len=129), parameter :: LocationLName =

             "threed sphere locations: lon, lat, vertical"

.. container:: indent1

   This is a **constant**. A parameter set to "threed sphere locations: lon, lat, vertical" used to identify this
   location module in output long name metadata.

|

Files
-----

========= =================================
filename  purpose
========= =================================
input.nml to read the location_mod namelist
========= =================================

References
----------

#. none

Private components
------------------

N/A
