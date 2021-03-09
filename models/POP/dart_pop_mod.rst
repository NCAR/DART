MODULE dart_pop_mod (POP)
=========================

Overview
--------

``dart_pop_mod`` provides a consistent collection of routines that are useful for multiple programs e.g.
``dart_to_pop``, ``pop_to_dart``, etc.

Namelist
--------

| There are no namelists unique to this module. It is necessary for this module to read some of the POP namelists, and
  so they are declared in this module. In one instance, DART will read the ``time_manager_nml`` namelist and **write**
  an updated version to control the length of the integration of POP. All other information is simply read from the
  namelists and is used in the same context as POP itself. The POP documentation should be consulted. **Only the
  variables of interest to DART are described in this document.**
| All namelists are read from a file named ``pop_in``.

.. container:: namelist

   ::

      namelist /time_manager_nml/  allow_leapyear, stop_count, stop_option

.. container:: indent1

   ``dart_to_pop`` controls the model advance of LANL/POP by creating a ``&time_manager_nml`` in ``pop_in.DART`` **IFF**
   the DART state being converted has the 'advance_to_time' record. The ``pop_in.DART`` must be concatenated with the
   other namelists needed by POP into a file called ``pop_in`` . We have chosen to store the other namelists (which
   contain static information) in a file called ``pop_in.part2``. Initially, the ``time_manager_nml`` is stored in a
   companion file called ``pop_in.part1`` and the two files are concatenated into the expected ``pop_in`` - then, during
   the course of an assimilation experiment, DART keeps writing out a new ``time_manager_nml`` with new integration
   information - which gets appended with the static information in ``pop_in.part2``

   +----------------+-------------------+-------------------------------------------------------------------------------+
   | Contents       | Type              | Description                                                                   |
   +================+===================+===============================================================================+
   | allow_leapyear | logical           | DART ignores the setting of this parameter. All observations must use a       |
   |                |                   | Gregorian calendar. There are pathological cases, but if you are doing data   |
   |                |                   | assimilation, just use the Gregorian calendar. *[default: .true.]*            |
   +----------------+-------------------+-------------------------------------------------------------------------------+
   | stop_count     | integer           | the number of model advance steps to take. *[default: 1]*                     |
   +----------------+-------------------+-------------------------------------------------------------------------------+
   | stop_option    | character(len=64) | The units for the number of model advance steps (``stop_count``) to take.     |
   |                |                   | *[default: 'ndays']*                                                          |
   +----------------+-------------------+-------------------------------------------------------------------------------+

| 

.. container:: namelist

   ::

      namelist /io_nml/  luse_pointer_files, pointer_filename

.. container:: indent1

   +--------------------+--------------------+--------------------------------------------------------------------------+
   | Contents           | Type               | Description                                                              |
   +====================+====================+==========================================================================+
   | luse_pointer_files | logical            | switch to indicate the use of pointer files or not. If ``.true.``, a     |
   |                    |                    | pointer file is used to contain the name of the restart file to be used. |
   |                    |                    | DART requires this to be ``.true``. *[default: .true.]*                  |
   +--------------------+--------------------+--------------------------------------------------------------------------+
   | pointer_filename   | character(len=100) | The name of the pointer file. All of the DART scripts presume and        |
   |                    |                    | require the use of the default. Each ensmeble member gets its own        |
   |                    |                    | pointer file. *[default: rpointer.ocn.[1-N].restart]*                    |
   +--------------------+--------------------+--------------------------------------------------------------------------+

| 

.. container:: namelist

   ::

      namelist /restart_nml/  restart_freq_opt, restart_freq

.. container:: indent1

   +--------------------+--------------------+--------------------------------------------------------------------------+
   | Contents           | Type               | Description                                                              |
   +====================+====================+==========================================================================+
   | luse_pointer_files | logical            | switch to indicate the use of pointer files or not. If ``.true.``, a     |
   |                    |                    | pointer file is used to contain the name of the restart file to be used. |
   |                    |                    | DART requires this to be ``.true``. *[default: .true.]*                  |
   +--------------------+--------------------+--------------------------------------------------------------------------+
   | pointer_filename   | character(len=100) | The name of the pointer file. All of the DART scripts presume and        |
   |                    |                    | require the use of the default. Each ensmeble member gets its own        |
   |                    |                    | pointer file. *[default: rpointer.ocn.[1-N].restart]*                    |
   +--------------------+--------------------+--------------------------------------------------------------------------+

| 

.. container:: namelist

   ::

      namelist /init_ts_nml/  init_ts_option, init_ts_file, init_ts_file_fmt

.. container:: indent1

   The ``dart_pop_mod:initialize_module()`` routine reads ``pop_in`` . There are several code stubs for future use that
   may allow for a more fully-supported POP namelist implementation. This namelist is one of them. Until further notice,
   the ``init_ts_nml`` is completely ignored by DART.

   +------------------+--------------------+----------------------------------------------------------------------------+
   | Contents         | Type               | Description                                                                |
   +==================+====================+============================================================================+
   | init_ts_option   | character(len=64)  | NOT USED by DART. All T,S information comes from a netCDF restart file     |
   |                  |                    | named ``pop.r.nc`` *[default: 'restart']*                                  |
   +------------------+--------------------+----------------------------------------------------------------------------+
   | init_ts_file     | character(len=100) | NOT USED by DART. All T,S information comes from ``pop.r.nc`` *[default:   |
   |                  |                    | 'pop.r']*                                                                  |
   +------------------+--------------------+----------------------------------------------------------------------------+
   | init_ts_file_fmt | character(len=64)  | NOT USED by DART. The file format is ``'nc'`` *[default: 'nc']*            |
   +------------------+--------------------+----------------------------------------------------------------------------+

| 

.. container:: namelist

   ::

      namelist /domain_nml/  ew_boundary_type

.. container:: indent1

   DART needs to know if the East-West domain is cyclic for spatial interpolations. Presently, DART has only been tested
   for the dipole grid, which is cyclic E-W and closed N-S.

   +------------------+-------------------+-----------------------------------------------------------------------------+
   | Contents         | Type              | Description                                                                 |
   +==================+===================+=============================================================================+
   | ew_boundary_type | character(len=64) | switch to indicate whether the East-West domain is cyclic or not. DART/POP  |
   |                  |                   | has not been tested in a regional configuration, so DART requires this to   |
   |                  |                   | be ``'cyclic'``. *[default: 'cyclic']*                                      |
   +------------------+-------------------+-----------------------------------------------------------------------------+

| 

.. container:: namelist

   ::

      namelist /grid_nml/  horiz_grid_opt,  vert_grid_opt,  topography_opt, &
                                     horiz_grid_file, vert_grid_file, topography_file

.. container:: indent1

   | The POP grid information comes in several files: horizontal grid lat/lons in one, the vertical grid spacing in
     another, and the topography (lowest valid vertical level) in a third.
   | Here is what we can get from the (binary) horizontal grid file:

   ::

      real(r8), dimension(:,:) :: ULAT,  &! latitude  (radians) of U points
      real(r8), dimension(:,:) :: ULON,  &! longitude (radians) of U points
      real(r8), dimension(:,:) :: HTN ,  &! length (cm) of north edge of T box
      real(r8), dimension(:,:) :: HTE ,  &! length (cm) of east  edge of T box
      real(r8), dimension(:,:) :: HUS ,  &! length (cm) of south edge of U box
      real(r8), dimension(:,:) :: HUW ,  &! length (cm) of west  edge of U box
      real(r8), dimension(:,:) :: ANGLE  &! angle

   The vertical grid file is ascii, with 3 columns/line:

   ::

      cell thickness(in cm)   cell center(in m)   cell bottom(in m)

   Here is what we can get from the topography file:

   ::

      integer, dimension(:,:), :: KMT    &! k index of deepest grid cell on T grid

   These must be derived or come from someplace else ...

   ::

      KMU               k index of deepest grid cell on U grid
      HT                real(r8) value of deepest valid T depth (in cm)
      HU                real(r8) value of deepest valid U depth (in cm)

   +-----------------------------------------------+--------------------+-----------------------------------------------+
   | Contents                                      | Type               | Description                                   |
   +===============================================+====================+===============================================+
   | horiz_grid_opt, vert_grid_opt, topography_opt | character(len=64)  | switch to indicate whether or not the grids   |
   |                                               |                    | will come from an external file or not. DART  |
   |                                               |                    | requires ALL of these to be ``'file'``.       |
   |                                               |                    | *[default: 'file']*                           |
   +-----------------------------------------------+--------------------+-----------------------------------------------+
   | horiz_grid_file                               | character(len=100) | The name of the binary file containing the    |
   |                                               |                    | values for the horizontal grid. The           |
   |                                               |                    | **dimensions** of the grid are read from      |
   |                                               |                    | ``pop.r.nc``. It would have been nice to      |
   |                                               |                    | include the actual grid information in the    |
   |                                               |                    | netCDF files. *[default:                      |
   |                                               |                    | 'horiz_grid.gx3v5.r8ieee.le']*                |
   +-----------------------------------------------+--------------------+-----------------------------------------------+
   | vert_grid_file                                | character(len=100) | The name of the ASCII file containing the     |
   |                                               |                    | values for the vertical grid. The file must   |
   |                                               |                    | contain three columns of data pertaining to   |
   |                                               |                    | the cell thickness (in cm), the cell center   |
   |                                               |                    | (in meters), and the cell bottom (in meters). |
   |                                               |                    | Again, it would have been nice to include the |
   |                                               |                    | vertical grid information in the netCDF       |
   |                                               |                    | files. *[default: 'vert_grid.gx3v5']*         |
   +-----------------------------------------------+--------------------+-----------------------------------------------+
   | topography_grid_file                          | character(len=100) | The name of the binary file containing the    |
   |                                               |                    | values for the topography information. The    |
   |                                               |                    | **dimensions** of the grid are read from      |
   |                                               |                    | ``pop.r.nc``. *[default:                      |
   |                                               |                    | 'topography.gx3v5.r8ieee.le']*                |
   +-----------------------------------------------+--------------------+-----------------------------------------------+

| 

Other modules used
------------------

::

   types_mod
   time_manager_mod
   utilities_mod
   typesizes
   netcdf

Public interfaces
-----------------

Only a select number of interfaces used are discussed here. Each module has its own discussion of their routines.

Interface routines
~~~~~~~~~~~~~~~~~~

========================== ========================
*use dart_pop_mod, only :* get_pop_calendar
\                          set_model_time_step
\                          get_horiz_grid_dims
\                          get_vert_grid_dim
\                          read_horiz_grid
\                          read_topography
\                          read_vert_grid
\                          write_pop_namelist
\                          get_pop_restart_filename
========================== ========================

Required interface routines
~~~~~~~~~~~~~~~~~~~~~~~~~~~

| 

.. container:: routine

   *call get_pop_calendar(calstring)*
   ::

      character(len=*), intent(out) :: calstring

.. container:: indent1

   Returns a string containing the type of calendar in use.

   ============= =====================================
   ``calstring`` DART/POP uses a 'gregorian' calendar.
   ============= =====================================

| 

.. container:: routine

   *poptimestep = set_model_time_step()*
   ::

      type(time_type), intent(out) :: poptimestep

.. container:: indent1

   ``set_model_time_step`` returns the model time step that was set in the restart_nml\ ``restart_freq``. This is the
   minimum amount of time DART thinks the POP model can advance. Indirectly, this specifies the minimum assimilation
   interval.

   =============== =================================
   ``poptimestep`` the minimum assimilation interval
   =============== =================================

| 

.. container:: routine

   *call get_horiz_grid_dims(Nx, Ny)*
   ::

      integer, intent(out) :: Nx, Ny

.. container:: indent1

   ``get_horiz_grid_dims`` reads ``pop.r.nc`` to determine the number of longitudes and latitudes.

   ====== =========================================================================================
   ``Nx`` the length of the 'i' dimension in the POP restart file. The number of longitudes in use.
   ``Ny`` the length of the 'j' dimension in the POP restart file. The number of latitudes in use.
   ====== =========================================================================================

| 

.. container:: routine

   *call get_vert_grid_dim( Nz )*
   ::

      integer, intent(out) :: Nz

.. container:: indent1

   ``get_vert_grid_dim`` reads ``pop.r.nc`` to determine the number of vertical levels in use.

   ====== ==============================================================================================
   ``Nz`` the length of the 'k' dimension in the POP restart file. The number of vertical levels in use.
   ====== ==============================================================================================

| 

.. container:: routine

   *call read_horiz_grid(nx, ny, ULAT, ULON, TLAT, TLON)*
   ::

      integer,                    intent(in)  :: nx, ny
      real(r8), dimension(nx,ny), intent(out) :: ULAT, ULON, TLAT, TLON

.. container:: indent1

   ``read_horiz_grid`` reads the direct access binary files containing the POP grid information. **The first record is
   REQUIRED to be 'ULAT', the second record is REQUIRED to be 'ULON'.**

   ======== ====================================================================================
   ``nx``   The number of longitudes in the grid.
   ``ny``   The number of latitudes in the grid.
   ``ULAT`` The matrix of latitudes for the UVEL and VVEL variables. Units are degrees [-90,90].
   ``ULON`` The matrix of longitudes for the UVEL and VVEL variables. Units are degrees. [0,360]
   ``TLAT`` The matrix of latitudes for the SALT and TEMP variables. Units are degrees [-90,90].
   ``TLON`` The matrix of longitudes for the SALT and TEMP variables. Units are degrees. [0,360]
   ======== ====================================================================================

| 

.. container:: routine

   *call read_topography(nx, ny, KMT, KMU)*
   ::

      integer,                   intent(in)  :: nx, ny
      integer, dimension(nx,ny), intent(out) :: KMT, KMU

.. container:: indent1

   ``read_topography`` reads the direct access binary files containing the POP topography information. **The first
   record is REQUIRED to be 'KMT'.** 'KMU' is calculated from 'KMT'.

   ======= =====================================================================
   ``nx``  The number of longitudes in the grid.
   ``ny``  The number of latitudes in the grid.
   ``KMT`` The matrix containing the lowest valid depth index at grid centroids.
   ``KMU`` The matrix containing the lowest valid depth index at grid corners.
   ======= =====================================================================

| 

.. container:: routine

   *call read_vert_grid(nz, ZC, ZG)*
   ::

      integer,                 intent(in)  :: nz
      real(r8), dimension(nz), intent(out) :: ZC, ZG

.. container:: indent1

   | ``read_vert_grid`` reads the ASCII file containing the information about the vertical levels. The file must contain
     three columns of data pertaining to; 1) the cell thickness (in cm),
   | 2) the cell center (in meters),
   | and 3) the cell bottom (in meters).

   ====== ==========================================
   ``nz`` The number of vertical levels.
   ``ZC`` The depth (in meters) at the grid centers.
   ``ZG`` The depth (in meters) at the grid edges.
   ====== ==========================================

| 

.. container:: routine

   *call write_pop_namelist(model_time, adv_to_time)*
   ::

      type(time_type), intent(in)  :: model_time
      type(time_type), intent(in)  :: adv_to_time

.. container:: indent1

   ``write_pop_namelist`` writes the POP namelist ``time_manager_nml`` with the information necessary to advance POP to
   the next assimilation time. The namelist is written to a file called ``pop_in.DART``. Presently, DART is configured
   to minimally advance POP for 86400 seconds - i.e. 1 day. The forecast length (the difference between 'model_time' and
   'adv_to_time') must be an integer number of days with the current setup. An error will result if it is not.

   =============== ============================================
   ``model_time``  The 'valid' time of the current model state.
   ``adv_to_time`` The time of the next assimilation.
   =============== ============================================

| 

.. container:: routine

   *call get_pop_restart_filename( filename )*
   ::

      character(len=*), intent(out) :: filename

.. container:: indent1

   ``get_pop_restart_filename`` returns the filename containing the POP restart information. At this point the filename
   is **hardwired** to ``pop.r.nc``, but may become more flexible in future versions. The filename may be derived from
   the ``restart_nml`` but is currently ignored.

   ============ =================================
   ``filename`` The name of the POP restart file.
   ============ =================================

| 

Files
-----

==================================== ============================================================
filename                             purpose
==================================== ============================================================
pop_in                               to read the POP namelists
pop.r.nc                             provides grid dimensions and 'valid_time' of the model state
``&grid_nml`` "horiz_grid_file"      contains the values of the horizontal grid
``&grid_nml`` "vert_grid_file"       contains the number and values of the vertical levels
``&grid_nml`` "topography_grid_file" contains the indices of the wet/dry cells
pop_in.DART                          to control the integration of the POP model advance
==================================== ============================================================

| 

References
----------

-  none

Private components
------------------

N/A
