MPAS OCN
========

.. attention::

   ``mpas_ocn`` was being developed with versions of DART *before* Manhattan (9.x.x) and has yet to be updated. If you are interested in
   using ``mpas_ocn`` with more recent versions of DART, contact DAReS staff to assess the feasibility of an update.
   Until that time, you should consider this documentation as out-of-date.


Overview
--------

| The **MPAS OCN** interface for **Data Assimilation Research Testbed (DART)** is under development.
| Since MPAS OCN uses netcdf files for their restart mechanism, a namelist-controlled set of variables is used to build
  the DART state vector. Each variable must also correspond to a DART "QUANTITY"; required for the DART interpolate
  routines. For example:

::

   &mpas_vars_nml
      mpas_state_variables = 'uReconstructZonal',      'QTY_U_WIND_COMPONENT',
                             'uReconstructMeridional', 'QTY_V_WIND_COMPONENT',
                             'w',                      'QTY_VERTICAL_VELOCITY',
                             'theta',                  'QTY_POTENTIAL_TEMPERATURE',
                             'qv',                     'QTY_VAPOR_MIXING_RATIO',
                             'qc',                     'QTY_CLOUDWATER_MIXING_RATIO',
                             'qr',                     'QTY_RAINWATER_MIXING_RATIO',
                             'qi',                     'QTY_ICE_MIXING_RATIO',
                             'qs',                     'QTY_SNOW_MIXING_RATIO',
                             'qg',                     'QTY_GRAUPEL_MIXING_RATIO',
                             'surface_pressure',       'QTY_SURFACE_PRESSURE'
      /

| These variables are then adjusted to be consistent with observations and stuffed back into the same netCDF analysis
  files. Since DART is an ensemble algorithm, there are multiple analysis files for a single analysis time: one for each
  ensemble member. Creating the initial ensemble of states is an area of active research.
| DART reads grid information from the MPAS OCN 'history' file, I have tried to keep the variable names the same.
  Internal to the DART code, the following variables exist:

============================= ===========================================================
integer :: nCells             the number of Cell Centers
integer :: nEdges             the number of Cell Edges
integer :: nVertices          the number of Cell Vertices
integer :: nVertLevels        the number of vertical level midpoints
integer :: nVertLevelsP1      the number of vertical level edges
integer :: nSoilLevels        the number of soil level ?midpoints?
real(r8) :: latCell(:)        the latitudes of the Cell Centers (-90,90)
real(r8) :: lonCell(:)        the longitudes of the Cell Centers [0, 360)
real(r8) :: zgrid(:,:)        cell center geometric height at cell centers (ncells,nvert)
integer :: CellsOnVertex(:,:) list of cell centers defining a triangle
============================= ===========================================================

model_mod variable storage
--------------------------

``input.nml``\ ``&mpas_vars_nml`` defines the list of MPAS variables used to build the DART state vector. Combined with
an MPAS analysis file, the information is used to determine the size of the DART state vector and derive the metadata.
To keep track of what variables are contained in the DART state vector, an array of a user-defined type called "progvar"
is available with the following components:

.. container:: unix

   ::

      type progvartype
         private
         character(len=NF90_MAX_NAME) :: varname
         character(len=NF90_MAX_NAME) :: long_name
         character(len=NF90_MAX_NAME) :: units
         character(len=NF90_MAX_NAME), dimension(NF90_MAX_VAR_DIMS) :: dimname
         integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
         integer :: xtype         ! netCDF variable type (NF90_double, etc.) 
         integer :: numdims       ! number of dims - excluding TIME
         integer :: numvertical   ! number of vertical levels in variable
         integer :: numcells      ! number of horizontal locations (typically cell centers)
         logical :: ZonHalf       ! vertical coordinate has dimension nVertLevels
         integer :: varsize       ! prod(dimlens(1:numdims))
         integer :: index1        ! location in dart state vector of first occurrence
         integer :: indexN        ! location in dart state vector of last  occurrence
         integer :: dart_kind
         character(len=paramname_length) :: kind_string
         logical  :: clamping     ! does variable need to be range-restricted before 
         real(r8) :: range(2)     ! being stuffed back into MPAS analysis file.
      end type progvartype

      type(progvartype), dimension(max_state_variables) :: progvar

The variables are simply read from the MPAS analysis file and stored in the DART state vector such that all quantities
for one variable are stored contiguously. Within each variable; they are stored vertically-contiguous for each
horizontal location. From a storage standpoint, this would be equivalent to a Fortran variable dimensioned
x(nVertical,nHorizontal,nVariables). The fastest-varying dimension is vertical, then horizontal, then variable ...
naturally, the DART state vector is 1D. Each variable is also stored this way in the MPAS analysis file.

.. container:: indent1

   .. rubric:: The DART interface for MPAS (atm)
      :name: the-dart-interface-for-mpas-atm

   | was compiled with the gfortran 4.2.3 compilers and run on a Mac.
   | The DART components were built with the following ``mkmf.template`` settings:

   ::

            FC = gfortran
            LD = gfortran
            NETCDF = /Users/thoar/GNU
            INCS = -I${NETCDF}/include
            LIBS = -L${NETCDF}/lib -lnetcdf -lcurl -lhdf5_hl -lhdf5 -lz  -lm
            FFLAGS = -O0 -fbounds-check -frecord-marker=4 -ffpe-trap=invalid $(INCS)
            LDFLAGS = $(FFLAGS) $(LIBS)
         

.. container:: indent1

   .. rubric:: Converting between DART files and MPAS analysis files
      :name: converting-between-dart-files-and-mpas-analysis-files

   | is relatively straighforward. Given the namelist mechanism for determining the state variables and the MPAS history
     netCDF files exist, - everything that is needed is readily determined.
   | There are two programs - both require the list of MPAS variables to use in the DART state vector: the
     ``mpas_vars_nml`` namelist in the ``input.nml`` file. The MPAS file name being read and/or written is - in all
     instances - specified by the ``model_nml:model_analysis_filename`` variable in the ``input.nml`` namelist file.

   +-------------------------------------------+-------------------------------------------------------------------------+
   | :doc:`./model_to_dart`                    | converts an MPAS analysis file (nominally named ``mpas_analysis.nc``)   |
   |                                           | into a DART-compatible file normally called ``dart_ics`` . We usually   |
   |                                           | wind up linking the actual analysis file to a static name that is used  |
   |                                           | by DART.                                                                |
   +-------------------------------------------+-------------------------------------------------------------------------+
   | `dart_to_model.f90 <dart_to_model.f90>`__ | inserts the DART output into an existing MPAS analysis netCDF file by   |
   |                                           | overwriting the variables in the analysis netCDF file. There are two    |
   |                                           | different types of DART output files, so there is a namelist option to  |
   |                                           | specify if the DART file has two time records or just one (if there are |
   |                                           | two, the first one is the 'advance_to' time, followed by the            |
   |                                           | 'valid_time' of the ensuing state). ``dart_to_model`` updates the MPAS  |
   |                                           | analysis file specified in                                              |
   |                                           | ``input.nml``\ ``model_nml:model_analysis_filename``. If the DART file  |
   |                                           | contains an 'advance_to' time, separate control information is written  |
   |                                           | to an auxiliary file that is used by the ``advance_model.csh`` script.  |
   +-------------------------------------------+-------------------------------------------------------------------------+

.. container:: indent1

   The header of an MPAS analysis file is presented below - simply for context. Keep in mind that **many** variables
   have been removed for clarity. Also keep in mind that the multi-dimensional arrays listed below have the dimensions
   reversed from the Fortran convention.

.. container:: unix

   ::

      366 mirage2:thoar% ncdump -h mpas_analysis.nc
      netcdf mpas_analysis {
      dimensions:
              StrLen = 64 ;
              Time = UNLIMITED ; // (1 currently)
              nCells = 10242 ;                                  available in DART
              nEdges = 30720 ;                                  available in DART
              maxEdges = 10 ;
              maxEdges2 = 20 ;
              nVertices = 20480 ;                               available in DART
              TWO = 2 ;
              THREE = 3 ;
              vertexDegree = 3 ;                                available in DART
              FIFTEEN = 15 ;
              TWENTYONE = 21 ;
              R3 = 3 ;
              nVertLevels = 41 ;                                available in DART
              nVertLevelsP1 = 42 ;                              available in DART
              nMonths = 12 ;
              nVertLevelsP2 = 43 ;
              nSoilLevels = 4 ;                                 available in DART
      variables:
              char xtime(Time, StrLen) ;                        available in DART
              double latCell(nCells) ;                          available in DART
              double lonCell(nCells) ;                          available in DART
              double latEdge(nEdges) ;
              double lonEdge(nEdges) ;
              int indexToEdgeID(nEdges) ;
              double latVertex(nVertices) ;
              double lonVertex(nVertices) ;
              int indexToVertexID(nVertices) ;
              int cellsOnEdge(nEdges, TWO) ;
              int nEdgesOnCell(nCells) ;
              int nEdgesOnEdge(nEdges) ;
              int edgesOnCell(nCells, maxEdges) ;
              int edgesOnEdge(nEdges, maxEdges2) ;
              double weightsOnEdge(nEdges, maxEdges2) ;
              double dvEdge(nEdges) ;
              double dcEdge(nEdges) ;
              double angleEdge(nEdges) ;
              double edgeNormalVectors(nEdges, R3) ;
              double cellTangentPlane(nEdges, TWO, R3) ;
              int cellsOnCell(nCells, maxEdges) ;
              int verticesOnCell(nCells, maxEdges) ;
              int verticesOnEdge(nEdges, TWO) ;
              int edgesOnVertex(nVertices, vertexDegree) ;
              int cellsOnVertex(nVertices, vertexDegree) ;      available in DART
              double kiteAreasOnVertex(nVertices, vertexDegree) ;
              double rainc(Time, nCells) ;
              double cuprec(Time, nCells) ;
              double cutop(Time, nCells) ;
              double cubot(Time, nCells) ;
              double relhum(Time, nCells, nVertLevels) ;
              double qsat(Time, nCells, nVertLevels) ;
              double graupelnc(Time, nCells) ;
              double snownc(Time, nCells) ;
              double rainnc(Time, nCells) ;
              double graupelncv(Time, nCells) ;
              double snowncv(Time, nCells) ;
              double rainncv(Time, nCells) ;
              double sr(Time, nCells) ;
              double surface_temperature(Time, nCells) ;
              double surface_pressure(Time, nCells) ;
              double coeffs_reconstruct(nCells, maxEdges, R3) ;
              double theta_base(Time, nCells, nVertLevels) ;
              double rho_base(Time, nCells, nVertLevels) ;
              double pressure_base(Time, nCells, nVertLevels) ;
              double exner_base(Time, nCells, nVertLevels) ;
              double exner(Time, nCells, nVertLevels) ;
              double h_divergence(Time, nCells, nVertLevels) ;
              double uReconstructMeridional(Time, nCells, nVertLevels) ;
              double uReconstructZonal(Time, nCells, nVertLevels) ;
              double uReconstructZ(Time, nCells, nVertLevels) ;
              double uReconstructY(Time, nCells, nVertLevels) ;
              double uReconstructX(Time, nCells, nVertLevels) ;
              double pv_cell(Time, nCells, nVertLevels) ;
              double pv_vertex(Time, nVertices, nVertLevels) ;
              double ke(Time, nCells, nVertLevels) ;
              double rho_edge(Time, nEdges, nVertLevels) ;
              double pv_edge(Time, nEdges, nVertLevels) ;
              double vorticity(Time, nVertices, nVertLevels) ;
              double divergence(Time, nCells, nVertLevels) ;
              double v(Time, nEdges, nVertLevels) ;
              double rh(Time, nCells, nVertLevels) ;
              double theta(Time, nCells, nVertLevels) ;
              double rho(Time, nCells, nVertLevels) ;
              double qv_init(nVertLevels) ;
              double t_init(nCells, nVertLevels) ;
              double u_init(nVertLevels) ;
              double pressure_p(Time, nCells, nVertLevels) ;
              double tend_theta(Time, nCells, nVertLevels) ;
              double tend_rho(Time, nCells, nVertLevels) ;
              double tend_w(Time, nCells, nVertLevelsP1) ;
              double tend_u(Time, nEdges, nVertLevels) ;
              double qv(Time, nCells, nVertLevels) ;
              double qc(Time, nCells, nVertLevels) ;
              double qr(Time, nCells, nVertLevels) ;
              double qi(Time, nCells, nVertLevels) ;
              double qs(Time, nCells, nVertLevels) ;
              double qg(Time, nCells, nVertLevels) ;
              double tend_qg(Time, nCells, nVertLevels) ;
              double tend_qs(Time, nCells, nVertLevels) ;
              double tend_qi(Time, nCells, nVertLevels) ;
              double tend_qr(Time, nCells, nVertLevels) ;
              double tend_qc(Time, nCells, nVertLevels) ;
              double tend_qv(Time, nCells, nVertLevels) ;
              double qnr(Time, nCells, nVertLevels) ;
              double qni(Time, nCells, nVertLevels) ;
              double tend_qnr(Time, nCells, nVertLevels) ;
              double tend_qni(Time, nCells, nVertLevels) ;

Namelist
--------

We adhere to the F90 standard of starting a namelist with an ampersand '&' and terminating with a slash '/' for all our
namelist input. Consider yourself forewarned that character strings that contain a '/' must be enclosed in quotes to
prevent them from prematurely terminating the namelist.

.. container:: namelist

   ::

      namelist /model_nml/  model_analysis_filename, &
                assimilation_period_days, assimilation_period_seconds, &
                model_perturbation_amplitude, output_state_vector, calendar, debug

.. container:: indent1

   This namelist is read in a file called ``input.nml``. This namelist provides control over the assimilation period for
   the model. All observations within (+/-) half of the assimilation period are assimilated. The assimilation period is
   the minimum amount of time the model can be advanced, and checks are performed to ensure that the assimilation window
   is a multiple of the model dynamical timestep. This also specifies the MPAS analysis file that will be read and/or
   written by the different program units.

   +---------------------------------------+---------------------------------------+---------------------------------------+
   | Contents                              | Type                                  | Description                           |
   +=======================================+=======================================+=======================================+
   | model_analysis_filename               | character(len=256)                    | Character string specifying the name  |
   |                                       | *[default: 'mpas_analysis.nc']*       | of the MPAS analysis file to be read  |
   |                                       |                                       | and/or written by the different       |
   |                                       |                                       | program units.                        |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | output_state_vector                   | logical *[default: .true.]*           | The switch to determine the form of   |
   |                                       |                                       | the state vector in the output netCDF |
   |                                       |                                       | files. If ``.true.`` the state vector |
   |                                       |                                       | will be output exactly as DART uses   |
   |                                       |                                       | it ... one long array. If             |
   |                                       |                                       | ``.false.``, the state vector is      |
   |                                       |                                       | parsed into prognostic variables and  |
   |                                       |                                       | output that way -- much easier to use |
   |                                       |                                       | with 'ncview', for example.           |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | assimilation_period_days              | integer *[default: 1]*                | The number of days to advance the     |
   |                                       |                                       | model for each assimilation.          |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | assimilation_period_seconds           | integer *[default: 0]*                | In addition to                        |
   |                                       |                                       | ``assimilation_period_days``, the     |
   |                                       |                                       | number of seconds to advance the      |
   |                                       |                                       | model for each assimilation.          |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | model_perturbation_amplitude          | real(r8) *[default: 0.2]*             | Reserved for future use.              |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | calendar                              | character(len=32)                     | Character string specifying the       |
   |                                       | *[default: 'Gregorian']*              | calendar being used by MPAS.          |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | debug                                 | integer *[default: 0]*                | The switch to specify the run-time    |
   |                                       |                                       | verbosity. ``0`` is as quiet as it    |
   |                                       |                                       | gets. ``> 1`` provides more run-time  |
   |                                       |                                       | messages. ``> 5`` provides ALL        |
   |                                       |                                       | run-time messages.                    |
   +---------------------------------------+---------------------------------------+---------------------------------------+

   .. rubric:: Example namelist
      :name: example-namelist

   ::

      &model_nml
         model_analysis_filename      = 'mpas_restart.nc';
         assimilation_period_days     = 0,
         assimilation_period_seconds  = 60,
         model_perturbation_amplitude = 0.2,
         output_state_vector          = .true.,
         calendar                     = 'Gregorian',
         debug                        = 0
         /

| 

.. container:: namelist

   ::

      namelist /mpas_vars_nml/ mpas_state_variables

.. container:: indent1

   This namelist is read from ``input.nml`` and contains the list of MPAS variables that make up the DART state vector.

   +---------------------------------------+---------------------------------------+---------------------------------------+
   | Contents                              | Type                                  | Description                           |
   +=======================================+=======================================+=======================================+
   | mpas_vars_nml                         | character(len=NF90_MAX_NAME)::        | The table that relates the GITM       |
   |                                       | dimension(160) *[default: see         | variables to use to build the DART    |
   |                                       | example]*                             | state vector, and the corresponding   |
   |                                       |                                       | DART kinds for those variables.       |
   +---------------------------------------+---------------------------------------+---------------------------------------+

   .. rubric:: Example
      :name: example
      :class: indent1

   The following mpas_vars_nml is just for demonstration purposes. You application will likely involve a different DART
   state vector.

   ::

      &mpas_vars_nml
         mpas_state_variables = 'theta',                 'QTY_POTENTIAL_TEMPERATURE',
                                'uReconstructZonal',     'QTY_U_WIND_COMPONENT',
                                'uReconstructMeridional','QTY_V_WIND_COMPONENT',
                                'w',                     'QTY_VERTICAL_VELOCITY',
                                'qv',                    'QTY_VAPOR_MIXING_RATIO',
                                'qc',                    'QTY_CLOUDWATER_MIXING_RATIO',
                                'qr',                    'QTY_RAINWATER_MIXING_RATIO',
                                'qi',                    'QTY_ICE_MIXING_RATIO',
                                'qs',                    'QTY_SNOW_MIXING_RATIO',
                                'qg',                    'QTY_GRAUPEL_MIXING_RATIO'
                                'surface_pressure',      'QTY_SURFACE_PRESSURE'
         /

   The variables are simply read from the MPAS analysis file and stored in the DART state vector such that all
   quantities for one variable are stored contiguously. Within each variable; they are stored vertically-contiguous for
   each horizontal location. From a storage standpoint, this would be equivalent to a Fortran variable dimensioned
   x(nVertical,nHorizontal,nVariables). The fastest-varying dimension is vertical, then horizontal, then variable ...
   naturally, the DART state vector is 1D. Each variable is also stored this way in the MPAS analysis file.

| 

Other modules used
------------------

::

   types_mod
   time_manager_mod
   threed_sphere/location_mod
   utilities_mod
   obs_kind_mod
   mpi_utilities_mod
   random_seq_mod

.. warning::

   DAReS staff began creating the MPAS_OCN interface to DART in preparation for the model's inclusion as the ocean
   component of the Community Earth System Model (CESM). The plans for including MPAS_OCN in CESM were abandoned and the
   Modular Ocean Model version 6 (MOM6) was included instead. Thus, the documentation on this page after this point
   describes an incomplete interface. Please contact DAReS staff by emailing dart@ucar.edu if you want to use DART with
   MPAS_OCN.

Public interfaces
-----------------

Only a select number of interfaces used are discussed here. Each module has its own discussion of their routines.

Required interface routines
~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Unique interface routines
~~~~~~~~~~~~~~~~~~~~~~~~~

======================= =========================
*use model_mod, only :* get_gridsize
\                       restart_file_to_sv
\                       sv_to_restart_file
\                       get_gitm_restart_filename
\                       get_base_time
\                       get_state_time
======================= =========================

+----------------------------+----------------------------------------------------------------------------------------+
| *use location_mod, only :* | `get_close_o                                                                           |
|                            | bs <../../assimilation_code/location/threed_sphere/location_mod.html#get_close_obs>`__ |
+----------------------------+----------------------------------------------------------------------------------------+

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

Interface routine descriptions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| 

.. container:: routine

   *model_size = get_model_size( )*
   ::

      integer :: get_model_size

.. container:: indent1

   Returns the length of the model state vector. Required.

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

   ``adv_1step`` is not used for the gitm model. Advancing the model is done through the ``advance_model`` script. This
   is a NULL_INTERFACE, provided only for compatibility with the DART requirements.

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

   ``get_state_meta_data`` returns metadata about a given element of the DART representation of the model state vector.
   Since the DART model state vector is a 1D array and the native model grid is multidimensional,
   ``get_state_meta_data`` returns information about the native model state vector representation. Things like the
   ``location``, or the type of the variable (for instance: temperature, u wind component, ...). The integer values used
   to indicate different variable types in ``var_type`` are themselves defined as public interfaces to model_mod if
   required.

   +--------------+------------------------------------------------------------------------------------------------------+
   | ``index_in`` | Index of state vector element about which information is requested.                                  |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``location`` | Returns the 3D location of the indexed state variable. The ``location_ type`` comes from             |
   |              | ``DART/assimilation_code/location/threed_sphere/location_mod.f90``. Note that the lat/lon are        |
   |              | specified in degrees by the user but are converted to radians internally.                            |
   +--------------+------------------------------------------------------------------------------------------------------+
   | *var_type*   | Returns the type of the indexed state variable as an optional argument. The type is one of the list  |
   |              | of supported observation types, found in the block of code starting                                  |
   |              | ``! Integer definitions for DART TYPES`` in                                                          |
   |              | ``DART/assimilation_code/modules/observations/obs_kind_mod.f90``                                     |
   +--------------+------------------------------------------------------------------------------------------------------+

   The list of supported variables in ``DART/assimilation_code/modules/observations/obs_kind_mod.f90`` is created by
   ``preprocess``.

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

   | Given a model state, ``model_interpolate`` returns the value of the desired observation type (which could be a
     state variable) that would be observed at the desired location. The interpolation method is either completely
     specified by the model, or uses some standard 2D or 3D scalar interpolation routines. Put another way,
     ``model_interpolate`` will apply the forward operator **H** to the model state to create an observation at the
     desired location.
   | If the interpolation is valid, ``istatus = 0``. In the case where the observation operator is not defined at the
     given location (e.g. the observation is below the lowest model level, above the top level, or 'dry'), interp_val is
     returned as 0.0 and istatus = 1.

   +---------------+-----------------------------------------------------------+
   | ``x``         | A model state vector.                                     |
   +---------------+-----------------------------------------------------------+
   | ``location``  | Location to which to interpolate.                         |
   +---------------+-----------------------------------------------------------+
   | ``itype``     | Integer indexing which type of observation is desired.    |
   +---------------+-----------------------------------------------------------+
   | ``obs_val``   | The interpolated value from the model.                    |
   +---------------+-----------------------------------------------------------+
   | ``istatus``   | Integer flag indicating the success of the interpolation. |
   |               | success == 0, failure == anything else                    |
   +---------------+-----------------------------------------------------------+

| 

.. container:: routine

   *var = get_model_time_step()*
   ::

      type(time_type) :: get_model_time_step

.. container:: indent1

   ``get_model_time_step`` returns the forecast length to be used as the "model base time step" in the filter. This is
   the minimum amount of time the model can be advanced by ``filter``. *This is also the assimilation window*. All
   observations within (+/-) one half of the forecast length are used for the assimilation. In the ``GITM`` case, this
   is set from the namelist values for
   ``input.nml``\ ``&model_nml:assimilation_period_days, assimilation_period_seconds``.

   ======= ============================
   ``var`` Smallest time step of model.
   ======= ============================

| 

.. container:: routine

   *call static_init_model()*

.. container:: indent1

   | ``static_init_model`` is called for runtime initialization of the model. The namelists are read to determine
     runtime configuration of the model, the grid coordinates, etc. There are no input arguments and no return values.
     The routine sets module-local private attributes that can then be queried by the public interface routines.
   | See the GITM documentation for all namelists in ``gitm_in`` . Be aware that DART reads the GITM ``&grid_nml``
     namelist to get the filenames for the horizontal and vertical grid information as well as the topography
     information.
   | The namelists (all mandatory) are:
   | ``input.nml``\ ``&model_mod_nml``,
   | ``gitm_in``\ ``&time_manager_nml``,
   | ``gitm_in``\ ``&io_nml``,
   | ``gitm_in``\ ``&init_ts_nml``,
   | ``gitm_in``\ ``&restart_nml``,
   | ``gitm_in``\ ``&domain_nml``, and
   | ``gitm_in``\ ``&grid_nml``.

| 

.. container:: routine

   *call end_model()*

.. container:: indent1

   ``end_model`` is used to clean up storage for the model, etc. when the model is no longer needed. There are no
   arguments and no return values. The grid variables are deallocated.

| 

.. container:: routine

   *call init_time(time)*
   ::

      type(time_type), intent(out) :: time

.. container:: indent1

   ``init_time`` returns the time at which the model will start if no input initial conditions are to be used. This is
   frequently used to spin-up models from rest, but is not meaningfully supported for the GITM model. The only time this
   routine would get called is if the ``input.nml``\ ``&perfect_model_obs_nml:start_from_restart`` is .false., which is
   not supported in the GITM model.

   ======== =====================================================================================================
   ``time`` the starting time for the model if no initial conditions are to be supplied. This is hardwired to 0.0
   ======== =====================================================================================================

| 

.. container:: routine

   *call init_conditions(x)*
   ::

      real(r8), dimension(:), intent(out) :: x

.. container:: indent1

   ``init_conditions`` returns default initial conditions for model; generally used for spinning up initial model
   states. For the GITM model it is just a stub because the initial state is always provided by the input files.

   ===== =============================================================
   ``x`` Initial conditions for state vector. This is hardwired to 0.0
   ===== =============================================================

| 

.. container:: routine

   *ierr = nc_write_model_atts(ncFileID)*
   ::

      integer             :: nc_write_model_atts
      integer, intent(in) :: ncFileID

.. container:: indent1

   ``nc_write_model_atts`` writes model-specific attributes to an opened netCDF file: In the GITM case, this includes
   information like the coordinate variables (the grid arrays: ULON, ULAT, TLON, TLAT, ZG, ZC, KMT, KMU), information
   from some of the namelists, and either the 1D state vector or the prognostic variables (SALT,TEMP,UVEL,VVEL,PSURF).
   All the required information (except for the netCDF file identifier) is obtained from the scope of the ``model_mod``
   module. Both the ``input.nml`` and ``gitm_in`` files are preserved in the netCDF file as variables ``inputnml`` and
   ``gitm_in``, respectively.

   ============ =========================================================
   ``ncFileID`` Integer file descriptor to previously-opened netCDF file.
   ``ierr``     Returns a 0 for successful completion.
   ============ =========================================================

   ``nc_write_model_atts`` is responsible for the model-specific attributes in the following DART-output netCDF files:
   ``true_state.nc``, ``preassim.nc``, and ``analysis.nc``.

| 

.. container:: routine

   *ierr = nc_write_model_vars(ncFileID, statevec, copyindex, timeindex)*
   ::

      integer,                intent(in) :: ncFileID
      real(r8), dimension(:), intent(in) :: statevec
      integer,                intent(in) :: copyindex
      integer,                intent(in) :: timeindex
      integer                            :: ierr

.. container:: indent1

   ``nc_write_model_vars`` writes a copy of the state variables to a NetCDF file. Multiple copies of the state for a
   given time are supported, allowing, for instance, a single file to include multiple ensemble estimates of the state.
   Whether the state vector is parsed into prognostic variables (SALT, TEMP, UVEL, VVEL, PSURF) or simply written as a
   1D array is controlled by ``input.nml``\ ``&model_mod_nml:output_state_vector``. If ``output_state_vector = .true.``
   the state vector is written as a 1D array (the simplest case, but hard to explore with the diagnostics). If
   ``output_state_vector = .false.`` the state vector is parsed into prognostic variables before being written.

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

   | Given a model state, ``pert_model_state`` produces a perturbed model state. This is used to generate ensemble
     initial conditions perturbed around some control trajectory state when one is preparing to spin-up ensembles. Since
     the DART state vector for the GITM model contains both 'wet' and 'dry' cells, it is imperative to provide an
     interface to perturb **just** the wet cells (``interf_provided == .true.``).
   | The magnitude of the perturbation is wholly determined by
     ``input.nml``\ ``&model_mod_nml:model_perturbation_amplitude`` and **utterly, completely fails**.
   | A more robust perturbation mechanism is needed. Until then, avoid using this routine by using your own ensemble of
     initial conditions. This is determined by setting ``input.nml``\ ``&filter_nml:start_from_restart = .false.``

   +---------------------+-----------------------------------------------------------------------------------------------+
   | ``state``           | State vector to be perturbed.                                                                 |
   +---------------------+-----------------------------------------------------------------------------------------------+
   | ``pert_state``      | The perturbed state vector.                                                                   |
   +---------------------+-----------------------------------------------------------------------------------------------+
   | ``interf_provided`` | Because of the 'wet/dry' issue discussed above, this is always ``.true.``, indicating a       |
   |                     | model-specific perturbation is available.                                                     |
   +---------------------+-----------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call get_close_maxdist_init(gc, maxdist)*
   ::

      type(get_close_type), intent(inout) :: gc
      real(r8),             intent(in)    :: maxdist

.. container:: indent1

   Pass-through to the 3-D sphere locations module. See
   `get_close_maxdist_init() <../../assimilation_code/location/threed_sphere/location_mod.html#get_close_maxdist_init>`__
   for the documentation of this subroutine.

| 

.. container:: routine

   *call get_close_obs_init(gc, num, obs)*
   ::

      type(get_close_type), intent(inout) :: gc
      integer,              intent(in)    :: num
      type(location_type),  intent(in)    :: obs(num)

.. container:: indent1

   Pass-through to the 3-D sphere locations module. See
   `get_close_obs_init() <../../assimilation_code/location/threed_sphere/location_mod.html#get_close_obs_init>`__ for
   the documentation of this subroutine.

| 

.. container:: routine

   *call get_close_obs(gc, base_obs_loc, base_obs_kind, obs, obs_kind, &
   num_close, close_ind [, dist])*
   ::

      type(get_close_type),              intent(in ) :: gc
      type(location_type),               intent(in ) :: base_obs_loc
      integer,                           intent(in ) :: base_obs_kind
      type(location_type), dimension(:), intent(in ) :: obs
      integer,             dimension(:), intent(in ) :: obs_kind
      integer,                           intent(out) :: num_close
      integer,             dimension(:), intent(out) :: close_ind
      real(r8), optional,  dimension(:), intent(out) :: dist

.. container:: indent1

   | Given a DART location (referred to as "base") and a set of locations, and a definition of 'close' - return a subset
     of locations that are 'close', as well as their distances to the DART location and their indices. This routine
     intentionally masks a routine of the same name in ``location_mod`` because we want to be able to discriminate
     against selecting 'dry land' locations.
   | Given a single location and a list of other locations, returns the indices of all the locations close to the single
     one along with the number of these and the distances for the close ones. The list of locations passed in via the
     ``obs`` argument must be identical to the list of ``obs`` passed into the most recent call to
     ``get_close_obs_init()``. If the list of locations of interest changes, ``get_close_obs_destroy()`` must be called
     and then the two initialization routines must be called before using ``get_close_obs()`` again.
   | For vertical distance computations, the general philosophy is to convert all vertical coordinates to a common
     coordinate. This coordinate type is defined in the namelist with the variable "vert_localization_coord".

   ================= =====================================================================================
   ``gc``            Structure to allow efficient identification of locations 'close' to a given location.
   ``base_obs_loc``  Single given location.
   ``base_obs_kind`` Kind of the single location.
   ``obs``           List of candidate locations.
   ``obs_kind``      Kind associated with candidate locations.
   ``num_close``     Number of locations close to the given location.
   ``close_ind``     Indices of those locations that are close.
   *dist*            Distance between given location and the close ones identified in close_ind.
   ================= =====================================================================================

| 

.. container:: routine

   *call ens_mean_for_model(ens_mean)*
   ::

      real(r8), dimension(:), intent(in) :: ens_mean

.. container:: indent1

   ``ens_mean_for_model`` normally saves a copy of the ensemble mean to module-local storage. This is a NULL_INTERFACE
   for the GITM model. At present there is no application which requires module-local storage of the ensemble mean. No
   storage is allocated.

   ============ ==========================================
   ``ens_mean`` State vector containing the ensemble mean.
   ============ ==========================================

| 

Unique interface routine descriptions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| 

.. container:: routine

   *call get_gridsize( num_x, num_y, num_z )*
   ::

      integer, intent(out) :: num_x, num_y, num_z

.. container:: indent1

   ``get_gridsize`` returns the dimensions of the compute domain. The horizontal gridsize is determined from
   ``gitm_restart.nc``.

   ========= ======================================
   ``num_x`` The number of longitudinal gridpoints.
   ``num_y`` The number of latitudinal gridpoints.
   ``num_z`` The number of vertical gridpoints.
   ========= ======================================

| 

.. container:: routine

   *call restart_file_to_sv(filename, state_vector, model_time)*
   ::

      character(len=*),       intent(in)    :: filename
      real(r8), dimension(:), intent(inout) :: state_vector
      type(time_type),        intent(out)   :: model_time

.. container:: indent1

   ``restart_file_to_sv`` Reads a GITM netCDF format restart file and packs the desired variables into a DART state
   vector. The desired variables are specified in the ``gitm_vars_nml`` namelist.

   ================ ======================================================================
   ``filename``     The name of the netCDF format GITM restart file.
   ``state_vector`` the 1D array containing the concatenated GITM variables.
   ``model_time``   the time of the model state. The last time in the netCDF restart file.
   ================ ======================================================================

| 

.. container:: routine

   *call sv_to_restart_file(state_vector, filename, statedate)*
   ::

      real(r8), dimension(:), intent(in) :: state_vector
      character(len=*),       intent(in) :: filename
      type(time_type),        intent(in) :: statedate

.. container:: indent1

   ``sv_to_restart_file`` updates the variables in the GITM restart file with values from the DART vector
   ``state_vector``. The last time in the file must match the ``statedate``.

   ================ ==================================================
   ``filename``     the netCDF-format GITM restart file to be updated.
   ``state_vector`` the 1D array containing the DART state vector.
   ``statedate``    the 'valid_time' of the DART state vector.
   ================ ==================================================

| 

.. container:: routine

   *call get_gitm_restart_filename( filename )*
   ::

      character(len=*), intent(out) :: filename

.. container:: indent1

   ``get_gitm_restart_filename`` returns the name of the gitm restart file - the filename itself is in private module
   storage.

   ============ ==================================
   ``filename`` The name of the GITM restart file.
   ============ ==================================

| 

.. container:: routine

   *time = get_base_time( filehandle )*
   ::

      integer,          intent(in) :: filehandle -OR-
      character(len=*), intent(in) :: filehandle
      type(time_type),  intent(out) :: time

.. container:: indent1

   ``get_base_time`` extracts the start time of the experiment as contained in the netCDF restart file. The file may be
   specified by either a character string or the integer netCDF fid.

| 

.. container:: routine

   *time = get_state_time( filehandle )*
   ::

      integer,          intent(in) :: filehandle -OR-
      character(len=*), intent(in) :: filehandle
      type(time_type),  intent(out) :: time

.. container:: indent1

   ``get_state_time`` extracts the time of the model state as contained in the netCDF restart file. In the case of
   multiple times in the file, the last time is the time returned. The file may be specified by either a character
   string or the integer netCDF fid.

| 

Files
-----

=========================== ===========================================================================
filename                    purpose
=========================== ===========================================================================
input.nml                   to read the model_mod namelist
gitm_vars.nml               to read the ``gitm_vars_nml`` namelist
gitm_restart.nc             provides grid dimensions, model state, and 'valid_time' of the model state
true_state.nc               the time-history of the "true" model state from an OSSE
preassim.nc                 the time-history of the model state before assimilation
analysis.nc                 the time-history of the model state after assimilation
dart_log.out [default name] the run-time diagnostic output
dart_log.nml [default name] the record of all the namelists actually USED - contains the default values
=========================== ===========================================================================

| 

References
----------

-  none

Private components
------------------

N/A
