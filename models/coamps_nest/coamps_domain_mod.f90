! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

!------------------------------
! MODULE:       coamps_domain_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Module containing the definition of a COAMPS domain data structure,
! which is an amalgamation of a grid (a map projection and how the domain
! connects to it), an arbitrary number of nests, and vertical structure.
!------------------------------ 
module coamps_domain_mod

    use coamps_intrinsic_mod, only : uvg2uv   

    use coamps_nest_mod,     only : coamps_nest,                 &
                                    nest_point,                  &
                                    initialize_nest,             &
                                    initialize_nest_latlon,      &
                                    dump_nest_info,              &
                                    decompose_nest,              &
                                    get_num_subnests,            &
                                    get_nest_id,                 &
                                    set_nest_id,                 &
                                    nest_point_to_coarse_point,  &
                                    coarse_point_to_nest_point,  &
                                    get_parent_nest_id,          &
                                    register_child_nest,         &
                                    register_parent_nest

    use coamps_vertical_mod, only : coamps_vertical,             &
                                    get_msigma,                  &
                                    get_wsigma,                  &
                                    get_dsigmaw,                 &
                                    get_num_levels,              &
                                    initialize_vertical,         &
                                    dump_vertical_info

    use coamps_map_mod,      only : coamps_grid,                 &
                                    initialize_grid,             &
                                    latlon_to_grid_point,        &
                                    grid_point_to_latlon,        &
                                    calc_grid_rotation,          &
                                    dump_grid_info

    use coamps_util_mod,     only : check_io_status,             &
                                    check_alloc_status,          &
                                    check_dealloc_status,        &
                                    generate_flat_file_name,     &
                                    read_datahd_file,            &
                                    DATAHD_LEN,                  &   
                                    trace_message,                    &
                                    DATAHD_NUM_NESTS 

    use location_mod,        only : get_location,                &
                                    location_type

    use utilities_mod,       only : do_output,                   &
                                    E_ERR,                       &
                                    E_WARN,                      &
                                    error_handler,               &
                                    get_unit,                    &
                                    register_module
    use types_mod,           only : r8

    implicit none

    private

    !------------------------------
    ! BEGIN PUBLIC INTERFACE
    !------------------------------

    ! COAMPS domain data structure and constructor
    public :: coamps_domain
    public :: initialize_domain

    !  ! Location conversion functions
    !  public :: gridpt_to_latlon
    public :: latlon_to_nest_point !!! FOR TESTING ONLY
    public :: nest_point_to_latlon
    public :: location_to_nest_point
    public :: grid_wind_to_earth_wind

    !  public :: grid_ij_to_vector_index
    public :: nest_ij_to_coarse_ij
    public :: coarse_ij_to_nest_ij
      
    !  ! Information about the grid size
    !  public :: get_grid_dims
    !  public :: get_grid_field_size
    !  public :: get_grid_num_levels
    !  public :: check_ij_within_grid
    !

    ! Information about vertical coordinates
    public :: get_domain_num_levels
    public :: get_domain_msigma
    public :: get_domain_wsigma
    public :: get_domain_dsigmaw

    ! Nest accessors
    public :: get_domain_nest
    public :: get_nest_count

    ! Domain decomposition tools 
    public :: get_num_subdomains
    public :: decompose_domain

    ! Debugging output
    public :: dump_domain_info

    !------------------------------
    ! END PUBLIC INTERFACE
    !------------------------------

    !------------------------------
    ! BEGIN EXTERNAL INTERFACES
    !------------------------------

    interface get_domain_msigma
        module procedure get_domain_msigma_array, get_domain_msigma_value
    end interface get_domain_msigma

    interface get_domain_wsigma
        module procedure get_domain_wsigma_array, get_domain_wsigma_value
    end interface get_domain_wsigma

    !------------------------------
    ! END EXTERNAL INTERFACES
    !------------------------------

    !------------------------------
    ! BEGIN TYPES AND CONSTANTS 
    !------------------------------

    type :: coamps_domain
        private

        type(coamps_grid)                        :: static_grid

        integer                                  :: nest_count
        type(coamps_nest), dimension(:), pointer :: nests

        type(coamps_vertical)                    :: vertical

        logical :: is_initialized = .false.
    end type coamps_domain

    integer, parameter :: COARSE_NEST = 1

    !------------------------------
    ! END TYPES AND CONSTANTS 
    !------------------------------

    !------------------------------
    ! BEGIN MODULE VARIABLES
    !------------------------------

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

    logical, save :: module_initialized = .false.

    !------------------------------
    ! END MODULE VARIABLES
    !------------------------------

contains
    !------------------------------
    ! BEGIN PUBLIC ROUTINES
    !------------------------------

    ! initialize_domain
    ! -----------------
    ! Constructor for a COAMPS domain possibly consisting of several nests
    !  PARAMETERS
    !   IN  dtg               COAMPS date-time-group (for filenames)
    !   OUT domain            Domain to initialize
    subroutine initialize_domain(dtg, domain)
        character(len=10),      intent(in) :: dtg
        type(coamps_domain), intent(inout) :: domain

        real(kind=r8), dimension(DATAHD_LEN) :: coamps_datahd

        if (.not. module_initialized) then
            call register_module(source, revision, revdate)
            module_initialized = .true.
        end if

        if (domain%is_initialized) return

        call read_datahd_file(dtg, coamps_datahd)

        call initialize_grid(coamps_datahd, domain%static_grid)

        call initialize_nests(dtg, coamps_datahd, domain)

        call initialize_vertical(coamps_datahd, domain%vertical)

        domain%is_initialized = .true.
    end subroutine initialize_domain

    ! get_domain_msigma_array
    ! -----------------------
    ! Returns the mass sigma levels on the COAMPS domain
    !  PARAMETERS
    !   IN  domain            COAMPS domain to pull from
    function get_domain_msigma_array(domain)
        type(coamps_domain),         intent(in)  :: domain
        real(kind=r8), dimension(:), pointer     :: get_domain_msigma_array    

        get_domain_msigma_array => get_msigma(domain%vertical)
    end function get_domain_msigma_array

    ! get_domain_msigma_value
    ! -----------------------
    ! Returns the mass sigma level on the COAMPS domain at a specified
    ! index
    !  PARAMETERS
    !   IN  domain          COAMPS domain to pull from
    !   IN  sigma_index     Sigma level index to isolate
    function get_domain_msigma_value(domain, sigma_index)
        type(coamps_domain), intent(in)  :: domain
        integer,             intent(in)  :: sigma_index
        real(kind=r8)                    :: get_domain_msigma_value

        get_domain_msigma_value = get_msigma(domain%vertical, sigma_index)
    end function get_domain_msigma_value

    ! get_domain_wsigma_array
    ! -----------------------
    ! Returns the mass sigma levels on the COAMPS domain
    !  PARAMETERS
    !   IN  domain            COAMPS domain to pull from
    function get_domain_wsigma_array(domain)
        type(coamps_domain),         intent(in)  :: domain
        real(kind=r8), dimension(:), pointer     :: get_domain_wsigma_array    

        get_domain_wsigma_array => get_wsigma(domain%vertical)
    end function get_domain_wsigma_array

    ! get_domain_wsigma_value
    ! -----------------------
    ! Returns the mass sigma level on the COAMPS domain at a specified
    ! index
    !  PARAMETERS
    !   IN  domain          COAMPS domain to pull from
    !   IN  sigma_index     Sigma level index to isolate
    function get_domain_wsigma_value(domain, sigma_index)
        type(coamps_domain), intent(in)  :: domain
        integer,             intent(in)  :: sigma_index
        real(kind=r8)                    :: get_domain_wsigma_value

        get_domain_wsigma_value = get_wsigma(domain%vertical, sigma_index)
    end function get_domain_wsigma_value

    ! get_domain_dsigmaw
    ! -----------------
    ! Returns the mass sigma levels on the COAMPS domain
    !  PARAMETERS
    !   IN  domain            COAMPS domain to pull from
    function get_domain_dsigmaw(domain)
        type(coamps_domain),         intent(in)  :: domain
        real(kind=r8), dimension(:), pointer     :: get_domain_dsigmaw

        get_domain_dsigmaw => get_dsigmaw(domain%vertical)
    end function get_domain_dsigmaw

    ! get_nest_count
    ! --------------
    ! Returns the number of nests for this domain
    !  PARAMETERS
    !   IN  domain          coamps_domain to query
    function get_nest_count(domain)
        type(coamps_domain), intent(in)  :: domain
        integer                          :: get_nest_count

        get_nest_count = domain%nest_count
    end function get_nest_count

    ! get_num_subdomains
    ! ------------------
    ! Return the number of subdomains computed by decompose_domain 
    ! Will throw an error if decompose_domain has not yet been run
    !  PARAMETERS
    !   IN  domain          coamps_domain to query
    function get_num_subdomains(domain)
        type(coamps_domain), intent(in)  :: domain
        integer                          :: get_num_subdomains

        ! The total number of subdomains is independent of the nest
        get_num_subdomains = get_num_subnests(domain%nests(COARSE_NEST))
    end function get_num_subdomains

    ! latlon_to_nest_point
    ! --------------------
    ! Given a COAMPS domain, convert a given (lat, lon) representation of
    ! a point to a nest_point (i, j, nest_num) representation.
    !  PARAMETERS
    !   IN  domain            coamps_domain to base conversion on
    !   IN  lat               point's latitude
    !   IN  lon               point's longitude
    !   OUT nest_pt           highest-resolution nest point at this lat/lon
    !   OUT within_domain     (Optional) True if the converted point is in
    !                         the boundaries of the domain
    subroutine latlon_to_nest_point(domain, lat, lon, nest_pt, within_domain) 
        type(coamps_domain), intent(in)  :: domain
        real(kind=r8),       intent(in)  :: lat
        real(kind=r8),       intent(in)  :: lon
        type(nest_point),    intent(out) :: nest_pt
        logical, optional,   intent(out) :: within_domain

        type(nest_point) :: coarse_pt
        real(kind=r8)    :: coarse_i
        real(kind=r8)    :: coarse_j

        logical          :: in_domain

        nullify(nest_pt%nest)

        ! COAMPS latlon conversion only works with the top-level grid
        call latlon_to_grid_point(domain%static_grid, lat, lon, &
                                  coarse_i, coarse_j)

        coarse_pt%ii   =  coarse_i
        coarse_pt%jj   =  coarse_j
        coarse_pt%nest => domain%nests(COARSE_NEST)

        call coarse_point_to_nest_point(coarse_pt, nest_pt, in_domain)

        if (present(within_domain)) within_domain = in_domain
    end subroutine latlon_to_nest_point

    ! nest_point_to_latlon
    ! --------------------
    ! Given a COAMPS domain, convert a given (i, j, nest) representation of
    ! a point to a nest_point (lat, lon) representation.
    !  PARAMETERS
    !   IN  domain            coamps_domain to base conversion on
    !   IN  nest_pt           nest point to convert
    !   IN  lat               point's latitude
    !   IN  lon               point's longitude
    subroutine nest_point_to_latlon(domain, nest_pt, lat, lon)
        type(coamps_domain), intent(in)  :: domain
        type(nest_point),    intent(in)  :: nest_pt
        real(kind=r8),       intent(out) :: lat
        real(kind=r8),       intent(out) :: lon

        type(nest_point) :: coarse_pt

        call nest_point_to_coarse_point(nest_pt, coarse_pt)

        ! COAMPS latlon conversion only works with the top-level grid
        call grid_point_to_latlon(domain%static_grid, coarse_pt%ii, coarse_pt%jj, &
                                  lat, lon)

    end subroutine nest_point_to_latlon

    ! location_to_nest_point
    ! ----------------------
    ! Converts a DART location to a point within a COAMPS domain.  The
    ! domain point consists of a vertical coordinate as well as a nest
    ! point (which supplies a nest number as well as an i/j coordinate)
    !  PARAMETERS
    !   IN  domain            coamps_domain for lat/lon conversion
    !   IN  loc               DART location structure to unpack
    !   OUT nest_pt           COAMPS nest point with horizontal location
    !   OUT vert_loc          point's vertical location
    !   OUT is_within_domain  (Optional) True if the converted point is
    !                         within the domain boundaries
    subroutine location_to_nest_point(domain, loc, nest_pt, vert_loc, &
                                      is_within_domain)
        type(coamps_domain), intent(in)  :: domain
        type(location_type), intent(in)  :: loc
        type(nest_point),    intent(out) :: nest_pt
        real(kind=r8),       intent(out) :: vert_loc
        logical, optional,   intent(out) :: is_within_domain

        integer, parameter :: DART_LOC_LON  = 1
        integer, parameter :: DART_LOC_LAT  = 2
        integer, parameter :: DART_LOC_VERT = 3

        ! The array returned by get_location is (lon, lat, vertical)
        real(kind=r8), dimension(3) :: loc_array

        logical :: in_domain

        loc_array = get_location(loc)

        ! All nests share the same vertical coordinates
        vert_loc = loc_array(DART_LOC_VERT)

        call latlon_to_nest_point(domain, loc_array(DART_LOC_LAT), &
                                  loc_array(DART_LOC_LON), nest_pt,&
                                  in_domain) 

        if (present(is_within_domain)) is_within_domain = in_domain
    end subroutine location_to_nest_point

    ! decompose_domain
    ! ----------------
    ! Decomposes each nest in the domain given the number of I/O
    ! processors, the number of processors in x and y, and the number
    ! of boundary points (nbnam in the COAMPS namelist).  Note that for
    ! N nests, the processors should be arrays of length N
    !  PARAMETERS
    ! INOUT domain            COAMPS domain to decompose
    !   IN  num_io_procs      number of I/O processors (0 or 1)
    !   IN  dom_proc_x        number of processors in the x direction
    !   IN  dom_proc_y        number of processors in the y direction
    !   IN  nbound            number of halo boundary points
    subroutine decompose_domain(domain, num_io_procs, dom_proc_x, &
                                dom_proc_y, n_bound)
        type(coamps_domain),   intent(inout) :: domain
        integer,               intent(in)    :: num_io_procs
        integer, dimension(:), intent(in)    :: dom_proc_x
        integer, dimension(:), intent(in)    :: dom_proc_y
        integer,               intent(in)    :: n_bound

        integer :: cur_nest

        do cur_nest = 1, domain%nest_count
        call decompose_nest(domain%nests(cur_nest), num_io_procs,       &
                            dom_proc_x(cur_nest), dom_proc_y(cur_nest), &
                            n_bound)
        end do
    end subroutine decompose_domain

    ! nest_ij_to_coarse_ij
    ! --------------------
    ! Converts an i/j point on a COAMPS nest to the corresponding i/j point
    ! on the coarsest nest
    !  PARAMETERS
    !   IN  domain            COAMPS domain that contains the nests
    !   IN  cur_nest          The nest the current i/j point is defined on
    !   IN  nest_i            x-location on the nest
    !   IN  nest_j            y-location on the nest
    !   OUT coarse_i          x-location on the coarse grid
    !   OUT coarse_j          y-location on the coarse grid
    subroutine nest_ij_to_coarse_ij(domain, cur_nest, nest_i, nest_j, &
                                    coarse_i, coarse_j)
        type(coamps_domain), intent(in)  :: domain
        integer,             intent(in)  :: cur_nest
        real(kind=r8),       intent(in)  :: nest_i
        real(kind=r8),       intent(in)  :: nest_j
        real(kind=r8),       intent(out) :: coarse_i
        real(kind=r8),       intent(out) :: coarse_j

        type(nest_point) :: nest_pt
        type(nest_point) :: coarse_pt

        nest_pt%nest => domain%nests(cur_nest)
        nest_pt%ii   =  nest_i
        nest_pt%jj   =  nest_j

        call nest_point_to_coarse_point(nest_pt, coarse_pt)

        ! Nest is implicitly COARSE_NEST
        coarse_i = coarse_pt%ii
        coarse_j = coarse_pt%jj
    end subroutine nest_ij_to_coarse_ij

    ! coarse_ij_to_nest_ij
    ! --------------------
    ! Based on an i/j point in the coarse mesh, finds the finest-resolution 
    ! nest containing that point and returns the i/j location on that nest as
    ! well as the nest number
    !  PARAMETERS
    !   IN  domain            COAMPS domain that contains the nests
    !   IN  coarse_i          x-location on coarse nest
    !   IN  coarse_j          y-location on coarse nest
    !   OUT nest_i            x-location on finest mesh
    !   OUT nest_j            y-location on finest mesh
    !   OUT nest_num          finest-resolution mesh available for this point
    subroutine coarse_ij_to_nest_ij(domain, coarse_i, coarse_j, nest_i, &
                                    nest_j, nest_num)
        type(coamps_domain), intent(in)  :: domain
        real(kind=r8),       intent(in)  :: coarse_i
        real(kind=r8),       intent(in)  :: coarse_j 
        real(kind=r8),       intent(out) :: nest_i
        real(kind=r8),       intent(out) :: nest_j
        integer,             intent(out) :: nest_num

        type(nest_point) :: coarse_pt
        type(nest_point) :: nest_pt

        coarse_pt%ii   = coarse_i
        coarse_pt%jj   = coarse_j
        coarse_pt%nest => domain%nests(COARSE_NEST)

        call coarse_point_to_nest_point(coarse_pt, nest_pt)

        nest_i   = nest_pt%ii
        nest_j   = nest_pt%jj
        nest_num = get_nest_id(nest_pt%nest)
    end subroutine coarse_ij_to_nest_ij

    ! get_domain_nest
    ! --------
    ! Given an index, returns that nest
    !  PARAMETERS
    !   IN  domain          COAMPS domain to get a nest from
    !   IN  nest_index      Number of the nest to return
    function get_domain_nest(domain, nest_index)
        type(coamps_domain),         intent(in)  :: domain
        integer,                     intent(in)  :: nest_index
        type(coamps_nest),   target              :: get_domain_nest

        if ((nest_index .lt. COARSE_NEST)        .or. &
            (nest_index .gt. domain%nest_count)) then
            call error_handler(E_ERR, 'get_domain_nest', 'Nest index is out ' // &
                               'of range', source, revision, revdate)
        end if

        get_domain_nest = domain%nests(nest_index)
    end function get_domain_nest 

    ! get_domain_num_levels
    ! ---------------------
    ! Return the number of (mass) sigma levels defined on the domain
    !  PARAMETERS
    !   IN  domain          COAMPS domain to query
    function get_domain_num_levels(domain)
        type(coamps_domain), intent(in)  :: domain
        integer                          :: get_domain_num_levels

        get_domain_num_levels = get_num_levels(domain%vertical)
    end function get_domain_num_levels

    ! dump_domain_info
    ! ----------------
    ! Dumps a human-readable summary of the supplied domain
    !  PARAMETERS
    !   IN  domain            COAMPS domain to print out
    !   IN  print_vertical    (Optional) True if this should print sigma data
    subroutine dump_domain_info(domain, print_vertical)
    type(coamps_domain), intent(in) :: domain
    logical, optional,   intent(in) :: print_vertical

    integer :: cur_nest

    if (do_output()) then
        write (*,*) "***************************************"
        write (*,*) "****** COAMPS DOMAIN INFORMATION ******"
        write (*,*) "***************************************"

        write (*,*)
        call dump_grid_info(domain%static_grid)
        write (*,*)

        write (*,*) "**** COAMPS NESTS ****"
        write (*,*) "Number of nests:          ", domain%nest_count

        do cur_nest = 1, domain%nest_count
        !call dump_nest_info(domain%nests(cur_nest), .true.)
          call dump_nest_info(domain%nests(cur_nest))
        end do

        if (present(print_vertical)) then
            if (print_vertical) then
                call dump_vertical_info(domain%vertical)
            end if
        end if

    end if
    end subroutine dump_domain_info

    ! grid_wind_to_earth_wind
    ! ----------------
    ! Rotates the grid relative winds to a true wind
    !  PARAMETERS
    !   IN  grid              coamps_grid structure
    !   IN  dart_loc          dart_location
    !   INOUT v_wind          v wind component
    subroutine grid_wind_to_earth_wind(u_wind, v_wind, domain, dart_loc)
      real(kind=r8),       intent(inout) :: u_wind
      real(kind=r8),       intent(inout) :: v_wind
      type(coamps_domain), intent(in)    :: domain
      type(location_type), intent(in)    :: dart_loc

      integer, parameter          :: SINGLE_POINT  = 1
      integer, parameter          :: DART_LOC_LON  = 1
      integer, parameter          :: DART_LOC_LAT  = 2

      real(kind=r8), dimension(3) :: location
      real(kind=r8)               :: longitude
      real(kind=r8)               :: grid_rotation(SINGLE_POINT)

      real(kind=r8)               :: u_tmp(SINGLE_POINT, SINGLE_POINT)
      real(kind=r8)               :: v_tmp(SINGLE_POINT, SINGLE_POINT)
      character(len=90)           :: uvstr

      u_tmp(SINGLE_POINT, SINGLE_POINT) = u_wind
      v_tmp(SINGLE_POINT, SINGLE_POINT) = v_wind

      location   = get_location(dart_loc)
      longitude  = location(DART_LOC_LON)
      longitude  = location(DART_LOC_LON)
      if(longitude<0.0_r8) longitude=longitude+360.0_r8

      grid_rotation(SINGLE_POINT) = calc_grid_rotation(domain%static_grid, longitude)
      call uvg2uv(u_tmp, v_tmp, SINGLE_POINT, SINGLE_POINT, grid_rotation, &
                  u_tmp, v_tmp)

      u_wind     = u_tmp(SINGLE_POINT, SINGLE_POINT)
      v_wind     = v_tmp(SINGLE_POINT, SINGLE_POINT)
    end subroutine grid_wind_to_earth_wind

!  ! gridpt_to_latlon
!  ! ----------------
!  ! Given a COAMPS grid structure, convert an (i,j) representation
!  ! of a point on that grid to a (lat, lon) representation.
!  ! Note that this routine will happily accept and convert (i,j)
!  ! points that are beyond the bounds of the domain.  For this 
!  ! direction of conversion, that's not too big a problem.
!  !  PARAMETERS
!  !   IN  grid              coamps_grid to base conversion on
!  !   IN  ii                point's x-direction index
!  !   IN  jj                point's y-direction index
!  !   OUT lat               point's latitude
!  !   OUT lon               point's longitude
!  subroutine gridpt_to_latlon(grid, ii, jj, lat, lon) 
!    type(coamps_grid), intent(in) :: grid
!    real(kind=r8), intent(in)     :: ii
!    real(kind=r8), intent(in)     :: jj
!    real(kind=r8), intent(out)    :: lat
!    real(kind=r8), intent(out)    :: lon
!
!    ! Work-around for not being able to pass in scalars as arrays
!    ! of length 1
!    real(kind=r8), dimension(SINGLE_POINT) :: ii_array,  jj_array
!    real(kind=r8), dimension(SINGLE_POINT) :: lat_array, lon_array
!
!    ii_array(SINGLE_POINT) = ii
!    jj_array(SINGLE_POINT) = jj
!
!    call ij2ll(int(grid%map_proj), grid%reflat,    grid%reflon,  &
!               int(grid%iref),     int(grid%jref), grid%stdlat1, &
!               grid%stdlat2,       grid%stdlon,    grid%delta_x, &
!               grid%delta_y,       ii_array,       jj_array,     &
!               SINGLE_POINT,       lat_array,      lon_array)  
!
!    lat = lat_array(SINGLE_POINT)
!    lon = lon_array(SINGLE_POINT)
!  end subroutine gridpt_to_latlon
!
!
!  ! grid_ij_to_vector_index
!  ! -----------------------
!  ! Given a COAMPS grid, converts an horizontal (ii,jj) location 
!  ! to an index in the vector representation
!  !  PARAMETERS
!  !   IN  grid              coamps_grid to use for the conversion
!  !   IN  ii                zonal index component
!  !   IN  jj                meridional index component
!  !   OUT index             position that the point (ii,jj) has in
!  !                         the field re-arranged in vector form
!  subroutine grid_ij_to_vector_index(grid, ii, jj, index)
!    type(coamps_grid), intent(in) :: grid
!    integer, intent(in)           :: ii
!    integer, intent(in)           :: jj
!    integer, intent(out)          :: index
!
!    index = (jj - 1) * grid%pts_x + ii
!  end subroutine grid_ij_to_vector_index
!
!  ! get_grid_dims
!  ! -------------
!  ! Given a COAMPS grid structure, calculate the i (i.e. x) and 
!  ! j (i.e. y) dimensions
!  !  PARAMETERS
!  !   IN  grid              coamps_grid to get limits from
!  !   OUT grid_i            number of points on the i dimension
!  !   OUT grid_j            number of points on the j dimension
!  subroutine get_grid_dims(grid, grid_i, grid_j)
!    type(coamps_grid), intent(in) :: grid
!    integer, intent(out)          :: grid_i
!    integer, intent(out)          :: grid_j
!
!    grid_i = int(grid%pts_x)
!    grid_j = int(grid%pts_y)
!  end subroutine get_grid_dims
!
!  ! get_grid_field_size
!  ! -------------------
!  ! Given a coamps grid structure, calculate the field size (the
!  ! total number of points in the i and j directions)
!  !  PARAMETERS
!  !   IN  grid              coamps_grid to get size from
!  !   OUT field_size        number of grid points on a single level
!  subroutine get_grid_field_size(grid, field_size)
!    type(coamps_grid), intent(in) :: grid
!    integer, intent(out)          :: field_size
!
!    field_size = int(grid%pts_x) * int(grid%pts_y)
!  end subroutine get_grid_field_size
!
!  ! get_grid_num_levels
!  ! -------------------
!  ! Returns the number of mass sigma levels on a given COAMPS grid
!  !  PARAMETERS
!  !   IN  grid              coamps_grid to get levels from
!  !   OUT num_sig_lvls      number of sigma levels on the grid
!  subroutine get_grid_num_levels(grid, num_sig_lvls)
!    type(coamps_grid), intent(in) :: grid
!    integer, intent(out)          :: num_sig_lvls
!    
!    num_sig_lvls = grid%sigm_lvls
!  end subroutine get_grid_num_levels
!
!  ! check_ij_within_grid
!  ! --------------------
!  ! Calculates if a given (i,j) index is within a given COAMPS grid
!  !  PARAMETERS
!  !   IN  grid              coamps_grid to use for checking
!  !   IN  ii                zonal index
!  !   IN  jj                meridional index
!  !   OUT is_within_grid    True if the (i,j) location listed is
!  !                         inside the given grid
!  subroutine check_ij_within_grid(grid, ii, jj, is_within_grid)
!    type(coamps_grid), intent(in) :: grid
!    integer, intent(in)           :: ii
!    integer, intent(in)           :: jj
!    logical, intent(out)          :: is_within_grid
!  
!    is_within_grid = .true.
!    if ( (ii < 1) .or. (ii > grid%pts_x) .or. &
!         (jj < 1) .or. (jj > grid%pts_y) ) then
!       is_within_grid = .false.
!    end if
!  end subroutine check_ij_within_grid
!

  !------------------------------
  ! END PUBLIC ROUTINES
  !------------------------------

  !------------------------------
  ! BEGIN PRIVATE ROUTINES
  !------------------------------

  ! initialize_nests
  ! ----------------
  ! Populate the nest portions of a COAMPS domain based on the given
  ! datahd record
  !  PARAMETERS
  !   IN  dtg               date-time-group of forecast
  !   IN  datahd            datahd record to source
  ! INOUT domain            COAMPS domain to populate
  subroutine initialize_nests(dtg, datahd, domain)
    character(len=10),           intent(in)    :: dtg
    real(kind=r8), dimension(:), intent(in)    :: datahd
    type(coamps_domain),         intent(inout) :: domain 

    character(len=*), parameter :: routine = 'initialize_nests'  
    integer                     :: alloc_status

    integer :: cur_nest_id, parent_nest_id

    domain%nest_count = int(datahd(DATAHD_NUM_NESTS))
    allocate(domain%nests(domain%nest_count), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'domain%nests')

    do cur_nest_id = 1, domain%nest_count
        call set_nest_id(domain%nests(cur_nest_id), cur_nest_id)

        call initialize_nest(domain%nests(cur_nest_id), dtg, datahd)    

        ! Interface this nest with the other nests
        parent_nest_id = get_parent_nest_id(domain%nests(cur_nest_id))
        call register_parent_nest(parent_nest = domain%nests(parent_nest_id),&
                                  child_nest  = domain%nests(cur_nest_id))
        call register_child_nest(child_nest  = domain%nests(cur_nest_id),  &
                                 parent_nest = domain%nests(parent_nest_id))
    end do

    do cur_nest_id = 1, domain%nest_count
        call initialize_nest_latlon(domain%nests(cur_nest_id), domain%static_grid)
    end do

  end subroutine initialize_nests

  !------------------------------
  ! END PRIVATE ROUTINES
  !------------------------------

end module coamps_domain_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
