! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

!------------------------------
! MODULE:       coamps_map_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Module containing data structures and routines for dealing
! with the map component of a coamps domain
!------------------------------ 

module coamps_map_mod

    use coamps_intrinsic_mod, only : ij2ll, ll2ij
    use types_mod, only : r8

    implicit none

    private

    !------------------------------
    ! BEGIN PUBLIC INTERFACE
    !------------------------------

    public :: coamps_grid
    public :: initialize_grid

    public :: grid_point_to_latlon
    public :: latlon_to_grid_point
    public :: calc_grid_rotation

    public :: dump_grid_info

    !------------------------------
    ! END PUBLIC INTERFACE
    !------------------------------

    !------------------------------
    ! BEGIN EXTERNAL INTERFACES
    !------------------------------
    ! [none]
    !------------------------------
    ! END EXTERNAL INTERFACES
    !------------------------------

    !------------------------------
    ! BEGIN TYPES AND CONSTANTS 
    !------------------------------

    type :: coamps_grid
        private

        integer       :: map_proj
        real(kind=r8) :: stdlat1
        real(kind=r8) :: stdlat2
        real(kind=r8) :: stdlon
        real(kind=r8) :: reflat
        real(kind=r8) :: reflon

        real(kind=r8) :: map_constant

        ! These are from the top-level grid 
        integer       :: iref
        integer       :: jref
        real(kind=r8) :: coarse_delta_x
        real(kind=r8) :: coarse_delta_y
    end type coamps_grid

    ! Static grid information in the COAMPS datahd file - note that
    ! the reference_i and reference_j points technically come from
    ! the coarse grid and are not stored with the rest of the static
    ! grid information
    integer, parameter :: DATAHD_MAP_PROJECTION = 3
    integer, parameter :: DATAHD_STANDARD_LAT_1 = 4
    integer, parameter :: DATAHD_STANDARD_LAT_2 = 5
    integer, parameter :: DATAHD_STANDARD_LON   = 6
    integer, parameter :: DATAHD_REFERENCE_LAT  = 7
    integer, parameter :: DATAHD_REFERENCE_LON  = 8
    integer, parameter :: DATAHD_COARSE_DELTA_X = 9
    integer, parameter :: DATAHD_COARSE_DELTA_Y = 10
    integer, parameter :: DATAHD_REFERENCE_I    = 34
    integer, parameter :: DATAHD_REFERENCE_J    = 35

    ! POSSIBLE COAMPS PROJECTIONS
    integer, parameter :: PROJ_MERCATOR                = 1
    integer, parameter :: PROJ_LAMBERT_CONIC_CONFORMAL = 2
    integer, parameter :: PROJ_POLAR_STEREOGRAPHIC     = 3
    integer, parameter :: PROJ_NONE                    = 4

    ! Easy recasting of scalar -> dimension(1) array
    integer, parameter :: SINGLE_POINT = 1

    real(kind=r8), parameter :: pi         = 3.141592741012573_r8
    real(kind=r8), parameter :: DEG_TO_RAD = pi/180.0_r8

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

    !------------------------------
    ! END MODULE VARIABLES
    !------------------------------

contains

    !------------------------------
    ! BEGIN PUBLIC ROUTINES
    !------------------------------

    ! initialize_grid
    ! ----------------------
    ! Populate the static grid portion of a COAMPS domain based on the
    ! given datahd record
    !  PARAMETERS
    !   IN  datahd            datahd record to source
    ! INOUT grid              COAMPS grid to populate
    subroutine initialize_grid(datahd, grid)
        real(kind=r8), dimension(:),  intent(in)    :: datahd
        type(coamps_grid),            intent(inout) :: grid

        grid%map_proj       = int(datahd(DATAHD_MAP_PROJECTION))
        grid%stdlat1        = datahd(DATAHD_STANDARD_LAT_1)
        grid%stdlat2        = datahd(DATAHD_STANDARD_LAT_2)
        grid%stdlon         = datahd(DATAHD_STANDARD_LON)
        grid%reflat         = datahd(DATAHD_REFERENCE_LAT)
        grid%reflon         = datahd(DATAHD_REFERENCE_LON) 
        grid%coarse_delta_x = datahd(DATAHD_COARSE_DELTA_X)
        grid%coarse_delta_y = datahd(DATAHD_COARSE_DELTA_Y)

        if(grid%stdlon < 0.0_r8) grid%stdlon = grid%stdlon + 360.0_r8
        if(grid%reflon < 0.0_r8) grid%reflon = grid%reflon + 360.0_r8

        ! Implicitly taken from grid 1
        grid%iref           = int(datahd(DATAHD_REFERENCE_I))
        grid%jref           = int(datahd(DATAHD_REFERENCE_J))

        grid%map_constant   = calc_map_constant(grid)
    end subroutine initialize_grid

    ! calc_grid_rotation
    ! ----------------
    ! Given a longitude and grid structure, calculates the 
    ! rotation angle between grid north and true north.
    !  PARAMETERS
    !   IN  grid              coamps_grid structure
    !   IN  longitude         longitude of point
    function calc_grid_rotation(grid, longitude) result(grid_rotation)
      type(coamps_grid),   intent(in)  :: grid
      real(kind=r8),       intent(in)  :: longitude
      real(kind=r8)                    :: grid_rotation
      real(kind=r8)                    :: long_loc
    
      if(longitude < 0 ) then
        long_loc = longitude+360.0_r8
      else
        long_loc = longitude
      endif

      select case(int(grid%map_proj))
        case(PROJ_LAMBERT_CONIC_CONFORMAL, PROJ_POLAR_STEREOGRAPHIC)

          grid_rotation=(grid%stdlon-long_loc)*grid%map_constant
          !case 1: standard longitude is east of Meridian, (between 0 and 90E)
          !and the grid longitude is west of Meridian (between 90W and Meridian)
          if (grid%stdlon >=   0.0_r8 .and. grid%stdlon <= 90.0_r8 .and.  &
              long_loc    >= 270.0_r8 .and.   long_loc  < 360.0_r8) then
                  grid_rotation=(grid%stdlon+360.0_r8-long_loc)*grid%map_constant
          !case 2: standard longitude is west of Meridian, (between 90W and Meridian)
          !and the grid longitude is east of Meridian (between 0 and 90E)
          elseif (grid%stdlon >= 270.0_r8 .and. grid%stdlon <  360.0_r8 .and. &
                  long_loc    >=   0.0_r8 .and. long_loc    <= 90.0_r8) then
              grid_rotation=(grid%stdlon-long_loc-360.0_r8)*grid%map_constant
          else
              grid_rotation=(grid%stdlon-long_loc)*grid%map_constant
          endif
        case default
          grid_rotation=0.0_r8
      end select
    end function calc_grid_rotation

    ! latlon_to_grid_point
    ! --------------------
    ! Converts a lat/lon point to its representation on the COAMPS grid
    !  PARAMETERS
    !   IN  grid            coamps_domain to base conversion on
    !   IN  lat             point's latitude
    !   IN  lon             point's longitude
    !   OUT grid_i          point's i-coordinate
    !   OUT grid_j          point's j-coordinate
    subroutine latlon_to_grid_point(grid, lat, lon, grid_i, grid_j)
        type(coamps_grid), intent(in)  :: grid
        real(kind=r8),     intent(in)  :: lat
        real(kind=r8),     intent(in)  :: lon
        real(kind=r8),     intent(out) :: grid_i
        real(kind=r8),     intent(out) :: grid_j

        ! Work-around for not being able to pass in scalars as arrays
        ! of length 1
        real(kind=r8), dimension(SINGLE_POINT) :: ii_array,  jj_array
        real(kind=r8), dimension(SINGLE_POINT) :: lat_array, lon_array

        lat_array(SINGLE_POINT) = lat
        lon_array(SINGLE_POINT) = lon

        call ll2ij(grid%map_proj, grid%reflat, grid%reflon, grid%iref, &
                   grid%jref, grid%stdlat1, grid%stdlat2, grid%stdlon, &
                   grid%coarse_delta_x, grid%coarse_delta_y, ii_array, &
                   jj_array, SINGLE_POINT, lat_array, lon_array)  

        grid_i = ii_array(SINGLE_POINT)
        grid_j = jj_array(SINGLE_POINT)

        !write (*,'(F9.5 X F9.5 X F9.5 X F9.5 X)') lat, lon, grid_i, grid_j
    end subroutine latlon_to_grid_point

    ! grid_point_to_latlon
    ! --------------------
    ! Converts a lat/lon point to its representation on the COAMPS grid
    !  PARAMETERS
    !   IN  grid            coamps_domain to base conversion on
    !   IN  grid_i          point's i-coordinate
    !   IN  grid_j          point's j-coordinate
    !   OUT lat             point's latitude
    !   OUT lon             point's longitude
    subroutine grid_point_to_latlon(grid, grid_i, grid_j, lat, lon)
        type(coamps_grid), intent(in)  :: grid
        real(kind=r8),     intent(in)  :: grid_i
        real(kind=r8),     intent(in)  :: grid_j
        real(kind=r8),     intent(out) :: lat
        real(kind=r8),     intent(out) :: lon

        ! Work-around for not being able to pass in scalars as arrays
        ! of length 1
        real(kind=r8), dimension(SINGLE_POINT) :: ii_array,  jj_array
        real(kind=r8), dimension(SINGLE_POINT) :: lat_array, lon_array

        ii_array(SINGLE_POINT) = grid_i
        jj_array(SINGLE_POINT) = grid_j

        call ij2ll(grid%map_proj, grid%reflat, grid%reflon, grid%iref, &
                   grid%jref, grid%stdlat1, grid%stdlat2, grid%stdlon, &
                   grid%coarse_delta_x, grid%coarse_delta_y, ii_array, &
                   jj_array, SINGLE_POINT, lat_array, lon_array)  

        lat = lat_array(SINGLE_POINT)
        lon = lon_array(SINGLE_POINT)

        !write (*,'(F9.5 X F9.5 X F9.5 X F9.5 X)') grid_i, grid_j, lat, lon
    end subroutine grid_point_to_latlon

    ! calc_map_constant
    ! ----------------
    ! Given a coamps grid, calculates the map constant for the projection
    !  PARAMETERS
    !   IN  grid              coamps_grid structure
    function calc_map_constant(grid) result(map_constant)
      type(coamps_grid), intent(in)  :: grid
      real(kind=r8)                  :: map_constant

    select case(int(grid%map_proj))
      case(PROJ_MERCATOR)
        map_constant=0.0_r8
      case(PROJ_LAMBERT_CONIC_CONFORMAL)
        if (grid%stdlat1 == grid%stdlat2) then
          map_constant=sin(abs(grid%stdlat1)*DEG_TO_RAD)
        else
          map_constant=( log(sin((90.0_r8-abs(grid%stdlat1))*DEG_TO_RAD))      &
                        -log(sin((90.0_r8-abs(grid%stdlat2))*DEG_TO_RAD)))     &
                       /(log(tan((90.0_r8-abs(grid%stdlat1))*DEG_TO_RAD*0.5_r8))  &
                        -log(tan((90.0_r8-abs(grid%stdlat2))*DEG_TO_RAD*0.5_r8)))
        endif
      case(PROJ_POLAR_STEREOGRAPHIC)
        map_constant=1.0_r8
      case default
        map_constant=1.0_r8
    end select
    end function calc_map_constant

    ! dump_grid_info
    ! --------------
    ! Dumps the COAMPS grid information in a human-readable form
    !  PARAMETERS
    !   IN  grid            COAMPS grid to dump
    subroutine dump_grid_info(grid)
        type(coamps_grid), intent(in) :: grid

        write (*,*) "**** COAMPS MAP ****"
        write (*,*) "Map projection type:      ", grid%map_proj
        write (*,*) "Latitude intersection 1:  ", grid%stdlat1
        write (*,*) "Latitude intersection 2:  ", grid%stdlat2
        write (*,*) "Longitude intersection:   ", grid%stdlon
        write (*,*) "Reference latitude:       ", grid%reflat
        write (*,*) "Reference longitude:      ", grid%reflon
        write (*,*) "i point at ref. lat/lon:  ", grid%iref
        write (*,*) "j point at ref. lat/lon:  ", grid%jref
    end subroutine dump_grid_info


    !------------------------------
    ! END PUBLIC ROUTINES
    !------------------------------

    !------------------------------
    ! BEGIN PRIVATE ROUTINES
    !------------------------------
    !------------------------------
    ! END PRIVATE ROUTINES
    !------------------------------

end module coamps_map_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
