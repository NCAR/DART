
   module mapper_types

   implicit none

   integer,parameter :: grid_max = 100
   integer :: grid_cnt = 0
   integer :: grid_ndx

   TYPE proj_info
      INTEGER          :: code     ! Integer code for projection TYPE
      INTEGER          :: nlat     ! For Gaussian -- number of latitude points 
                                   !  north of the equator 
      INTEGER          :: nlon     !
                                   !
      INTEGER          :: ixdim    ! For Rotated Lat/Lon -- number of mass points
                                   !  in an odd row
      INTEGER          :: jydim    ! For Rotated Lat/Lon -- number of rows
      INTEGER          :: stagger  ! For Rotated Lat/Lon -- mass or velocity grid 
      REAL             :: phi      ! For Rotated Lat/Lon -- domain half-extent in 
                                   !  degrees latitude
      REAL             :: lambda   ! For Rotated Lat/Lon -- domain half-extend in
                                   !  degrees longitude
      REAL             :: lat1     ! SW latitude (1,1) in degrees (-90->90N)
      REAL             :: lon1     ! SW longitude (1,1) in degrees (-180->180E)
      REAL             :: lat0     ! For Cassini, latitude of projection pole
      REAL             :: lon0     ! For Cassini, longitude of projection pole
      REAL             :: dx       ! Grid spacing in meters at truelats, used
                                   !  only for ps, lc, and merc projections
      REAL             :: dy       ! Grid spacing in meters at truelats, used
                                   !  only for ps, lc, and merc projections
      REAL             :: latinc   ! Latitude increment for cylindrical lat/lon
      REAL             :: loninc   ! Longitude increment for cylindrical lat/lon
                                   !  also the lon increment for Gaussian grid
      REAL             :: dlat     ! Lat increment for lat/lon grids
      REAL             :: dlon     ! Lon increment for lat/lon grids
      REAL             :: stdlon   ! Longitude parallel to y-axis (-180->180E)
      REAL             :: truelat1 ! First true latitude (all projections)
      REAL             :: truelat2 ! Second true lat (LC only)
      REAL             :: hemi     ! 1 for NH, -1 for SH
      REAL             :: cone     ! Cone factor for LC projections
      REAL             :: polei    ! Computed i-location of pole point
      REAL             :: polej    ! Computed j-location of pole point
      REAL             :: rsw      ! Computed radius to SW corner
      REAL             :: rebydx   ! Earth radius divided by dx
      REAL             :: knowni   ! X-location of known lat/lon
      REAL             :: knownj   ! Y-location of known lat/lon
      REAL             :: re_m     ! Radius of spherical earth, meters
      REAL             :: rho0     ! For Albers equal area
      REAL             :: nc       ! For Albers equal area
      REAL             :: bigc     ! For Albers equal area
      LOGICAL          :: init     ! Flag to indicate if this struct is 
                                   !  ready for use
      LOGICAL          :: wrap     ! For Gaussian -- flag to indicate wrapping 
                                   !  around globe?
      LOGICAL          :: comp_ll  ! Work in computational lat/lon space for Cassini
      REAL, POINTER, DIMENSION(:) :: gauss_lat  ! Latitude array for Gaussian grid
   END TYPE proj_info

   TYPE area_type
     integer :: lon_s, lon_e
     integer :: lat_s, lat_e
     integer :: active_dcell_cnt
     integer :: total_dcell_cnt
     integer :: interior_dcell_cnt
     integer :: partial_dcell_cnt
     integer, allocatable :: dcell_lon_ndx(:)
     integer, allocatable :: dcell_lat_ndx(:)
     real, allocatable    :: wght(:)
     logical :: has_data
   END TYPE area_type

   TYPE grid_type
     integer :: nlons
     integer :: nlats
     integer, allocatable :: ix(:,:,:)                    ! index used by bilin interpolation
     integer, allocatable :: jy(:,:,:)                    ! index used by bilin interpolation
     real, allocatable :: lon(:)
     real, allocatable :: lat(:)
     real, allocatable :: ax(:,:,:)                        ! weight coef for bilin interpolation
     real, allocatable :: by(:,:,:)                        ! weight coef for bilin interpolation
     real(8), allocatable :: xedge(:)
     real(8), allocatable :: yedge(:)
     logical :: has_area_map
     logical :: has_lon_shift
     logical :: reorder_lons
     logical :: reorder_lats
     type(area_type), allocatable :: model_area_type(:,:)
   END TYPE grid_type

   type(grid_type), target :: grid_specs(grid_max)

   end module mapper_types
