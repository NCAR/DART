! This code is not protected by the DART copyright agreement.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE CONSTANTS_MODULE
!
! This module defines constants that are used by other modules 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module constants_module

   use     types_mod, only : r8, digits12

   real (kind=r8), parameter :: PI = 3.141592653589793_r8
   real (kind=r8), parameter :: OMEGA_E = 0.00007292_r8 ! Angular rotation rate of the earth

   real (kind=r8), parameter :: DEG_PER_RAD = 180.0_r8 / PI
   real (kind=r8), parameter :: RAD_PER_DEG = PI / 180.0_r8

   ! Mean Earth Radius in m.  The value below is consistent
   ! with NCEP's routines and grids.
   real (kind=r8), parameter :: A_WGS84  = 6378137.0_r8
   real (kind=r8), parameter :: B_WGS84  = 6356752.314_r8
   real (kind=r8), parameter :: RE_WGS84 = A_WGS84
   real (kind=r8), parameter :: E_WGS84  = 0.081819192_r8

   real (kind=r8), parameter :: A_NAD83  = 6378137.0_r8
   real (kind=r8), parameter :: RE_NAD83 = A_NAD83
   real (kind=r8), parameter :: E_NAD83  = 0.0818187034_r8

   real (kind=r8), parameter :: EARTH_RADIUS_M = 6370000.0_r8   ! same as MM5 system
   real (kind=r8), parameter :: EARTH_CIRC_M = 2.0_r8 * PI * EARTH_RADIUS_M

end module constants_module


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE MISC_DEFINITIONS_MODULE
!
! This module defines various non-meteorological constants that are used 
!   by other modules for readability.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module misc_definitions_module

   use     types_mod, only : r8

   real (kind=r8), parameter :: NAN=1.E20_r8

   real (kind=r8), parameter :: NOT_MASKED   = -2.0_r8,  &
                      MASKED_BOTH  = -1.0_r8,  &
                      MASKED_WATER =  0.0_r8,  &
                      MASKED_LAND  =  1.0_r8

   integer, parameter :: OUTSIDE_DOMAIN=1E8, NOT_PROCESSED=1E9, INVALID=1E9

   integer, parameter :: SIXTEEN_POINT=1, FOUR_POINT=2, N_NEIGHBOR=3, &
                         AVERAGE4=4, AVERAGE16=5, W_AVERAGE4=6, W_AVERAGE16=7, &
                         SEARCH=8

   integer, parameter :: BOTTOM_TOP=1, TOP_BOTTOM=2

   integer, parameter :: CONTINUOUS=0, CATEGORICAL=1, SP_CONTINUOUS=2

   integer, parameter :: M=1, U=2, V=3, HH=4, VV=5

   integer, parameter :: ONETWOONE=1, SMTHDESMTH=2, SMTHDESMTH_SPECIAL=3

   integer, parameter :: BINARY=1, NETCDF=2, GRIB1=3, HDF=4

   integer, parameter :: BIG_ENDIAN=0, LITTLE_ENDIAN=1

   ! Projection codes for proj_info structure:
   INTEGER, PUBLIC, PARAMETER  :: PROJ_LATLON       = 0
   INTEGER, PUBLIC, PARAMETER  :: PROJ_LC           = 1
   INTEGER, PUBLIC, PARAMETER  :: PROJ_PS           = 2
   INTEGER, PUBLIC, PARAMETER  :: PROJ_PS_WGS84     = 102
   INTEGER, PUBLIC, PARAMETER  :: PROJ_MERC         = 3
   INTEGER, PUBLIC, PARAMETER  :: PROJ_GAUSS        = 4
   INTEGER, PUBLIC, PARAMETER  :: PROJ_CYL          = 5
   INTEGER, PUBLIC, PARAMETER  :: PROJ_CASSINI      = 6
   INTEGER, PUBLIC, PARAMETER  :: PROJ_ALBERS_NAD83 = 105 
   INTEGER, PUBLIC, PARAMETER  :: PROJ_ROTLL        = 203

end module misc_definitions_module


MODULE map_utils

! Module that defines constants, data structures, and
! subroutines used to convert grid indices to lat/lon
! and vice versa.   
!
! SUPPORTED PROJECTIONS
! ---------------------
! Cylindrical Lat/Lon (code = PROJ_LATLON)
! Mercator (code = PROJ_MERC)
! Lambert Conformal (code = PROJ_LC)
! Gaussian (code = PROJ_GAUSS)
! Polar Stereographic (code = PROJ_PS)
! Rotated Lat/Lon (code = PROJ_ROTLL)
!
! REMARKS
! -------
! The routines contained within were adapted from routines
! obtained from NCEP's w3 library.  The original NCEP routines were less
! flexible (e.g., polar-stereo routines only supported truelat of 60N/60S)
! than what we needed, so modifications based on equations in Hoke, Hayes, and
! Renninger (AFGWC/TN/79-003) were added to improve the flexibility.  
! Additionally, coding was improved to F90 standards and the routines were
! combined into this module.  
!
! ASSUMPTIONS
! -----------
!  Grid Definition:
!    For mercator, lambert conformal, and polar-stereographic projections,
!    the routines within assume the following:
!
!       1.  Grid is dimensioned (i,j) where i is the East-West direction,
!           positive toward the east, and j is the north-south direction,
!           positive toward the north.
!       2.  Origin is at (1,1) and is located at the southwest corner,
!           regardless of hemispere.
!       3.  Grid spacing (dx) is always positive.
!       4.  Values of true latitudes must be positive for NH domains
!           and negative for SH domains.
!
!     For the latlon and Gaussian projection, the grid origin may be at any
!     of the corners, and the deltalat and deltalon values can be signed to
!     account for this using the following convention:
!       Origin Location        Deltalat Sign      Deltalon Sign
!       ---------------        -------------      -------------
!        SW Corner                  +                   +
!        NE Corner                  -                   -
!        NW Corner                  -                   +
!        SE Corner                  +                   -
!
!  Data Definitions:
!       1. Any arguments that are a latitude value are expressed in
!          degrees north with a valid range of -90 -> 90
!       2. Any arguments that are a longitude value are expressed in
!          degrees east with a valid range of -180 -> 180.
!       3. Distances are in meters and are always positive.
!       4. The standard longitude (stdlon) is defined as the longitude
!          line which is parallel to the grid's y-axis (j-direction), along
!          which latitude increases (NOT the absolute value of latitude, but
!          the actual latitude, such that latitude increases continuously
!          from the south pole to the north pole) as j increases.
!       5. One true latitude value is required for polar-stereographic and
!          mercator projections, and defines at which latitude the
!          grid spacing is true.  For lambert conformal, two true latitude
!          values must be specified, but may be set equal to each other to
!          specify a tangent projection instead of a secant projection.
!
! USAGE
! -----
! To use the routines in this module, the calling routines must have the
! following statement at the beginning of its declaration block:
!   USE map_utils
!
! The use of the module not only provides access to the necessary routines,
! but also defines a structure of TYPE (proj_info) that can be used
! to declare a variable of the same type to hold your map projection
! information.  It also defines some integer parameters that contain
! the projection codes so one only has to use those variable names rather
! than remembering the acutal code when using them.  The basic steps are
! as follows:
!
!   1.  Ensure the "USE map_utils" is in your declarations.
!   2.  Declare the projection information structure as type(proj_info):
!         TYPE(proj_info) :: proj
!   3.  Populate your structure by calling the map_set routine:
!         CALL map_set(code,lat1,lon1,knowni,knownj,dx,stdlon,truelat1,truelat2,proj)
!       where:
!         code (input) = one of PROJ_LATLON, PROJ_MERC, PROJ_LC, PROJ_PS,
!                        PROJ_GAUSS, or PROJ_ROTLL
!         lat1 (input) = Latitude of grid origin point (i,j)=(1,1)
!                         (see assumptions!)
!         lon1 (input) = Longitude of grid origin
!         knowni (input) = origin point, x-location
!         knownj (input) = origin point, y-location
!         dx (input) = grid spacing in meters (ignored for LATLON projections)
!         stdlon (input) = Standard longitude for PROJ_PS and PROJ_LC,
!               deltalon (see assumptions) for PROJ_LATLON,
!               ignored for PROJ_MERC
!         truelat1 (input) = 1st true latitude for PROJ_PS, PROJ_LC, and
!                PROJ_MERC, deltalat (see assumptions) for PROJ_LATLON
!         truelat2 (input) = 2nd true latitude for PROJ_LC,
!                ignored for all others.
!         proj (output) = The structure of type (proj_info) that will be fully
!                populated after this call
!
!   4.  Now that the proj structure is populated, you may call either
!       of the following routines:
!
!       latlon_to_ij(proj, lat, lon, i, j)
!       ij_to_latlon(proj, i, j, lat, lon)
!
!       It is incumbent upon the calling routine to determine whether or
!       not the values returned are within your domain's bounds.  All values
!       of i, j, lat, and lon are REAL values.
!
!
! REFERENCES
! ----------
!  Hoke, Hayes, and Renninger, "Map Projections and Grid Systems for
!       Meteorological Applications." AFGWC/TN-79/003(Rev), Air Weather
!       Service, 1985.
!
!  NCAR MM5v3 Modeling System, REGRIDDER program, module_first_guess_map.F
!  NCEP routines w3fb06, w3fb07, w3fb08, w3fb09, w3fb11, w3fb12
!
! HISTORY
! -------
! 27 Mar 2001 - Original Version
!               Brent L. Shaw, NOAA/FSL (CSU/CIRA)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! other modules included in this .f90 file
   use constants_module
   use misc_definitions_module

   ! external modules
   use     types_mod, only : r8
!   use utilities_mod, only : register_module

   ! Define some private constants
   INTEGER, PRIVATE, PARAMETER :: HIGH = digits12

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
      REAL(r8)             :: phi      ! For Rotated Lat/Lon -- domain half-extent in 
                                   !  degrees latitude
      REAL(r8)             :: lambda   ! For Rotated Lat/Lon -- domain half-extend in
                                   !  degrees longitude
      REAL(r8)             :: lat1     ! SW latitude (1,1) in degrees (-90->90N)
      REAL(r8)             :: lon1     ! SW longitude (1,1) in degrees (-180->180E)
      REAL(r8)             :: lat0     ! For Cassini, latitude of projection pole
      REAL(r8)             :: lon0     ! For Cassini, longitude of projection pole
      REAL(r8)             :: dx       ! Grid spacing in meters at truelats, used
                                   !  only for ps, lc, and merc projections
      REAL(r8)             :: dy       ! Grid spacing in meters at truelats, used
                                   !  only for ps, lc, and merc projections
      REAL(r8)             :: latinc   ! Latitude increment for cylindrical lat/lon
      REAL(r8)             :: loninc   ! Longitude increment for cylindrical lat/lon
                                   !  also the lon increment for Gaussian grid
      REAL(r8)             :: dlat     ! Lat increment for lat/lon grids
      REAL(r8)             :: dlon     ! Lon increment for lat/lon grids
      REAL(r8)             :: stdlon   ! Longitude parallel to y-axis (-180->180E)
      REAL(r8)             :: truelat1 ! First true latitude (all projections)
      REAL(r8)             :: truelat2 ! Second true lat (LC only)
      REAL(r8)             :: hemi     ! 1 for NH, -1 for SH
      REAL(r8)             :: cone     ! Cone factor for LC projections
      REAL(r8)             :: polei    ! Computed i-location of pole point
      REAL(r8)             :: polej    ! Computed j-location of pole point
      REAL(r8)             :: rsw      ! Computed radius to SW corner
      REAL(r8)             :: rebydx   ! Earth radius divided by dx
      REAL(r8)             :: knowni   ! X-location of known lat/lon
      REAL(r8)             :: knownj   ! Y-location of known lat/lon
      REAL(r8)             :: re_m     ! Radius of spherical earth, meters
      REAL(r8)             :: rho0     ! For Albers equal area
      REAL(r8)             :: nc       ! For Albers equal area
      REAL(r8)             :: bigc     ! For Albers equal area
      LOGICAL          :: init     ! Flag to indicate if this struct is 
                                   !  ready for use
      LOGICAL          :: wrap     ! For Gaussian -- flag to indicate wrapping 
                                   !  around globe?
      REAL(r8), POINTER, DIMENSION(:) :: gauss_lat  ! Latitude array for Gaussian grid

   END TYPE proj_info

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CONTAINS
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE map_init(proj)
      ! Initializes the map projection structure to missing values

      IMPLICIT NONE
      TYPE(proj_info), INTENT(INOUT)  :: proj

      proj%lat1     = -999.9_r8
      proj%lon1     = -999.9_r8
      proj%lat0     = -999.9_r8
      proj%lon0     = -999.9_r8
      proj%dx       = -999.9_r8
      proj%dy       = -999.9_r8
      proj%latinc   = -999.9_r8
      proj%loninc   = -999.9_r8
      proj%stdlon   = -999.9_r8
      proj%truelat1 = -999.9_r8
      proj%truelat2 = -999.9_r8
      proj%phi      = -999.9_r8
      proj%lambda   = -999.9_r8
      proj%ixdim    = -999
      proj%jydim    = -999
      proj%stagger  = HH
      proj%nlat     = 0
      proj%nlon     = 0
      proj%hemi     = 0.0_r8
      proj%cone     = -999.9_r8
      proj%polei    = -999.9_r8
      proj%polej    = -999.9_r8
      proj%rsw      = -999.9_r8
      proj%knowni   = -999.9_r8
      proj%knownj   = -999.9_r8
      proj%re_m     = EARTH_RADIUS_M
      proj%init     = .FALSE.
      proj%wrap     = .FALSE.
      proj%rho0     = 0.0_r8
      proj%nc       = 0.0_r8
      proj%bigc     = 0.0_r8
      nullify(proj%gauss_lat)

   END SUBROUTINE map_init


!nc -- Global WRF assumes proj_code = PROJ_CASSINI (=6 in misc_definitions_module)

   SUBROUTINE map_set(proj_code, proj, lat1, lon1, lat0, lon0, knowni, knownj, dx, latinc, &
                      loninc, stdlon, truelat1, truelat2, nlat, nlon, ixdim, jydim, &
                      stagger, phi, lambda, r_earth)
      ! Given a partially filled proj_info structure, this routine computes
      ! polei, polej, rsw, and cone (if LC projection) to complete the 
      ! structure.  This allows us to eliminate redundant calculations when
      ! calling the coordinate conversion routines multiple times for the
      ! same map.
      ! This will generally be the first routine called when a user wants
      ! to be able to use the coordinate conversion routines, and it
      ! will call the appropriate subroutines based on the 
      ! proj%code which indicates which projection type this is.

      IMPLICIT NONE

      ! Declare arguments
      INTEGER, INTENT(IN)               :: proj_code
      INTEGER, INTENT(IN), OPTIONAL     :: nlat
      INTEGER, INTENT(IN), OPTIONAL     :: nlon
      INTEGER, INTENT(IN), OPTIONAL     :: ixdim
      INTEGER, INTENT(IN), OPTIONAL     :: jydim
      INTEGER, INTENT(IN), OPTIONAL     :: stagger
      REAL(r8), INTENT(IN), OPTIONAL        :: latinc
      REAL(r8), INTENT(IN), OPTIONAL        :: loninc
      REAL(r8), INTENT(IN), OPTIONAL        :: lat1
      REAL(r8), INTENT(IN), OPTIONAL        :: lon1
      REAL(r8), INTENT(IN), OPTIONAL        :: lat0
      REAL(r8), INTENT(IN), OPTIONAL        :: lon0
      REAL(r8), INTENT(IN), OPTIONAL        :: dx
      REAL(r8), INTENT(IN), OPTIONAL        :: stdlon
      REAL(r8), INTENT(IN), OPTIONAL        :: truelat1
      REAL(r8), INTENT(IN), OPTIONAL        :: truelat2
      REAL(r8), INTENT(IN), OPTIONAL        :: knowni
      REAL(r8), INTENT(IN), OPTIONAL        :: knownj
      REAL(r8), INTENT(IN), OPTIONAL        :: phi
      REAL(r8), INTENT(IN), OPTIONAL        :: lambda
      REAL(r8), INTENT(IN), OPTIONAL        :: r_earth
      TYPE(proj_info), INTENT(OUT)      :: proj

      INTEGER :: iter
      REAL(r8) :: dummy_lon1
      REAL(r8) :: dummy_lon0
      REAL(r8) :: dummy_stdlon

      ! First, verify that mandatory parameters are present for the specified proj_code
      IF ( proj_code == PROJ_LC ) THEN
         IF ( .NOT.PRESENT(truelat1) .OR. &
              .NOT.PRESENT(truelat2) .OR. &
              .NOT.PRESENT(lat1) .OR. &
              .NOT.PRESENT(lon1) .OR. &
              .NOT.PRESENT(knowni) .OR. &
              .NOT.PRESENT(knownj) .OR. &
              .NOT.PRESENT(stdlon) .OR. &
              .NOT.PRESENT(dx) ) THEN
            PRINT '(A,I2)', 'The following are mandatory parameters for projection code : ', proj_code
            PRINT '(A)', ' truelat1, truelat2, lat1, lon1, knowni, knownj, stdlon, dx'
            write(6,*) 'ERROR: MAP_INIT'
         END IF
      ELSE IF ( proj_code == PROJ_PS ) THEN
         IF ( .NOT.PRESENT(truelat1) .OR. &
              .NOT.PRESENT(lat1) .OR. &
              .NOT.PRESENT(lon1) .OR. &
              .NOT.PRESENT(knowni) .OR. &
              .NOT.PRESENT(knownj) .OR. &
              .NOT.PRESENT(stdlon) .OR. &
              .NOT.PRESENT(dx) ) THEN
            PRINT '(A,I2)', 'The following are mandatory parameters for projection code : ', proj_code
            PRINT '(A)', ' truelat1, lat1, lon1, knonwi, knownj, stdlon, dx'
            write(6,*) 'ERROR: MAP_INIT'
         END IF
      ELSE IF ( proj_code == PROJ_PS_WGS84 ) THEN
         IF ( .NOT.PRESENT(truelat1) .OR. &
              .NOT.PRESENT(lat1) .OR. &
              .NOT.PRESENT(lon1) .OR. &
              .NOT.PRESENT(knowni) .OR. &
              .NOT.PRESENT(knownj) .OR. &
              .NOT.PRESENT(stdlon) .OR. &
              .NOT.PRESENT(dx) ) THEN
            PRINT '(A,I2)', 'The following are mandatory parameters for projection code : ', proj_code
            PRINT '(A)', ' truelat1, lat1, lon1, knonwi, knownj, stdlon, dx'
            write(6,*) 'ERROR: MAP_INIT'
         END IF
      ELSE IF ( proj_code == PROJ_ALBERS_NAD83 ) THEN
         IF ( .NOT.PRESENT(truelat1) .OR. &
              .NOT.PRESENT(truelat2) .OR. &
              .NOT.PRESENT(lat1) .OR. &
              .NOT.PRESENT(lon1) .OR. &
              .NOT.PRESENT(knowni) .OR. &
              .NOT.PRESENT(knownj) .OR. &
              .NOT.PRESENT(stdlon) .OR. &
              .NOT.PRESENT(dx) ) THEN
            PRINT '(A,I2)', 'The following are mandatory parameters for projection code : ', proj_code
            PRINT '(A)', ' truelat1, truelat2, lat1, lon1, knonwi, knownj, stdlon, dx'
            write(6,*) 'ERROR: MAP_INIT'
         END IF
      ELSE IF ( proj_code == PROJ_MERC ) THEN
         IF ( .NOT.PRESENT(truelat1) .OR. &
              .NOT.PRESENT(lat1) .OR. &
              .NOT.PRESENT(lon1) .OR. &
              .NOT.PRESENT(knowni) .OR. &
              .NOT.PRESENT(knownj) .OR. &
              .NOT.PRESENT(dx) ) THEN
            PRINT '(A,I2)', 'The following are mandatory parameters for projection code : ', proj_code
            PRINT '(A)', ' truelat1, lat1, lon1, knowni, knownj, dx'
            write(6,*) 'ERROR: MAP_INIT'
         END IF
      ELSE IF ( proj_code == PROJ_LATLON ) THEN
         IF ( .NOT.PRESENT(latinc) .OR. &
              .NOT.PRESENT(loninc) .OR. &
              .NOT.PRESENT(knowni) .OR. &
              .NOT.PRESENT(knownj) .OR. &
              .NOT.PRESENT(lat1) .OR. &
              .NOT.PRESENT(lon1) ) THEN
            PRINT '(A,I2)', 'The following are mandatory parameters for projection code : ', proj_code
            PRINT '(A)', ' latinc, loninc, knowni, knownj, lat1, lon1'
            write(6,*) 'ERROR: MAP_INIT'
         END IF
      ELSE IF ( proj_code == PROJ_CYL ) THEN
         IF ( .NOT.PRESENT(latinc) .OR. &
              .NOT.PRESENT(loninc) .OR. &
              .NOT.PRESENT(stdlon) ) THEN
            PRINT '(A,I2)', 'The following are mandatory parameters for projection code : ', proj_code
            PRINT '(A)', ' latinc, loninc, stdlon'
            write(6,*) 'ERROR: MAP_INIT'
         END IF
      ELSE IF ( proj_code == PROJ_CASSINI ) THEN
         IF ( .NOT.PRESENT(latinc) .OR. &
              .NOT.PRESENT(loninc) .OR. &
              .NOT.PRESENT(lat1) .OR. &
              .NOT.PRESENT(lon1) .OR. &
              .NOT.PRESENT(lat0) .OR. &
              .NOT.PRESENT(lon0) .OR. &
              .NOT.PRESENT(knowni) .OR. &
              .NOT.PRESENT(knownj) .OR. &
              .NOT.PRESENT(stdlon) ) THEN
            PRINT '(A,I2)', 'The following are mandatory parameters for projection code : ', proj_code
            PRINT '(A)', ' latinc, loninc, lat1, lon1, knowni, knownj, lat0, lon0, stdlon'
            write(6,*) 'ERROR: MAP_INIT'
         END IF
      ELSE IF ( proj_code == PROJ_GAUSS ) THEN
         IF ( .NOT.PRESENT(nlat) .OR. &
              .NOT.PRESENT(lat1) .OR. &
              .NOT.PRESENT(lon1) .OR. &
              .NOT.PRESENT(loninc) ) THEN
            PRINT '(A,I2)', 'The following are mandatory parameters for projection code : ', proj_code
            PRINT '(A)', ' nlat, lat1, lon1, loninc'
            write(6,*) 'ERROR: MAP_INIT'
         END IF
      ELSE IF ( proj_code == PROJ_ROTLL ) THEN
         IF ( .NOT.PRESENT(ixdim) .OR. &
              .NOT.PRESENT(jydim) .OR. &
              .NOT.PRESENT(phi) .OR. &
              .NOT.PRESENT(lambda) .OR. &
              .NOT.PRESENT(lat1) .OR. &
              .NOT.PRESENT(lon1) .OR. &
              .NOT.PRESENT(stagger) ) THEN
            PRINT '(A,I2)', 'The following are mandatory parameters for projection code : ', proj_code
            PRINT '(A)', ' ixdim, jydim, phi, lambda, lat1, lon1, stagger'
            write(6,*) 'ERROR: MAP_INIT'
         END IF
      ELSE
         PRINT '(A,I2)', 'Unknown projection code: ', proj_code
         write(6,*) 'ERROR: MAP_INIT'
      END IF

      ! Check for validity of mandatory variables in proj
      IF ( PRESENT(lat1) ) THEN
         IF ( ABS(lat1) .GT. 90.0_r8 ) THEN
            PRINT '(A)', 'Latitude of origin corner required as follows:'
            PRINT '(A)', '    -90N <= lat1 < = 90.N'
            write(6,*) 'ERROR: MAP_INIT'
         ENDIF
      ENDIF

      IF ( PRESENT(lon1) ) THEN
         dummy_lon1 = lon1
         IF ( ABS(dummy_lon1) .GT. 180.0_r8 ) THEN
            iter = 0 
            DO WHILE (ABS(dummy_lon1) > 180.0_r8 .AND. iter < 10)
               IF (dummy_lon1 < -180.0_r8) dummy_lon1 = dummy_lon1 + 360.0_r8
               IF (dummy_lon1 > 180.0_r8) dummy_lon1 = dummy_lon1 - 360.0_r8
               iter = iter + 1
            END DO
            IF (abs(dummy_lon1) > 180.0_r8) THEN
               PRINT '(A)', 'Longitude of origin required as follows:'
               PRINT '(A)', '   -180E <= lon1 <= 180W'
               write(6,*) 'ERROR: MAP_INIT'
            ENDIF
         ENDIF
      ENDIF

      IF ( PRESENT(lon0) ) THEN
         dummy_lon0 = lon0
         IF ( ABS(dummy_lon0) .GT. 180.0_r8) THEN
            iter = 0 
            DO WHILE (ABS(dummy_lon0) > 180.0_r8 .AND. iter < 10)
               IF (dummy_lon0 < -180.0_r8) dummy_lon0 = dummy_lon0 + 360.0_r8
               IF (dummy_lon0 > 180.0_r8) dummy_lon0 = dummy_lon0 - 360.0_r8
               iter = iter + 1
            END DO
            IF (abs(dummy_lon0) > 180.0_r8) THEN
               PRINT '(A)', 'Longitude of pole required as follows:'
               PRINT '(A)', '   -180E <= lon0 <= 180W'
               write(6,*) 'ERROR: MAP_INIT'
            ENDIF
         ENDIF
      ENDIF

      IF ( PRESENT(dx) ) THEN
         IF ((dx .LE. 0.).AND.(proj_code .NE. PROJ_LATLON)) THEN
            PRINT '(A)', 'Require grid spacing (dx) in meters be positive!'
            write(6,*) 'ERROR: MAP_INIT'
         ENDIF
      ENDIF

      IF ( PRESENT(stdlon) ) THEN
         dummy_stdlon = stdlon
         IF ((ABS(dummy_stdlon) > 180.0_r8).AND.(proj_code /= PROJ_MERC)) THEN
            iter = 0 
            DO WHILE (ABS(dummy_stdlon) > 180.0_r8 .AND. iter < 10)
               IF (dummy_stdlon < -180.0_r8) dummy_stdlon = dummy_stdlon + 360.0_r8
               IF (dummy_stdlon > 180.0_r8) dummy_stdlon = dummy_stdlon - 360.0_r8
               iter = iter + 1
            END DO
            IF (abs(dummy_stdlon) > 180.0_r8) THEN
               PRINT '(A)', 'Need orientation longitude (stdlon) as: '
               PRINT '(A)', '   -180E <= stdlon <= 180W' 
               write(6,*) 'ERROR: MAP_INIT'
            ENDIF
         ENDIF
      ENDIF

      IF ( PRESENT(truelat1) ) THEN
         IF (ABS(truelat1).GT.90.0_r8) THEN
            PRINT '(A)', 'Set true latitude 1 for all projections!'
            write(6,*) 'ERROR: MAP_INIT'
         ENDIF
      ENDIF

      CALL map_init(proj) 
      proj%code  = proj_code
      IF ( PRESENT(lat1) )     proj%lat1     = lat1
      IF ( PRESENT(lon1) )     proj%lon1     = dummy_lon1
      IF ( PRESENT(lat0) )     proj%lat0     = lat0
      IF ( PRESENT(lon0) )     proj%lon0     = dummy_lon0
      IF ( PRESENT(latinc) )   proj%latinc   = latinc
      IF ( PRESENT(loninc) )   proj%loninc   = loninc
      IF ( PRESENT(knowni) )   proj%knowni   = knowni
      IF ( PRESENT(knownj) )   proj%knownj   = knownj
      IF ( PRESENT(dx) )       proj%dx       = dx
      IF ( PRESENT(stdlon) )   proj%stdlon   = dummy_stdlon
      IF ( PRESENT(truelat1) ) proj%truelat1 = truelat1
      IF ( PRESENT(truelat2) ) proj%truelat2 = truelat2
      IF ( PRESENT(nlat) )     proj%nlat     = nlat
      IF ( PRESENT(nlon) )     proj%nlon     = nlon
      IF ( PRESENT(ixdim) )    proj%ixdim    = ixdim
      IF ( PRESENT(jydim) )    proj%jydim    = jydim
      IF ( PRESENT(stagger) )  proj%stagger  = stagger
      IF ( PRESENT(phi) )      proj%phi      = phi
      IF ( PRESENT(lambda) )   proj%lambda   = lambda
      IF ( PRESENT(r_earth) )  proj%re_m     = r_earth

      IF ( PRESENT(dx) ) THEN 
         IF ( (proj_code == PROJ_LC) .OR. (proj_code == PROJ_PS) .OR. &
              (proj_code == PROJ_PS_WGS84) .OR. (proj_code == PROJ_ALBERS_NAD83) .OR. &
              (proj_code == PROJ_MERC) ) THEN
            proj%dx = dx
            IF (truelat1 .LT. 0.0_r8) THEN
               proj%hemi = -1.0_r8 
            ELSE
               proj%hemi = 1.0_r8
            ENDIF
            proj%rebydx = proj%re_m / dx
         ENDIF
      ENDIF

      pick_proj: SELECT CASE(proj%code)

         CASE(PROJ_PS)
            CALL set_ps(proj)

         CASE(PROJ_PS_WGS84)
            CALL set_ps_wgs84(proj)

         CASE(PROJ_ALBERS_NAD83)
            CALL set_albers_nad83(proj)

         CASE(PROJ_LC)
            IF (ABS(proj%truelat2) .GT. 90.0_r8) THEN
               proj%truelat2=proj%truelat1
            ENDIF
            CALL set_lc(proj)

         CASE (PROJ_MERC)
            CALL set_merc(proj)

         CASE (PROJ_LATLON)

         CASE (PROJ_GAUSS)
            CALL set_gauss(proj)

         CASE (PROJ_CYL)
            CALL set_cyl(proj)

         CASE (PROJ_CASSINI)
            CALL set_cassini(proj)

         CASE (PROJ_ROTLL)

      END SELECT pick_proj
      proj%init = .TRUE.

      RETURN

   END SUBROUTINE map_set


   SUBROUTINE latlon_to_ij(proj, lat, lon, i, j)
      ! Converts input lat/lon values to the cartesian (i,j) value
      ! for the given projection. 

      IMPLICIT NONE
      TYPE(proj_info), INTENT(IN)          :: proj
      REAL(r8), INTENT(IN)                     :: lat
      REAL(r8), INTENT(IN)                     :: lon
      REAL(r8), INTENT(OUT)                    :: i
      REAL(r8), INTENT(OUT)                    :: j

      IF (.NOT.proj%init) THEN
         PRINT '(A)', 'You have not called map_set for this projection!'
         write(6,*) 'ERROR: LATLON_TO_IJ'
      ENDIF

      SELECT CASE(proj%code)

         CASE(PROJ_LATLON)
            CALL llij_latlon(lat,lon,proj,i,j)

         CASE(PROJ_MERC)
            CALL llij_merc(lat,lon,proj,i,j)

         CASE(PROJ_PS)
            CALL llij_ps(lat,lon,proj,i,j)

         CASE(PROJ_PS_WGS84)
            CALL llij_ps_wgs84(lat,lon,proj,i,j)

         CASE(PROJ_ALBERS_NAD83)
            CALL llij_albers_nad83(lat,lon,proj,i,j)

         CASE(PROJ_LC)
            CALL llij_lc(lat,lon,proj,i,j)

         CASE(PROJ_GAUSS)
            CALL llij_gauss(lat,lon,proj,i,j)

         CASE(PROJ_CYL)
            CALL llij_cyl(lat,lon,proj,i,j)

         CASE(PROJ_CASSINI)
            CALL llij_cassini(lat,lon,proj,i,j)

         CASE(PROJ_ROTLL)
            CALL llij_rotlatlon(lat,lon,proj,i,j)

         CASE DEFAULT
            PRINT '(A,I2)', 'Unrecognized map projection code: ', proj%code
            write(6,*) 'ERROR: LATLON_TO_IJ'

      END SELECT

      RETURN

   END SUBROUTINE latlon_to_ij


   SUBROUTINE ij_to_latlon(proj, i, j, lat, lon)
      ! Computes geographical latitude and longitude for a given (i,j) point
      ! in a grid with a projection of proj

      IMPLICIT NONE
      TYPE(proj_info),INTENT(IN)          :: proj
      REAL(r8), INTENT(IN)                    :: i
      REAL(r8), INTENT(IN)                    :: j
      REAL(r8), INTENT(OUT)                   :: lat
      REAL(r8), INTENT(OUT)                   :: lon

      IF (.NOT.proj%init) THEN
         PRINT '(A)', 'You have not called map_set for this projection!'
         write(6,*) 'ERROR: IJ_TO_LATLON'
      ENDIF
      SELECT CASE (proj%code)

         CASE (PROJ_LATLON)
            CALL ijll_latlon(i, j, proj, lat, lon)

         CASE (PROJ_MERC)
            CALL ijll_merc(i, j, proj, lat, lon)

         CASE (PROJ_PS)
            CALL ijll_ps(i, j, proj, lat, lon)

         CASE (PROJ_PS_WGS84)
            CALL ijll_ps_wgs84(i, j, proj, lat, lon)

         CASE (PROJ_ALBERS_NAD83)
            CALL ijll_albers_nad83(i, j, proj, lat, lon)

         CASE (PROJ_LC)
            CALL ijll_lc(i, j, proj, lat, lon)

         CASE (PROJ_CYL)
            CALL ijll_cyl(i, j, proj, lat, lon)

         CASE (PROJ_CASSINI)
            CALL ijll_cassini(i, j, proj, lat, lon)

         CASE (PROJ_ROTLL)
            CALL ijll_rotlatlon(i, j, proj, lat, lon)

         CASE DEFAULT
            PRINT '(A,I2)', 'Unrecognized map projection code: ', proj%code
            write(6,*) 'ERROR: IJ_TO_LATLON'

      END SELECT
      RETURN
   END SUBROUTINE ij_to_latlon


   SUBROUTINE set_ps(proj)
      ! Initializes a polar-stereographic map projection from the partially
      ! filled proj structure. This routine computes the radius to the
      ! southwest corner and computes the i/j location of the pole for use
      ! in llij_ps and ijll_ps.
      IMPLICIT NONE

      ! Declare args
      TYPE(proj_info), INTENT(INOUT)    :: proj

      ! Local vars
      REAL(r8)                              :: ala1
      REAL(r8)                              :: alo1
      REAL(r8)                              :: reflon
      REAL(r8)                              :: scale_top

      ! Executable code

      ! Thanks to Kevin Manning for the 'cone' fix.  It must be set to 1.0
      !  to fix wind rotation for polar stereographic projection
      reflon = proj%stdlon + 90.0_r8
      proj%cone = 1.0_r8

      ! Compute numerator term of map scale factor
      scale_top = 1.0_r8 + proj%hemi * SIN(proj%truelat1 * rad_per_deg)

      ! Compute radius to lower-left (SW) corner
      ala1 = proj%lat1 * rad_per_deg
      proj%rsw = proj%rebydx*COS(ala1)*scale_top/(1.0_r8+proj%hemi*SIN(ala1))

      ! Find the pole point
      alo1 = (proj%lon1 - reflon) * rad_per_deg
      proj%polei = proj%knowni - proj%rsw * COS(alo1)
      proj%polej = proj%knownj - proj%hemi * proj%rsw * SIN(alo1)

      RETURN

   END SUBROUTINE set_ps


   SUBROUTINE llij_ps(lat,lon,proj,i,j)
      ! Given latitude (-90 to 90), longitude (-180 to 180), and the
      ! standard polar-stereographic projection information via the 
      ! public proj structure, this routine returns the i/j indices which
      ! if within the domain range from 1->nx and 1->ny, respectively.

      IMPLICIT NONE

      ! Delcare input arguments
      REAL(r8), INTENT(IN)               :: lat
      REAL(r8), INTENT(IN)               :: lon
      TYPE(proj_info),INTENT(IN)     :: proj

      ! Declare output arguments     
      REAL(r8), INTENT(OUT)              :: i !(x-index)
      REAL(r8), INTENT(OUT)              :: j !(y-index)

      ! Declare local variables

      REAL(r8)                           :: reflon
      REAL(r8)                           :: scale_top
      REAL(r8)                           :: ala
      REAL(r8)                           :: alo
      REAL(r8)                           :: rm

      ! BEGIN CODE

      reflon = proj%stdlon + 90.0_r8

      ! Compute numerator term of map scale factor

      scale_top = 1.0_r8 + proj%hemi * SIN(proj%truelat1 * rad_per_deg)

      ! Find radius to desired point
      ala = lat * rad_per_deg
      rm = proj%rebydx * COS(ala) * scale_top/(1.0_r8 + proj%hemi *SIN(ala))
      alo = (lon - reflon) * rad_per_deg
      i = proj%polei + rm * COS(alo)
      j = proj%polej + proj%hemi * rm * SIN(alo)

      RETURN

   END SUBROUTINE llij_ps


   SUBROUTINE ijll_ps(i, j, proj, lat, lon)

      ! This is the inverse subroutine of llij_ps.  It returns the 
      ! latitude and longitude of an i/j point given the projection info 
      ! structure.  

      IMPLICIT NONE

      ! Declare input arguments
      REAL(r8), INTENT(IN)                    :: i    ! Column
      REAL(r8), INTENT(IN)                    :: j    ! Row
      TYPE (proj_info), INTENT(IN)        :: proj

      ! Declare output arguments
      REAL(r8), INTENT(OUT)                   :: lat     ! -90 -> 90 north
      REAL(r8), INTENT(OUT)                   :: lon     ! -180 -> 180 East

      ! Local variables
      REAL(r8)                                :: reflon
      REAL(r8)                                :: scale_top
      REAL(r8)                                :: xx,yy
      REAL(r8)                                :: gi2, r2
      REAL(r8)                                :: arccos

      ! Begin Code

      ! Compute the reference longitude by rotating 90 degrees to the east
      ! to find the longitude line parallel to the positive x-axis.
      reflon = proj%stdlon + 90.0_r8

      ! Compute numerator term of map scale factor
      scale_top = 1.0_r8 + proj%hemi * SIN(proj%truelat1 * rad_per_deg)

      ! Compute radius to point of interest
      xx = i - proj%polei
      yy = (j - proj%polej) * proj%hemi
      r2 = xx**2 + yy**2

      ! Now the magic code
      IF (r2 .EQ. 0.0_r8) THEN 
         lat = proj%hemi * 90.0_r8
         lon = reflon
      ELSE
         gi2 = (proj%rebydx * scale_top)**2.
         lat = deg_per_rad * proj%hemi * ASIN((gi2-r2)/(gi2+r2))
         arccos = ACOS(xx/SQRT(r2))
         IF (yy .GT. 0.0_r8) THEN
            lon = reflon + deg_per_rad * arccos
         ELSE
            lon = reflon - deg_per_rad * arccos
         ENDIF
      ENDIF

      ! Convert to a -180 -> 180 East convention
      IF (lon .GT. 180.0_r8) lon = lon - 360.0_r8
      IF (lon .LT. -180.0_r8) lon = lon + 360.0_r8

      RETURN

   END SUBROUTINE ijll_ps


   SUBROUTINE set_ps_wgs84(proj)
      ! Initializes a polar-stereographic map projection (WGS84 ellipsoid) 
      ! from the partially filled proj structure. This routine computes the 
      ! radius to the southwest corner and computes the i/j location of the 
      ! pole for use in llij_ps and ijll_ps.

      IMPLICIT NONE

      ! Arguments
      TYPE(proj_info), INTENT(INOUT)    :: proj

      ! Local variables
      real(r8) :: h, mc, tc, t, rho

      h = proj%hemi

      mc = cos(h*proj%truelat1*rad_per_deg)/sqrt(1.0_r8-(E_WGS84*sin(h*proj%truelat1*rad_per_deg))**2.0)
      tc = sqrt(((1.0_r8-sin(h*proj%truelat1*rad_per_deg))/(1.0_r8+sin(h*proj%truelat1*rad_per_deg)))* &
                (((1.0_r8+E_WGS84*sin(h*proj%truelat1*rad_per_deg))/(1.0_r8 - &
                   E_WGS84*sin(h*proj%truelat1*rad_per_deg)))**E_WGS84 ))

      ! Find the i/j location of reference lat/lon with respect to the pole of the projection
      t = sqrt(((1.0_r8-sin(h*proj%lat1*rad_per_deg))/(1.0_r8+sin(h*proj%lat1*rad_per_deg)))* &
               (((1.0_r8+E_WGS84*sin(h*proj%lat1*rad_per_deg))/(1.0_r8 - &
                  E_WGS84*sin(h*proj%lat1*rad_per_deg)) )**E_WGS84 ) )
      rho = h * (A_WGS84 / proj%dx) * mc * t / tc
      proj%polei = rho * sin((h*proj%lon1 - h*proj%stdlon)*rad_per_deg)
      proj%polej = -rho * cos((h*proj%lon1 - h*proj%stdlon)*rad_per_deg)

      RETURN

   END SUBROUTINE set_ps_wgs84


   SUBROUTINE llij_ps_wgs84(lat,lon,proj,i,j)
      ! Given latitude (-90 to 90), longitude (-180 to 180), and the
      ! standard polar-stereographic projection information via the 
      ! public proj structure, this routine returns the i/j indices which
      ! if within the domain range from 1->nx and 1->ny, respectively.

      IMPLICIT NONE

      ! Arguments
      REAL(r8), INTENT(IN)               :: lat
      REAL(r8), INTENT(IN)               :: lon
      REAL(r8), INTENT(OUT)              :: i !(x-index)
      REAL(r8), INTENT(OUT)              :: j !(y-index)
      TYPE(proj_info),INTENT(IN)     :: proj

      ! Local variables
      real(r8) :: h, mc, tc, t, rho

      h = proj%hemi

      mc = cos(h*proj%truelat1*rad_per_deg)/sqrt(1.0_r8-(E_WGS84*sin(h*proj%truelat1*rad_per_deg))**2.0)
      tc = sqrt(((1.0_r8-sin(h*proj%truelat1*rad_per_deg))/(1.0_r8+sin(h*proj%truelat1*rad_per_deg)))* &
                (((1.0_r8+E_WGS84*sin(h*proj%truelat1*rad_per_deg))/(1.0_r8 - &
                   E_WGS84*sin(h*proj%truelat1*rad_per_deg)))**E_WGS84 ))

      t = sqrt(((1.0_r8-sin(h*lat*rad_per_deg))/(1.0_r8+sin(h*lat*rad_per_deg))) * &
               (((1.0_r8+E_WGS84*sin(h*lat*rad_per_deg))/(1.0_r8 - &
                  E_WGS84*sin(h*lat*rad_per_deg)))**E_WGS84))

      ! Find the x/y location of the requested lat/lon with respect to the pole of the projection
      rho = (A_WGS84 / proj%dx) * mc * t / tc
      i = h *  rho * sin((h*lon - h*proj%stdlon)*rad_per_deg)
      j = h *(-rho)* cos((h*lon - h*proj%stdlon)*rad_per_deg)

      ! Get i/j relative to reference i/j
      i = proj%knowni + (i - proj%polei)
      j = proj%knownj + (j - proj%polej)

      RETURN

   END SUBROUTINE llij_ps_wgs84


   SUBROUTINE ijll_ps_wgs84(i, j, proj, lat, lon)

      ! This is the inverse subroutine of llij_ps.  It returns the 
      ! latitude and longitude of an i/j point given the projection info 
      ! structure.  

      implicit none

      ! Arguments
      REAL(r8), INTENT(IN)                    :: i    ! Column
      REAL(r8), INTENT(IN)                    :: j    ! Row
      REAL(r8), INTENT(OUT)                   :: lat     ! -90 -> 90 north
      REAL(r8), INTENT(OUT)                   :: lon     ! -180 -> 180 East
      TYPE (proj_info), INTENT(IN)        :: proj

      ! Local variables
      real(r8) :: h, mc, tc, t, rho, x, y
      real(r8) :: chi, a, b, c, d

      h = proj%hemi
      x = (i - proj%knowni + proj%polei)
      y = (j - proj%knownj + proj%polej)

      mc = cos(h*proj%truelat1*rad_per_deg)/sqrt(1.0_r8-(E_WGS84*sin(h*proj%truelat1*rad_per_deg))**2.0)
      tc = sqrt(((1.0_r8-sin(h*proj%truelat1*rad_per_deg))/(1.0_r8+sin(h*proj%truelat1*rad_per_deg))) * &
                (((1.0_r8+E_WGS84*sin(h*proj%truelat1*rad_per_deg))/(1.0_r8 - &
                   E_WGS84*sin(h*proj%truelat1*rad_per_deg)))**E_WGS84 ))

      rho = sqrt((x*proj%dx)**2.0 + (y*proj%dx)**2.0)
      t = rho * tc / (A_WGS84 * mc) 

      lon = h*proj%stdlon + h*atan2(h*x,h*(-y))

      chi = PI/2.0_r8-2.0_r8*atan(t)
      a = (1.0_r8/2.0_r8)*E_WGS84**2. + (5.0_r8/24.0_r8)*E_WGS84**4. + (1.0_r8/40.0_r8)*E_WGS84**6. + &
              (73.0_r8/2016.0_r8)*E_WGS84**8.
      b = (7.0_r8/24.0_r8)*E_WGS84**4. + (29.0_r8/120.0_r8)*E_WGS84**6. + &
              (54113.0_r8/40320.0_r8)*E_WGS84**8.
      c = (7.0_r8/30.0_r8)*E_WGS84**6. + (81.0_r8/280.0_r8)*E_WGS84**8.
      d = (4279.0_r8/20160.0_r8)*E_WGS84**8.

      lat = chi + sin(2.0_r8*chi)*(a + cos(2.0_r8*chi)*(b + cos(2.0_r8*chi)*(c + d*cos(2.0_r8*chi))))
      lat = h * lat

      lat = lat*deg_per_rad
      lon = lon*deg_per_rad

      RETURN

   END SUBROUTINE ijll_ps_wgs84


   SUBROUTINE set_albers_nad83(proj)
      ! Initializes an Albers equal area map projection (NAD83 ellipsoid) 
      ! from the partially filled proj structure. This routine computes the 
      ! radius to the southwest corner and computes the i/j location of the 
      ! pole for use in llij_albers_nad83 and ijll_albers_nad83.

      IMPLICIT NONE

      ! Arguments
      TYPE(proj_info), INTENT(INOUT)    :: proj

      ! Local variables
      real(r8) :: h, m1, m2, q1, q2, theta, q, sinphi

      h = proj%hemi

      m1 = cos(h*proj%truelat1*rad_per_deg)/sqrt(1.0_r8-(E_NAD83*sin(h*proj%truelat1*rad_per_deg))**2.0)
      m2 = cos(h*proj%truelat2*rad_per_deg)/sqrt(1.0_r8-(E_NAD83*sin(h*proj%truelat2*rad_per_deg))**2.0)

      sinphi = sin(proj%truelat1*rad_per_deg)
      q1 = (1.0_r8-E_NAD83**2.0) * &
           ((sinphi/(1.0_r8-(E_NAD83*sinphi)**2.0)) - 1.0_r8/(2.0_r8*E_NAD83) * &
             log((1.0_r8-E_NAD83*sinphi)/(1.0_r8+E_NAD83*sinphi)))

      sinphi = sin(proj%truelat2*rad_per_deg)
      q2 = (1.0_r8-E_NAD83**2.0) * &
           ((sinphi/(1.0_r8-(E_NAD83*sinphi)**2.0)) - 1.0_r8/(2.0_r8*E_NAD83) * &
             log((1.0_r8-E_NAD83*sinphi)/(1.0_r8+E_NAD83*sinphi)))

      if (proj%truelat1 == proj%truelat2) then
         proj%nc = sin(proj%truelat1*rad_per_deg)
      else
         proj%nc = (m1**2.0 - m2**2.0) / (q2 - q1)
      end if

      proj%bigc = m1**2.0 + proj%nc*q1

      ! Find the i/j location of reference lat/lon with respect to the pole of the projection
      sinphi = sin(proj%lat1*rad_per_deg)
      q = (1.0_r8-E_NAD83**2.0) * &
           ((sinphi/(1.0_r8-(E_NAD83*sinphi)**2.0)) - 1.0_r8/(2.0_r8*E_NAD83) * &
             log((1.0_r8-E_NAD83*sinphi)/(1.0_r8+E_NAD83*sinphi)))

      proj%rho0 = h * (A_NAD83 / proj%dx) * sqrt(proj%bigc - proj%nc * q) / proj%nc 
      theta = proj%nc*(proj%lon1 - proj%stdlon)*rad_per_deg

      proj%polei = proj%rho0 * sin(h*theta) 
      proj%polej = proj%rho0 - proj%rho0 * cos(h*theta)

      RETURN

   END SUBROUTINE set_albers_nad83


   SUBROUTINE llij_albers_nad83(lat,lon,proj,i,j)
      ! Given latitude (-90 to 90), longitude (-180 to 180), and the
      ! standard projection information via the 
      ! public proj structure, this routine returns the i/j indices which
      ! if within the domain range from 1->nx and 1->ny, respectively.

      IMPLICIT NONE

      ! Arguments
      REAL(r8), INTENT(IN)               :: lat
      REAL(r8), INTENT(IN)               :: lon
      REAL(r8), INTENT(OUT)              :: i !(x-index)
      REAL(r8), INTENT(OUT)              :: j !(y-index)
      TYPE(proj_info),INTENT(IN)     :: proj

      ! Local variables
      real(r8) :: h, q, rho, theta, sinphi

      h = proj%hemi

      sinphi = sin(h*lat*rad_per_deg)

      ! Find the x/y location of the requested lat/lon with respect to the pole of the projection
      q = (1.0_r8-E_NAD83**2.0) * &
           ((sinphi/(1.0_r8-(E_NAD83*sinphi)**2.0)) - 1.0_r8/(2.0_r8*E_NAD83) * &
             log((1.0_r8-E_NAD83*sinphi)/(1.0_r8+E_NAD83*sinphi)))

      rho = h * (A_NAD83 / proj%dx) * sqrt(proj%bigc - proj%nc * q) / proj%nc
      theta = proj%nc * (h*lon - h*proj%stdlon)*rad_per_deg

      i = h*rho*sin(theta)
      j = h*proj%rho0 - h*rho*cos(theta)

      ! Get i/j relative to reference i/j
      i = proj%knowni + (i - proj%polei)
      j = proj%knownj + (j - proj%polej)

      RETURN

   END SUBROUTINE llij_albers_nad83


   SUBROUTINE ijll_albers_nad83(i, j, proj, lat, lon)

      ! This is the inverse subroutine of llij_albers_nad83.  It returns the 
      ! latitude and longitude of an i/j point given the projection info 
      ! structure.  

      implicit none

      ! Arguments
      REAL(r8), INTENT(IN)                    :: i    ! Column
      REAL(r8), INTENT(IN)                    :: j    ! Row
      REAL(r8), INTENT(OUT)                   :: lat     ! -90 -> 90 north
      REAL(r8), INTENT(OUT)                   :: lon     ! -180 -> 180 East
      TYPE (proj_info), INTENT(IN)        :: proj

      ! Local variables
      real(r8) :: h, q, rho, theta, beta, x, y
      real(r8) :: a, b, c

      h = proj%hemi

      x = (i - proj%knowni + proj%polei)
      y = (j - proj%knownj + proj%polej)

      rho = sqrt(x**2.0 + (proj%rho0 - y)**2.0)
      theta = atan2(x, proj%rho0-y)

      q = (proj%bigc - (rho*proj%nc*proj%dx/A_NAD83)**2.0) / proj%nc

      beta = asin(q/(1.0_r8 - &
                  log((1.0_r8-E_NAD83)/(1.0_r8+E_NAD83))*(1.0_r8-E_NAD83**2.0)/(2.0_r8*E_NAD83)))
      a = (1.0_r8/3.0_r8)*E_NAD83**2 + (31.0_r8/180.0_r8)*E_NAD83**4 + (517.0_r8/5040.0_r8)*E_NAD83**6
      b = (23.0_r8/360.0_r8)*E_NAD83**4 + (251.0_r8/3780.0_r8)*E_NAD83**6
      c = (761.0_r8/45360.0_r8)*E_NAD83**6

      lat = beta + a*sin(2.0_r8*beta) + b*sin(4.0_r8*beta) + c*sin(6.0_r8*beta)

      lat = h*lat*deg_per_rad
      lon = proj%stdlon + theta*deg_per_rad/proj%nc

      RETURN

   END SUBROUTINE ijll_albers_nad83


   SUBROUTINE set_lc(proj)
      ! Initialize the remaining items in the proj structure for a
      ! lambert conformal grid.

      IMPLICIT NONE

      TYPE(proj_info), INTENT(INOUT)     :: proj

      REAL(r8)                               :: arg
      REAL(r8)                               :: deltalon1
      REAL(r8)                               :: tl1r
      REAL(r8)                               :: ctl1r

      ! Compute cone factor
      CALL lc_cone(proj%truelat1, proj%truelat2, proj%cone)

      ! Compute longitude differences and ensure we stay out of the
      ! forbidden "cut zone"
      deltalon1 = proj%lon1 - proj%stdlon
      IF (deltalon1 .GT. +180.0_r8) deltalon1 = deltalon1 - 360.0_r8
      IF (deltalon1 .LT. -180.0_r8) deltalon1 = deltalon1 + 360.0_r8

      ! Convert truelat1 to radian and compute COS for later use
      tl1r = proj%truelat1 * rad_per_deg
      ctl1r = COS(tl1r)

      ! Compute the radius to our known lower-left (SW) corner
      proj%rsw = proj%rebydx * ctl1r/proj%cone * &
             (TAN((90.0_r8*proj%hemi-proj%lat1)*rad_per_deg/2.0_r8) / &
              TAN((90.0_r8*proj%hemi-proj%truelat1)*rad_per_deg/2.0_r8))**proj%cone

      ! Find pole point
      arg = proj%cone*(deltalon1*rad_per_deg)
      proj%polei = proj%hemi*proj%knowni - proj%hemi * proj%rsw * SIN(arg)
      proj%polej = proj%hemi*proj%knownj + proj%rsw * COS(arg)  

      RETURN

   END SUBROUTINE set_lc                             


   SUBROUTINE lc_cone(truelat1, truelat2, cone)

   ! Subroutine to compute the cone factor of a Lambert Conformal projection

      IMPLICIT NONE

      ! Input Args
      REAL(r8), INTENT(IN)             :: truelat1  ! (-90 -> 90 degrees N)
      REAL(r8), INTENT(IN)             :: truelat2  !   "   "  "   "     "

      ! Output Args
      REAL(r8), INTENT(OUT)            :: cone

      ! Locals

      ! BEGIN CODE

      ! First, see if this is a secant or tangent projection.  For tangent
      ! projections, truelat1 = truelat2 and the cone is tangent to the 
      ! Earth's surface at this latitude.  For secant projections, the cone
      ! intersects the Earth's surface at each of the distinctly different
      ! latitudes
      IF (ABS(truelat1-truelat2) .GT. 0.10_r8) THEN
         cone = LOG10(COS(truelat1*rad_per_deg)) - &
                LOG10(COS(truelat2*rad_per_deg))
         cone = cone /(LOG10(TAN((45.0_r8 - ABS(truelat1)/2.0_r8) * rad_per_deg)) - &
                LOG10(TAN((45.0_r8 - ABS(truelat2)/2.0_r8) * rad_per_deg)))        
      ELSE
         cone = SIN(ABS(truelat1)*rad_per_deg )  
      ENDIF

      RETURN

   END SUBROUTINE lc_cone


   SUBROUTINE ijll_lc( i, j, proj, lat, lon)

   ! Subroutine to convert from the (i,j) cartesian coordinate to the 
   ! geographical latitude and longitude for a Lambert Conformal projection.

   ! History:
   ! 25 Jul 01: Corrected by B. Shaw, NOAA/FSL
   ! 
      IMPLICIT NONE

      ! Input Args
      REAL(r8), INTENT(IN)              :: i        ! Cartesian X coordinate
      REAL(r8), INTENT(IN)              :: j        ! Cartesian Y coordinate
      TYPE(proj_info),INTENT(IN)    :: proj     ! Projection info structure

      ! Output Args                 
      REAL(r8), INTENT(OUT)             :: lat      ! Latitude (-90->90 deg N)
      REAL(r8), INTENT(OUT)             :: lon      ! Longitude (-180->180 E)

      ! Locals 
      REAL(r8)                          :: inew
      REAL(r8)                          :: jnew
      REAL(r8)                          :: r
      REAL(r8)                          :: chi,chi1,chi2
      REAL(r8)                          :: r2
      REAL(r8)                          :: xx
      REAL(r8)                          :: yy

      ! BEGIN CODE

      chi1 = (90.0_r8 - proj%hemi*proj%truelat1)*rad_per_deg
      chi2 = (90.0_r8 - proj%hemi*proj%truelat2)*rad_per_deg

      ! See if we are in the southern hemispere and flip the indices
      ! if we are. 
      inew = proj%hemi * i
      jnew = proj%hemi * j

      ! Compute radius**2 to i/j location
      xx = inew - proj%polei
      yy = proj%polej - jnew
      r2 = (xx*xx + yy*yy)
      r = SQRT(r2)/proj%rebydx

      ! Convert to lat/lon
      IF (r2 .EQ. 0.0_r8) THEN
         lat = proj%hemi * 90.0_r8
         lon = proj%stdlon
      ELSE

         ! Longitude
         lon = proj%stdlon + deg_per_rad * ATAN2(proj%hemi*xx,yy)/proj%cone
         lon = MOD(lon+360.0_r8, 360.0_r8)

         ! Latitude.  Latitude determined by solving an equation adapted 
         ! from:
         !  Maling, D.H., 1973: Coordinate Systems and Map Projections
         ! Equations #20 in Appendix I.  

         IF (chi1 .EQ. chi2) THEN
            chi = 2.0_r8*ATAN( ( r/TAN(chi1) )**(1.0_r8/proj%cone) * TAN(chi1*0.5_r8) )
         ELSE
            chi = 2.0_r8*ATAN( (r*proj%cone/SIN(chi1))**(1.0_r8/proj%cone) * TAN(chi1*0.5_r8)) 
         ENDIF
         lat = (90.0_r8-chi*deg_per_rad)*proj%hemi

      ENDIF

      IF (lon .GT. +180.0_r8) lon = lon - 360.0_r8
      IF (lon .LT. -180.0_r8) lon = lon + 360.0_r8

      RETURN

   END SUBROUTINE ijll_lc


   SUBROUTINE llij_lc( lat, lon, proj, i, j)

   ! Subroutine to compute the geographical latitude and longitude values
   ! to the cartesian x/y on a Lambert Conformal projection.

      IMPLICIT NONE

      ! Input Args
      REAL(r8), INTENT(IN)              :: lat      ! Latitude (-90->90 deg N)
      REAL(r8), INTENT(IN)              :: lon      ! Longitude (-180->180 E)
      TYPE(proj_info),INTENT(IN)      :: proj     ! Projection info structure

      ! Output Args                 
      REAL(r8), INTENT(OUT)             :: i        ! Cartesian X coordinate
      REAL(r8), INTENT(OUT)             :: j        ! Cartesian Y coordinate

      ! Locals 
      REAL(r8)                          :: arg
      REAL(r8)                          :: deltalon
      REAL(r8)                          :: tl1r
      REAL(r8)                          :: rm
      REAL(r8)                          :: ctl1r


      ! BEGIN CODE

      ! Compute deltalon between known longitude and standard lon and ensure
      ! it is not in the cut zone
      deltalon = lon - proj%stdlon
      IF (deltalon .GT. +180.0_r8) deltalon = deltalon - 360.0_r8
      IF (deltalon .LT. -180.0_r8) deltalon = deltalon + 360.0_r8

      ! Convert truelat1 to radian and compute COS for later use
      tl1r = proj%truelat1 * rad_per_deg
      ctl1r = COS(tl1r)     

      ! Radius to desired point
      rm = proj%rebydx * ctl1r/proj%cone * &
           (TAN((90.0_r8*proj%hemi-lat)*rad_per_deg/2.0_r8) / &
            TAN((90.0_r8*proj%hemi-proj%truelat1)*rad_per_deg/2.0_r8))**proj%cone

      arg = proj%cone*(deltalon*rad_per_deg)
      i = proj%polei + proj%hemi * rm * SIN(arg)
      j = proj%polej - rm * COS(arg)

      ! Finally, if we are in the southern hemisphere, flip the i/j
      ! values to a coordinate system where (1,1) is the SW corner
      ! (what we assume) which is different than the original NCEP
      ! algorithms which used the NE corner as the origin in the 
      ! southern hemisphere (left-hand vs. right-hand coordinate?)
      i = proj%hemi * i  
      j = proj%hemi * j

      RETURN
   END SUBROUTINE llij_lc


   SUBROUTINE set_merc(proj)

      ! Sets up the remaining basic elements for the mercator projection

      IMPLICIT NONE
      TYPE(proj_info), INTENT(INOUT)       :: proj
      REAL(r8)                                 :: clain


      !  Preliminary variables

      clain = COS(rad_per_deg*proj%truelat1)
      proj%dlon = proj%dx / (proj%re_m * clain)

      ! Compute distance from equator to origin, and store in the 
      ! proj%rsw tag.

      proj%rsw = 0.0_r8
      IF (proj%lat1 .NE. 0.0_r8) THEN
         proj%rsw = (LOG(TAN(0.50_r8*((proj%lat1+90.0_r8)*rad_per_deg))))/proj%dlon
      ENDIF

      RETURN

   END SUBROUTINE set_merc


   SUBROUTINE llij_merc(lat, lon, proj, i, j)

      ! Compute i/j coordinate from lat lon for mercator projection

      IMPLICIT NONE
      REAL(r8), INTENT(IN)              :: lat
      REAL(r8), INTENT(IN)              :: lon
      TYPE(proj_info),INTENT(IN)    :: proj
      REAL(r8),INTENT(OUT)              :: i
      REAL(r8),INTENT(OUT)              :: j
      REAL(r8)                          :: deltalon, i2

      deltalon = lon - proj%lon1
      IF (deltalon .LT. -180.0_r8) deltalon = deltalon + 360.0_r8
      IF (deltalon .GT. 180.0_r8) deltalon = deltalon - 360.0_r8
      i = proj%knowni + (deltalon/(proj%dlon*deg_per_rad))
      i2 = proj%knowni + ((deltalon+360.0_r8)/(proj%dlon*deg_per_rad))
      if ( i < 0.0_r8 ) i = i2
      j = proj%knownj + (LOG(TAN(0.50_r8*((lat + 90.0_r8) * rad_per_deg)))) / &
             proj%dlon - proj%rsw

      RETURN

   END SUBROUTINE llij_merc


   SUBROUTINE ijll_merc(i, j, proj, lat, lon)

      ! Compute the lat/lon from i/j for mercator projection

      IMPLICIT NONE
      REAL(r8),INTENT(IN)               :: i
      REAL(r8),INTENT(IN)               :: j    
      TYPE(proj_info),INTENT(IN)    :: proj
      REAL(r8), INTENT(OUT)             :: lat
      REAL(r8), INTENT(OUT)             :: lon 


      lat = 2.00_r8*ATAN(EXP(proj%dlon*(proj%rsw + j-proj%knownj)))*deg_per_rad - 90.0_r8
      lon = (i-proj%knowni)*proj%dlon*deg_per_rad + proj%lon1
      IF (lon.GT.180.0_r8) lon = lon - 360.0_r8
      IF (lon.LT.-180.0_r8) lon = lon + 360.0_r8
      RETURN

   END SUBROUTINE ijll_merc


   SUBROUTINE llij_latlon(lat, lon, proj, i, j)

      ! Compute the i/j location of a lat/lon on a LATLON grid.
      IMPLICIT NONE
      REAL(r8), INTENT(IN)             :: lat
      REAL(r8), INTENT(IN)             :: lon
      TYPE(proj_info), INTENT(IN)  :: proj
      REAL(r8), INTENT(OUT)            :: i
      REAL(r8), INTENT(OUT)            :: j

      REAL(r8)                         :: deltalat
      REAL(r8)                         :: deltalon

      ! Compute deltalat and deltalon as the difference between the input 
      ! lat/lon and the origin lat/lon
      deltalat = lat - proj%lat1
      deltalon = lon - proj%lon1      

      ! Compute i/j
      i = deltalon/proj%loninc
      j = deltalat/proj%latinc

      i = i + proj%knowni
      j = j + proj%knownj

      RETURN

   END SUBROUTINE llij_latlon


   SUBROUTINE ijll_latlon(i, j, proj, lat, lon)

      ! Compute the lat/lon location of an i/j on a LATLON grid.
      IMPLICIT NONE
      REAL(r8), INTENT(IN)             :: i
      REAL(r8), INTENT(IN)             :: j
      TYPE(proj_info), INTENT(IN)  :: proj
      REAL(r8), INTENT(OUT)            :: lat
      REAL(r8), INTENT(OUT)            :: lon

      REAL(r8)                         :: i_work, j_work
      REAL(r8)                         :: deltalat
      REAL(r8)                         :: deltalon

      i_work = i - proj%knowni
      j_work = j - proj%knownj

      ! Compute deltalat and deltalon 
      deltalat = j_work*proj%latinc
      deltalon = i_work*proj%loninc

      lat = proj%lat1 + deltalat
      lon = proj%lon1 + deltalon

      RETURN

   END SUBROUTINE ijll_latlon


   SUBROUTINE set_cyl(proj)

      implicit none

      ! Arguments
      type(proj_info), intent(inout) :: proj

      proj%hemi = 1.00_r8

   END SUBROUTINE set_cyl


   SUBROUTINE llij_cyl(lat, lon, proj, i, j)

      implicit none

      ! Arguments
      real(r8), intent(in) :: lat, lon
      real(r8), intent(out) :: i, j
      type(proj_info), intent(in) :: proj

      ! Local variables
      real(r8) :: deltalat
      real(r8) :: deltalon

      ! Compute deltalat and deltalon as the difference between the input
      ! lat/lon and the origin lat/lon
      deltalat = lat - proj%lat1
!      deltalon = lon - proj%stdlon
      deltalon = lon - proj%lon1

      if (deltalon <   0.0_r8) deltalon = deltalon + 360.0_r8
      if (deltalon > 360.0_r8) deltalon = deltalon - 360.0_r8

      ! Compute i/j
      i = deltalon/proj%loninc
      j = deltalat/proj%latinc

      if (i <= 0.0_r8)              i = i + 360.0_r8/proj%loninc
      if (i > 360.0_r8/proj%loninc) i = i - 360.0_r8/proj%loninc

      i = i + proj%knowni
!nc -- typo in original I am assuming 
!nc -- orig -->     j = j + proj%knowni
      j = j + proj%knownj

   END SUBROUTINE llij_cyl


   SUBROUTINE ijll_cyl(i, j, proj, lat, lon)

      implicit none

      ! Arguments
      real(r8), intent(in) :: i, j
      real(r8), intent(out) :: lat, lon
      type(proj_info), intent(in) :: proj

      ! Local variables
      real(r8) :: deltalat
      real(r8) :: deltalon
      real(r8) :: i_work, j_work

      i_work = i - proj%knowni 
      j_work = j - proj%knownj

      if (i_work < 0.0_r8)                i_work = i_work + 360.0_r8/proj%loninc
      if (i_work >= 360.0_r8/proj%loninc) i_work = i_work - 360.0_r8/proj%loninc

      ! Compute deltalat and deltalon
      deltalat = j_work*proj%latinc
      deltalon = i_work*proj%loninc

      lat = deltalat + proj%lat1
!      lon = deltalon + proj%stdlon
      lon = deltalon + proj%lon1

      if (lon < -180.0_r8) lon = lon + 360.0_r8
      if (lon >  180.0_r8) lon = lon - 360.0_r8

   END SUBROUTINE ijll_cyl


   SUBROUTINE set_cassini(proj)

      implicit none

      ! Arguments
      type(proj_info), intent(inout) :: proj

      ! Local variables
      real(r8) :: comp_lat, comp_lon
      logical :: global_domain

      proj%hemi = 1.00_r8

      ! Try to determine whether this domain has global coverage
      if (abs(proj%lat1 - proj%latinc/2.0_r8 + 90.0_r8) < 0.001_r8 .and. &
          abs(mod(proj%lon1 - proj%loninc/2.0_r8 - proj%stdlon,360.0_r8)) < 0.001_r8) then
         global_domain = .true.
      else
         global_domain = .false.
      end if

      if (abs(proj%lat0) /= 90.0_r8 .and. .not.global_domain) then
         call rotate_coords(proj%lat1,proj%lon1,comp_lat,comp_lon,proj%lat0,proj%lon0,proj%stdlon,-1)
         proj%lat1 = comp_lat
         proj%lon1 = comp_lon
      end if

   END SUBROUTINE set_cassini


   SUBROUTINE llij_cassini(lat, lon, proj, i, j)

      implicit none

      ! Arguments
      real(r8), intent(in) :: lat, lon
      real(r8), intent(out) :: i, j
      type(proj_info), intent(in) :: proj

      ! Local variables
      real(r8) :: comp_lat, comp_lon

      ! Convert geographic to computational lat/lon
      if (abs(proj%lat0) /= 90.0_r8) then
         call rotate_coords(lat,lon,comp_lat,comp_lon,proj%lat0,proj%lon0,proj%stdlon,-1)
      else
         comp_lat = lat
         comp_lon = lon
      end if

      ! Convert computational lat/lon to i/j
      call llij_cyl(comp_lat, comp_lon, proj, i, j)

   END SUBROUTINE llij_cassini


   SUBROUTINE ijll_cassini(i, j, proj, lat, lon)

      implicit none

      ! Arguments
      real(r8), intent(in) :: i, j
      real(r8), intent(out) :: lat, lon
      type(proj_info), intent(in) :: proj

      ! Local variables
      real(r8) :: comp_lat, comp_lon

      ! Convert i/j to computational lat/lon
      call ijll_cyl(i, j, proj, comp_lat, comp_lon)

      ! Convert computational to geographic lat/lon
      if (abs(proj%lat0) /= 90.0_r8) then
         call rotate_coords(comp_lat,comp_lon,lat,lon,proj%lat0,proj%lon0,proj%stdlon,1)
      else
         lat = comp_lat
         lon = comp_lon
      end if

   END SUBROUTINE ijll_cassini


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Purpose: Converts between computational and geographic lat/lon for Cassini
   !          
   ! Notes: This routine was provided by Bill Skamarock, 2007-03-27
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   SUBROUTINE rotate_coords(ilat,ilon,olat,olon,lat_np,lon_np,lon_0,direction)

      IMPLICIT NONE

      REAL(r8), INTENT(IN   ) :: ilat, ilon
      REAL(r8), INTENT(  OUT) :: olat, olon
      REAL(r8), INTENT(IN   ) :: lat_np, lon_np, lon_0
      INTEGER, INTENT(IN  ), OPTIONAL :: direction
      ! >=0, default : computational -> geographical
      ! < 0          : geographical  -> computational

      REAL(r8) :: rlat, rlon
      REAL(r8) :: phi_np, lam_np, lam_0, dlam
      REAL(r8) :: sinphi, cosphi, coslam, sinlam

      ! Convert all angles to radians
      phi_np = lat_np * rad_per_deg
      lam_np = lon_np * rad_per_deg
      lam_0  = lon_0  * rad_per_deg
      rlat = ilat * rad_per_deg
      rlon = ilon * rad_per_deg

      IF (PRESENT(direction) .AND. (direction < 0)) THEN
         ! The equations are exactly the same except for one small difference
         ! with respect to longitude ...
         dlam = PI - lam_0
      ELSE
         dlam = lam_np
      END IF
      sinphi = COS(phi_np)*COS(rlat)*COS(rlon-dlam) + SIN(phi_np)*SIN(rlat)
      cosphi = SQRT(1.0_r8-sinphi*sinphi)
      coslam = SIN(phi_np)*COS(rlat)*COS(rlon-dlam) - COS(phi_np)*SIN(rlat)
      sinlam = COS(rlat)*SIN(rlon-dlam)
      IF ( cosphi /= 0.0_r8 ) THEN
         coslam = coslam/cosphi
         sinlam = sinlam/cosphi
      END IF
      olat = deg_per_rad*ASIN(sinphi)
      olon = deg_per_rad*(ATAN2(sinlam,coslam)-dlam-lam_0+lam_np)
      ! Both of my F90 text books prefer the DO-EXIT form, and claim it is faster
      ! when optimization is turned on (as we will always do...)
      DO
         IF (olon >= -180.0_r8) EXIT
         olon = olon + 360.0_r8
      END DO
      DO
         IF (olon <=  180.0_r8) EXIT
         olon = olon - 360.0_r8
      END DO

   END SUBROUTINE rotate_coords


   SUBROUTINE llij_rotlatlon(lat, lon, proj, i, j)

      IMPLICIT NONE

      ! Arguments
      REAL(r8), INTENT(IN) :: lat, lon
      REAL(r8), INTENT(OUT) :: i, j
      TYPE (proj_info), INTENT(IN) :: proj

      ! Local variables -- pi is a parameter in constants_module! -- local value will 
      !   be "ppi"
      REAL(KIND=HIGH) :: dphd,dlmd !Grid increments, degrees
      INTEGER :: ii,jj,jmt,ncol,nrow
      REAL(KIND=HIGH) :: glatd  !Geographic latitude, positive north
      REAL(KIND=HIGH) :: glond  !Geographic longitude, positive west
      REAL(KIND=HIGH) :: col,d1,d2,d2r,dlm,dlm1,dlm2,dph,glat,glon,    &
              ppi,r2d,row,tlat,tlat1,tlat2,              &
              tlon,tlon1,tlon2,tph0,tlm0,x,y,z

      glatd = REAL(lat,HIGH)                                        ! HIGH
      glond = -REAL(lon,HIGH)                                       ! HIGH

      dphd = REAL(proj%phi,HIGH)/REAL((proj%jydim-1)/2,HIGH)        ! HIGH
      dlmd = REAL(proj%lambda,HIGH)/REAL(proj%ixdim-1,HIGH)         ! HIGH

      ppi = ACOS(-1.0_HIGH)                                          ! HIGH (same name as module variable)
      d2r = ppi/180.0_HIGH                                           ! HIGH
      r2d = 1.0_HIGH/d2r                                            ! HIGH

      jmt = proj%jydim/2+1                                          ! INTEGER

      glat = glatd*d2r                                              ! HIGH
      glon = glond*d2r                                              ! HIGH
      dph = dphd*d2r                                                ! HIGH
      dlm = dlmd*d2r                                                ! HIGH
      tph0 = REAL(proj%lat1,HIGH)*d2r                               ! HIGH
      tlm0 = -REAL(proj%lon1,HIGH)*d2r                              ! HIGH

      x = COS(tph0)*COS(glat)*COS(glon-tlm0)+SIN(tph0)*SIN(glat)    ! HIGH
      y = -COS(glat)*SIN(glon-tlm0)                                 ! HIGH
      z = COS(tph0)*SIN(glat)-SIN(tph0)*COS(glat)*COS(glon-tlm0)    ! HIGH
      tlat = r2d*ATAN(z/SQRT(x*x+y*y))                              ! HIGH
      tlon = r2d*ATAN(y/x)                                          ! HIGH

      row = tlat/dphd+REAL(jmt,HIGH)                                ! HIGH
      col = tlon/dlmd+REAL(proj%ixdim,HIGH)                         ! HIGH

      nrow = NINT(row)                                              ! INTEGER
      ncol = NINT(col)                                              ! INTEGER

      tlat = tlat*d2r                                               ! HIGH
      tlon = tlon*d2r                                               ! HIGH

      IF (proj%stagger == HH) THEN

         IF ((abs(MOD(nrow,2)) == 1 .AND. abs(MOD(ncol,2)) == 1) .OR. &
             (MOD(nrow,2) == 0 .AND. MOD(ncol,2) == 0)) THEN

            tlat1 = REAL((nrow-jmt),HIGH)*dph
            tlat2 = tlat1+dph
            tlon1 = REAL((ncol-proj%ixdim),HIGH)*dlm
            tlon2 = tlon1+dlm

            dlm1 = tlon-tlon1
            dlm2 = tlon-tlon2
            d1 = ACOS(COS(tlat)*COS(tlat1)*COS(dlm1)+SIN(tlat)*SIN(tlat1))
            d2 = ACOS(COS(tlat)*COS(tlat2)*COS(dlm2)+SIN(tlat)*SIN(tlat2))

            IF (d1 > d2) THEN
               nrow = nrow+1
               ncol = ncol+1
            END IF

         ELSE
            tlat1 = REAL((nrow+1-jmt),HIGH)*dph
            tlat2 = tlat1-dph
            tlon1 = REAL((ncol-proj%ixdim),HIGH)*dlm
            tlon2 = tlon1+dlm
            dlm1 = tlon-tlon1
            dlm2 = tlon-tlon2
            d1 = ACOS(COS(tlat)*COS(tlat1)*COS(dlm1)+SIN(tlat)*SIN(tlat1))
            d2 = ACOS(COS(tlat)*COS(tlat2)*COS(dlm2)+SIN(tlat)*SIN(tlat2))

            IF (d1 < d2) THEN
               nrow = nrow+1
            ELSE
               ncol = ncol+1
            END IF
         END IF

      ELSE IF (proj%stagger == VV) THEN

         IF ((MOD(nrow,2) == 0 .AND. abs(MOD(ncol,2)) == 1) .OR. &
             (abs(MOD(nrow,2)) == 1 .AND. MOD(ncol,2) == 0)) THEN
            tlat1 = REAL((nrow-jmt),HIGH)*dph
            tlat2 = tlat1+dph
            tlon1 = REAL((ncol-proj%ixdim),HIGH)*dlm
            tlon2 = tlon1+dlm
            dlm1 = tlon-tlon1
            dlm2 = tlon-tlon2
            d1 = ACOS(COS(tlat)*COS(tlat1)*COS(dlm1)+SIN(tlat)*SIN(tlat1))
            d2 = ACOS(COS(tlat)*COS(tlat2)*COS(dlm2)+SIN(tlat)*SIN(tlat2))

            IF (d1 > d2) THEN
               nrow = nrow+1
               ncol = ncol+1
            END IF

         ELSE
            tlat1 = REAL((nrow+1-jmt),HIGH)*dph
            tlat2 = tlat1-dph
            tlon1 = REAL((ncol-proj%ixdim),HIGH)*dlm
            tlon2 = tlon1+dlm
            dlm1 = tlon-tlon1
            dlm2 = tlon-tlon2
            d1 = ACOS(COS(tlat)*COS(tlat1)*COS(dlm1)+SIN(tlat)*SIN(tlat1))
            d2 = ACOS(COS(tlat)*COS(tlat2)*COS(dlm2)+SIN(tlat)*SIN(tlat2))

            IF (d1 < d2) THEN
               nrow = nrow+1
            ELSE
               ncol = ncol+1
            END IF
         END IF
      END IF


!!! Added next line as a Kludge - not yet understood why needed
      if (ncol .le. 0) ncol=ncol-1

      jj = nrow
      ii = ncol/2

      IF (proj%stagger == HH) THEN
         IF (abs(MOD(jj,2)) == 1) ii = ii+1
      ELSE IF (proj%stagger == VV) THEN
         IF (MOD(jj,2) == 0) ii=ii+1
      END IF

      i = REAL(ii,r8)
      j = REAL(jj,r8)

   END SUBROUTINE llij_rotlatlon


   SUBROUTINE ijll_rotlatlon(i, j, proj, lat,lon)

      IMPLICIT NONE

      ! Arguments
      REAL(R8), INTENT(IN) :: i, j
      REAL(R8), INTENT(OUT) :: lat, lon
      TYPE (proj_info), INTENT(IN) :: proj

      ! Local variables
      INTEGER :: ih,jh
      INTEGER :: midcol,midrow,ncol
      REAL(KIND=HIGH) :: dphd,dlmd !Grid increments, degrees
      REAL(KIND=HIGH) :: arg1,arg2,d2r,fctr,glatr,glatd,glond,ppi, &
              r2d,tlatd,tlond,tlatr,tlonr,tlm0,tph0

      ih = NINT(i)
      jh = NINT(j)

      dphd = REAL(proj%phi,HIGH)/REAL((proj%jydim-1)/2,HIGH)
      dlmd = REAL(proj%lambda,HIGH)/REAL(proj%ixdim-1,HIGH)

      ppi = ACOS(-1.0_HIGH)
      d2r = ppi/180.0_HIGH
      r2d = 1.0_HIGH/d2r
      tph0 = REAL(proj%lat1,HIGH)*d2r
      tlm0 = -REAL(proj%lon1,HIGH)*d2r

      midrow = (proj%jydim+1)/2
      midcol = proj%ixdim

!      IF (proj%stagger == HH) THEN

!!!         ncol = 2*ih-1+MOD(jh+1,2)
         ncol = 2*ih-1+abs(MOD(jh+1,2))
         tlatd = REAL((jh-midrow),HIGH)*dphd
         tlond = REAL((ncol-midcol),HIGH)*dlmd
       IF (proj%stagger == VV) THEN
          if (mod(jh,2) .eq. 0) then
             tlond = tlond - DLMD
          else
             tlond = tlond + DLMD
          end if
       END IF

      tlatr = tlatd*d2r
      tlonr = tlond*d2r
      arg1 = SIN(tlatr)*COS(tph0)+COS(tlatr)*SIN(tph0)*COS(tlonr)
      glatr = ASIN(arg1)

      glatd = glatr*r2d

      arg2 = COS(tlatr)*COS(tlonr)/(COS(glatr)*COS(tph0))-TAN(glatr)*TAN(tph0)
      IF (ABS(arg2) > 1.0_HIGH) arg2 = ABS(arg2)/arg2
      fctr = 1.0_HIGH
      IF (tlond > 0.0_HIGH) fctr = -1.0_HIGH

      glond = tlm0*r2d+fctr*ACOS(arg2)*r2d

      lat = REAL(glatd,r8)
      lon = -REAL(glond,r8)

      IF (lon .GT. +180.0_r8) lon = lon - 360.0_r8
      IF (lon .LT. -180.0_r8) lon = lon + 360.0_r8

   END SUBROUTINE ijll_rotlatlon


   SUBROUTINE set_gauss(proj)

      IMPLICIT NONE

      ! Argument
      type (proj_info), intent(inout) :: proj

      !  Initialize the array that will hold the Gaussian latitudes.

      IF ( ASSOCIATED( proj%gauss_lat ) ) THEN
         DEALLOCATE ( proj%gauss_lat )
      END IF

      !  Get the needed space for our array.

      ALLOCATE ( proj%gauss_lat(proj%nlat*2) )

      !  Compute the Gaussian latitudes.

      CALL gausll( proj%nlat*2 , proj%gauss_lat )

      !  Now, these could be upside down from what we want, so let's check.
      !  We take advantage of the equatorial symmetry to remove any sort of
      !  array re-ordering.

      IF ( ABS(proj%gauss_lat(1) - proj%lat1) .GT. 0.01_r8 ) THEN
         proj%gauss_lat = -1.0_r8 * proj%gauss_lat
      END IF

      !  Just a sanity check.

      IF ( ABS(proj%gauss_lat(1) - proj%lat1) .GT. 0.01_r8 ) THEN
         PRINT '(A)','Oops, something is not right with the Gaussian latitude computation.'
         PRINT '(A,F8.3,A)','The input data gave the starting latitude as ',proj%lat1,'.'
         PRINT '(A,F8.3,A)','This routine computed the starting latitude as +-',ABS(proj%gauss_lat(1)),'.'
         PRINT '(A,F8.3,A)','The difference is larger than 0.01 degrees, which is not expected.'
         write(6,*) 'ERROR: Gaussian_latitude_computation'
      END IF

   END SUBROUTINE set_gauss


   SUBROUTINE gausll ( nlat , lat_sp )

      IMPLICIT NONE

      INTEGER                            :: nlat , i
      REAL (KIND=HIGH) , PARAMETER       :: ppi = 3.141592653589793_HIGH
      REAL (KIND=HIGH) , DIMENSION(nlat) :: cosc , gwt , sinc , colat , wos2 , lat
      REAL(R8)             , DIMENSION(nlat) :: lat_sp

      CALL lggaus(nlat, cosc, gwt, sinc, colat, wos2)

      DO i = 1, nlat
         lat(i) = ACOS(sinc(i)) * 180.0_HIGH / ppi
         IF (i.gt.nlat/2) lat(i) = -lat(i)
      END DO

      lat_sp = REAL(lat,r8)

   END SUBROUTINE gausll


   SUBROUTINE lggaus( nlat, cosc, gwt, sinc, colat, wos2 )

      IMPLICIT NONE

      !  LGGAUS finds the Gaussian latitudes by finding the roots of the
      !  ordinary Legendre polynomial of degree NLAT using Newton's
      !  iteration method.

      !  On entry:
      integer NLAT ! the number of latitudes (degree of the polynomial)

      !  On exit: for each Gaussian latitude
      !     COSC   - cos(colatitude) or sin(latitude)
      !     GWT    - the Gaussian weights
      !     SINC   - sin(colatitude) or cos(latitude)
      !     COLAT  - the colatitudes in radians
      !     WOS2   - Gaussian weight over sin**2(colatitude)

      REAL (KIND=HIGH) , DIMENSION(nlat) :: cosc , gwt , sinc , colat  , wos2 
      REAL (KIND=HIGH) , PARAMETER       :: ppi = 3.141592653589793_HIGH

      !  Convergence criterion for iteration of cos latitude

      REAL(KIND=HIGH) , PARAMETER :: xlim  = 1.0E-14_HIGH

      INTEGER :: nzero, i, j
      REAL (KIND=HIGH) :: fi, fi1, a, b, g, gm, gp, gt, delta, c, d

      !  The number of zeros between pole and equator

      nzero = nlat/2

      !  Set first guess for cos(colat)

      DO i=1,nzero
         cosc(i) = SIN( (REAL(i,HIGH)-0.5_HIGH)*ppi/REAL(nlat,HIGH) + ppi*0.5_HIGH )
      END DO

      !  Constants for determining the derivative of the polynomial
      fi  = nlat
      fi1 = fi + 1.0_HIGH
      a   = fi*fi1 / SQRT(4.0_HIGH * fi1*fi1 - 1.0_HIGH)
      b   = fi1*fi / SQRT(4.0_HIGH * fi*fi - 1.0_HIGH)

      !  Loop over latitudes, iterating the search for each root

      DO i=1,nzero
         j=0

         !  Determine the value of the ordinary Legendre polynomial for
         !  the current guess root

         DO
            CALL lgord( g, cosc(i), nlat )

            !  Determine the derivative of the polynomial at this point

            CALL lgord( gm, cosc(i), nlat-1 )
            CALL lgord( gp, cosc(i), nlat+1 )
            gt = (cosc(i)*cosc(i)-1.0_HIGH) / (a*gp-b*gm)

            !  Update the estimate of the root

            delta   = g*gt
            cosc(i) = cosc(i) - delta

            !  If convergence criterion has not been met, keep trying

            j = j+1
            IF( ABS(delta).GT.xlim ) CYCLE

            !  Determine the Gaussian weights

            c      = 2.0_HIGH *( 1.0_HIGH - cosc(i)*cosc(i) )
            CALL lgord( d, cosc(i), nlat-1 )
            d      = d*d*fi*fi
            gwt(i) = c *( fi-0.5_HIGH ) / d
            EXIT

         END DO

      END DO

      !  Determine the colatitudes and sin(colat) and weights over sin**2

      DO i=1,nzero
         colat(i)= ACOS(cosc(i))
         sinc(i) = SIN(colat(i))
         wos2(i) = gwt(i) /( sinc(i)*sinc(i) )
      END DO

      !  If NLAT is odd, set values at the equator

      IF( MOD(nlat,2) .NE. 0 ) THEN
         i       = nzero+1
         cosc(i) = 0.0_HIGH
         c       = 2.0_HIGH
         CALL lgord( d, cosc(i), nlat-1 )
         d       = d*d*fi*fi
         gwt(i)  = c *( fi-0.5_HIGH ) / d
         colat(i)= ppi*0.5_HIGH
         sinc(i) = 1.0_HIGH
         wos2(i) = gwt(i)
      END IF

      !  Determine the southern hemisphere values by symmetry

      DO i=nlat-nzero+1,nlat
         cosc(i) =-cosc(nlat+1-i)
         gwt(i)  = gwt(nlat+1-i)
         colat(i)= ppi-colat(nlat+1-i)
         sinc(i) = sinc(nlat+1-i)
         wos2(i) = wos2(nlat+1-i)
      END DO

   END SUBROUTINE lggaus


   SUBROUTINE lgord( f, cosc, n )

      IMPLICIT NONE

      !  LGORD calculates the value of an ordinary Legendre polynomial at a
      !  specific latitude.

      !  On entry:
      !     cosc - COS(colatitude)
      !     n      - the degree of the polynomial

      !  On exit:
      !     f      - the value of the Legendre polynomial of degree N at
      !              latitude ASIN(cosc)

      REAL (KIND=HIGH) :: s1, c4, a, b, fk, f, cosc, colat, c1, fn, ang
      INTEGER :: n, k

      !  Determine the colatitude

      colat = ACOS(cosc)

      c1 = SQRT(2.0_HIGH)
      DO k=1,n
         c1 = c1 * SQRT( 1.0_HIGH - 1.0_HIGH/REAL((4*k*k),HIGH) )
      END DO

      fn = n
      ang= fn * colat
      s1 = 0.0_HIGH
      c4 = 1.0_HIGH
      a  =-1.0_HIGH
      b  = 0.0_HIGH
      DO k=0,n,2
         IF (k.eq.n) c4 = 0.5_HIGH * c4
         s1 = s1 + c4 * COS(ang)
         a  = a + 2.0_HIGH
         b  = b + 1.0_HIGH
         fk = k
         ang= colat * (fn-fk-2.0_HIGH)
         c4 = ( a * (fn-b+1.0_HIGH) / ( b * (fn+fn-a) ) ) * c4
      END DO 

      f = s1 * c1

   END SUBROUTINE lgord


   SUBROUTINE llij_gauss (lat, lon, proj, i, j) 

      IMPLICIT NONE

      REAL(R8)    , INTENT(IN)  :: lat, lon
      REAL(R8)    , INTENT(OUT) :: i, j
      TYPE (proj_info), INTENT(IN) :: proj

      INTEGER :: n , n_low
      LOGICAL :: found = .FALSE.
      REAL(R8)    :: diff_1 , diff_nlat

      !  The easy one first, get the x location.  The calling routine has already made
      !  sure that the necessary assumptions concerning the sign of the deltalon and the
      !  relative east/west'ness of the longitude and the starting longitude are consistent
      !  to allow this easy computation.

      i = ( lon - proj%lon1 ) / proj%loninc + 1.0_r8

      !  Since this is a global data set, we need to be concerned about wrapping the
      !  fields around the globe.

!      IF      ( ( proj%loninc .GT. 0 ) .AND. &
!                ( FLOOR((lon-proj%lon1)/proj%loninc) + 1 .GE. proj%ixdim ) .AND. &
!                ( lon + proj%loninc .GE. proj%lon1 + 360 ) ) THEN
!! BUG: We may need to set proj%wrap, but proj is intent(in)
!! WHAT IS THIS USED FOR?
!!        proj%wrap = .TRUE.
!         i = proj%ixdim
!      ELSE IF ( ( proj%loninc .LT. 0 ) .AND. &
!                ( FLOOR((lon-proj%lon1)/proj%loninc) + 1 .GE. proj%ixdim ) .AND. &
!                ( lon + proj%loninc .LE. proj%lon1 - 360 ) ) THEN
! ! BUG: We may need to set proj%wrap, but proj is intent(in)
! ! WHAT IS THIS USED FOR?
! !        proj%wrap = .TRUE.
!         i = proj%ixdim
!      END IF

      !  Yet another quicky test, can we find bounding values?  If not, then we may be
      !  dealing with putting data to a polar projection, so just give them them maximal
      !  value for the location.  This is an OK assumption for the interpolation across the
      !  top of the pole, given how close the longitude lines are.

      IF ( ABS(lat) .GT. ABS(proj%gauss_lat(1)) ) THEN

         diff_1    = lat - proj%gauss_lat(1)
         diff_nlat = lat - proj%gauss_lat(proj%nlat*2)

         IF ( ABS(diff_1) .LT. ABS(diff_nlat) ) THEN
            j = 1
         ELSE
            j = proj%nlat*2
         END IF

      !  If the latitude is between the two bounding values, we have to search and interpolate.

      ELSE

         DO n = 1 , proj%nlat*2 -1
            IF ( ( proj%gauss_lat(n) - lat ) * ( proj%gauss_lat(n+1) - lat ) .LE. 0.0_r8 ) THEN
               found = .TRUE.
               n_low = n
               EXIT
            END IF
         END DO

         !  Everything still OK?

         IF ( .NOT. found ) THEN
            PRINT '(A)','Troubles in river city.  No bounding values of latitude found in the Gaussian routines.'
            write(6,*) 'ERROR: Gee_no_bounding_lats_Gaussian'
         END IF

         j = ( ( proj%gauss_lat(n_low) - lat                     ) * REAL( (n_low + 1),r8 ) + &
               ( lat                   - proj%gauss_lat(n_low+1) ) * REAL( n_low,r8     ) ) / &
               ( proj%gauss_lat(n_low) - proj%gauss_lat(n_low+1) )

      END IF

   END SUBROUTINE llij_gauss 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!nc -- This subroutine was in the original module_map_utils.f90 (before PROJ_CASSINI), and
!   the model_mod.f90 would still like it to be, hence we are including it.
  SUBROUTINE gridwind_to_truewind(lon,proj,ugrid,vgrid,utrue,vtrue)
  
    ! Subroutine to convert a wind from grid north to true north.

    IMPLICIT NONE
    
    ! Arguments
    REAL(r8), INTENT(IN)       :: lon     ! Longitude of point in degrees
    TYPE(proj_info),INTENT(IN) :: proj    ! Projection info structure 
    REAL(r8), INTENT(IN)       :: ugrid   ! U-component, grid-relative
    REAL(r8), INTENT(IN)       :: vgrid   ! V-component, grid-relative
    REAL(r8), INTENT(OUT)      :: utrue   ! U-component, earth-relative
    REAL(r8), INTENT(OUT)      :: vtrue   ! V-component, earth-relative

    ! Locals
    REAL(r8)                   :: alpha
    REAL(r8)                   :: diff

    IF ((proj%code .EQ. PROJ_PS).OR.(proj%code .EQ. PROJ_LC))THEN
      diff = lon - proj%stdlon
      IF (diff .GT. 180.0_r8) diff = diff - 360.0_r8
      IF (diff .LT.-180.0_r8) diff = diff + 360.0_r8
      alpha = diff * proj%cone * rad_per_deg * SIGN(1.0_r8,proj%truelat1)
      utrue = vgrid * SIN(alpha) + ugrid * COS(alpha)
      vtrue = vgrid * COS(alpha) - ugrid * SIN(alpha)
!nc -- we added in a case structure for CASSINI and CYL
    ELSEIF ((proj%code .EQ. PROJ_MERC).OR.(proj%code .EQ. PROJ_LATLON) &
            .OR.(proj%code .EQ. PROJ_CASSINI).OR.(proj%code .EQ. PROJ_CYL))THEN
      utrue = ugrid
      vtrue = vgrid
    ELSE
      PRINT '(A)', 'Unrecognized projection.'
      STOP 'GRIDWIND_TO_TRUEWIND'
    ENDIF

    RETURN
  END SUBROUTINE gridwind_to_truewind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!nc -- This subroutine was in the original module_map_utils.f90 (before PROJ_CASSINI), and
!   the model_mod.f90 would still like it to be, hence we are including it.
  SUBROUTINE truewind_to_gridwind(lon,proj,utrue,vtrue,ugrid,vgrid)
      
    ! Subroutine to compute grid-relative u/v wind components from the earth-
    ! relative values for a given projection.

    IMPLICIT NONE

    ! Arguments
    REAL(r8), INTENT(IN)       :: lon     ! Longitude of point in degrees
    TYPE(proj_info),INTENT(IN) :: proj    ! Projection info structure
    REAL(r8), INTENT(IN)       :: utrue   ! U-component, earth-relative
    REAL(r8), INTENT(IN)       :: vtrue   ! V-component, earth-relative
    REAL(r8), INTENT(OUT)      :: ugrid   ! U-component, grid-relative
    REAL(r8), INTENT(OUT)      :: vgrid   ! V-component, grid-relative

    ! Locals
    REAL(r8)                   :: alpha
    REAL(r8)                   :: diff

    IF ((proj%code .EQ. PROJ_PS).OR.(proj%code .EQ. PROJ_LC))THEN

      diff = proj%stdlon - lon
      IF (diff .GT. 180.0_r8) diff = diff - 360.0_r8
      IF (diff .LT.-180.0_r8) diff = diff + 360.0_r8
      alpha = diff * proj%cone * rad_per_deg * SIGN(1.0_r8,proj%truelat1)
      ugrid = vtrue * SIN(alpha) + utrue * COS(alpha)
      vgrid = vtrue * COS(alpha) - utrue * SIN(alpha)
!nc -- we added in a case structure for CASSINI and CYL
    ELSEIF ((proj%code .EQ. PROJ_MERC).OR.(proj%code .EQ. PROJ_LATLON) &
            .OR.(proj%code .EQ. PROJ_CASSINI).OR.(proj%code .EQ. PROJ_CYL))THEN
      ugrid = utrue
      vgrid = vtrue
    ELSE
      PRINT '(A)', 'Unrecognized map projection.'
      STOP 'TRUEWIND_TO_GRIDWIND'
    ENDIF
    RETURN
  END SUBROUTINE truewind_to_gridwind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE map_utils

