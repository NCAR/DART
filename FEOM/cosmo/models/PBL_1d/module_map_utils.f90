!dis
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis

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
! Polar Stereographic (code = PROJ_PS)
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
!     For the latlon projection, the grid origin may be at any of the
!     corners, and the deltalat and deltalon values can be signed to 
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
!         code (input) = one of PROJ_LATLON, PROJ_MERC, PROJ_LC, or PROJ_PS
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
!  Hoke, Hayes, and Renninger, "Map Preojections and Grid Systems for
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
! DART $Id$
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$


  use         types_mod, only : r8

  IMPLICIT NONE
  private

  public proj_info, map_init, map_set, latlon_to_ij, &
       PROJ_LATLON, PROJ_MERC, PROJ_LC, PROJ_PS, gridwind_to_truewind

  ! Define some private constants
  REAL (kind=r8), PARAMETER    :: pi = 3.1415927_r8
  REAL (kind=r8), PARAMETER    :: deg_per_rad = 180.0_r8/pi
  REAL (kind=r8), PARAMETER    :: rad_per_deg = pi / 180.0_r8

  ! Mean Earth Radius in m.  The value below is consistent
  ! with NCEP's routines and grids.
! REAL (kind=r8), PARAMETER    :: earth_radius_m = 6371200.0_r8 ! from Brent
  REAL (kind=r8), PARAMETER    :: earth_radius_m = 6370000.0_r8 ! same as MM5 system
  REAL (kind=r8), PARAMETER    :: radians_per_degree = pi / 180.0_r8

  ! Projection codes for proj_info structure:
  INTEGER, PARAMETER  :: PROJ_LATLON = 0
  INTEGER, PARAMETER  :: PROJ_MERC = 1
  INTEGER, PARAMETER  :: PROJ_LC = 3
  INTEGER, PARAMETER  :: PROJ_PS = 5


  ! Define data structures to define various projections

  TYPE proj_info

     INTEGER          :: code     ! Integer code for projection type
     REAL(r8)         :: lat1     ! SW latitude (1,1) in degrees (-90->90N)
     REAL(r8)         :: lon1     ! SW longitude (1,1) in degrees (-180->180E)
     REAL(r8)         :: dx       ! Grid spacing in meters at truelats, used
                                  ! only for ps, lc, and merc projections
     REAL(r8)         :: dlat     ! Lat increment for lat/lon grids
     REAL(r8)         :: dlon     ! Lon increment for lat/lon grids
     REAL(r8)         :: stdlon   ! Longitude parallel to y-axis (-180->180E)
     REAL(r8)         :: truelat1 ! First true latitude (all projections)
     REAL(r8)         :: truelat2 ! Second true lat (LC only)
     REAL(r8)         :: hemi     ! 1 for NH, -1 for SH
     REAL(r8)         :: cone     ! Cone factor for LC projections
     REAL(r8)         :: polei    ! Computed i-location of pole point
     REAL(r8)         :: polej    ! Computed j-location of pole point
     REAL(r8)         :: rsw      ! Computed radius to SW corner
     REAL(r8)         :: rebydx   ! Earth radius divided by dx
     REAL(r8)         :: knowni   ! X-location of known lat/lon
     REAL(r8)         :: knownj   ! Y-location of known lat/lon
     LOGICAL          :: init     ! Flag to indicate if this struct is 
                                  ! ready for use
  END TYPE proj_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE map_init(proj)
    ! Initializes the map projection structure to missing values

    IMPLICIT NONE
    TYPE(proj_info), INTENT(INOUT)  :: proj

    proj%lat1 =    -999.9_r8
    proj%lon1 =    -999.9_r8
    proj%dx    =    -999.9_r8
    proj%stdlon =   -999.9_r8
    proj%truelat1 = -999.9_r8
    proj%truelat2 = -999.9_r8
    proj%hemi     = 0.0_r8
    proj%cone     = -999.9_r8
    proj%polei    = -999.9_r8
    proj%polej    = -999.9_r8
    proj%rsw      = -999.9_r8
    proj%knowni   = -999.9_r8
    proj%knownj   = -999.9_r8
    proj%init     = .FALSE.
  
  END SUBROUTINE map_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE map_set(proj_code,lat1,lon1,knowni,knownj,dx,stdlon,truelat1,truelat2,proj)
    ! Given a partially filled proj_info structure, this routine computes
    ! polei, polej, rsw, and cone (if LC projection) to complete the 
    ! structure.  This allows us to eliminate redundant calculations when
    ! calling the coordinate conversion routines multiple times for the
    ! same map.
    ! This will generally be the first routine called when a user wants
    ! to be able to use the coordinate conversion routines, and it
    ! will call the appropriate subroutines based on the 
    ! proj%code which indicates which projection type  this is.

    IMPLICIT NONE
    
    ! Declare arguments
    INTEGER, INTENT(IN)               :: proj_code
    REAL(r8), INTENT(IN)              :: lat1
    REAL(r8), INTENT(IN)              :: lon1
    REAL(r8), INTENT(IN)              :: dx
    REAL(r8), INTENT(IN)              :: stdlon
    REAL(r8), INTENT(IN)              :: truelat1
    REAL(r8), INTENT(IN)              :: truelat2
    REAL(r8), INTENT(IN)              :: knowni , knownj
    TYPE(proj_info), INTENT(OUT)      :: proj

    ! Local variables


    ! Executable code

    ! First, check for validity of mandatory variables in proj
    IF ( ABS(lat1) .GT. 90.0_r8 ) THEN
      PRINT '(A)', 'Latitude of origin corner required as follows:'
      PRINT '(A)', '    -90N <= lat1 < = 90.N'
      STOP 'MAP_INIT'
    ENDIF
    IF ( ABS(lon1) .GT. 180.0_r8) THEN
      PRINT '(A)', 'Longitude of origin required as follows:'
      PRINT '(A)', '   -180E <= lon1 <= 180W'
      STOP 'MAP_INIT'
    ENDIF
    IF ((dx .LE. 0.0_r8).AND.(proj_code .NE. PROJ_LATLON)) THEN
      PRINT '(A)', 'Require grid spacing (dx) in meters be positive!'
      STOP 'MAP_INIT'
    ENDIF
    IF ((ABS(stdlon) .GT. 180.0_r8).AND.(proj_code .NE. PROJ_MERC)) THEN
      PRINT '(A)', 'Need orientation longitude (stdlon) as: '
      PRINT '(A)', '   -180E <= lon1 <= 180W' 
      STOP 'MAP_INIT'
    ENDIF
    IF (ABS(truelat1).GT.90.0_r8) THEN
      PRINT '(A)', 'Set true latitude 1 for all projections!'
      STOP 'MAP_INIT'
    ENDIF
   
    CALL map_init(proj) 
    proj%code  = proj_code
    proj%lat1 = lat1
    proj%lon1 = lon1
    proj%knowni = knowni
    proj%knownj = knownj
    proj%dx    = dx
    proj%stdlon = stdlon
    proj%truelat1 = truelat1
    proj%truelat2 = truelat2
    IF (proj%code .NE. PROJ_LATLON) THEN
      proj%dx = dx
      IF (truelat1 .LT. 0.0_r8) THEN
        proj%hemi = -1.0_r8
      ELSE
        proj%hemi = 1.0_r8
      ENDIF
      proj%rebydx = earth_radius_m / dx
    ENDIF
    pick_proj: SELECT CASE(proj%code)

      CASE(PROJ_PS)
        PRINT '(A)', 'Setting up POLAR STEREOGRAPHIC map...'
        CALL set_ps(proj)

      CASE(PROJ_LC)
        PRINT '(A)', 'Setting up LAMBERT CONFORMAL map...'
        IF (ABS(proj%truelat2) .GT. 90.0_r8) THEN
          PRINT '(A)', 'Second true latitude not set, assuming a tangent'
          PRINT '(A,F10.3)', 'projection at truelat1: ', proj%truelat1
          proj%truelat2=proj%truelat1
        ENDIF
        CALL set_lc(proj)
   
      CASE (PROJ_MERC)
        PRINT '(A)', 'Setting up MERCATOR map...'
        CALL set_merc(proj)
   
      CASE (PROJ_LATLON)
        PRINT '(A)', 'Setting up CYLINDRICAL EQUIDISTANT LATLON map...'
        ! Convert lon1 to 0->360 notation
        IF (proj%lon1 .LT. 0.0_r8) proj%lon1 = proj%lon1 + 360.0_r8
   
      CASE DEFAULT
        PRINT '(A,I2)', 'Unknown projection code: ', proj%code
        STOP 'MAP_INIT'
    
    END SELECT pick_proj
    proj%init = .TRUE.
    RETURN
  END SUBROUTINE map_set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE latlon_to_ij(proj, lat, lon, i, j)
    ! Converts input lat/lon values to the cartesian (i,j) value
    ! for the given projection. 

    IMPLICIT NONE
    TYPE(proj_info), INTENT(IN)          :: proj
    REAL(r8), INTENT(IN)                 :: lat
    REAL(r8), INTENT(IN)                 :: lon
    REAL(r8), INTENT(OUT)                :: i
    REAL(r8), INTENT(OUT)                :: j

    IF (.NOT.proj%init) THEN
      PRINT '(A)', 'You have not called map_set for this projection!'
      STOP 'LATLON_TO_IJ'
    ENDIF

    SELECT CASE(proj%code)
 
      CASE(PROJ_LATLON)
        CALL llij_latlon(lat,lon,proj,i,j)

      CASE(PROJ_MERC)
        CALL llij_merc(lat,lon,proj,i,j)

      CASE(PROJ_PS)
        CALL llij_ps(lat,lon,proj,i,j)
      
      CASE(PROJ_LC)
        CALL llij_lc(lat,lon,proj,i,j)

      CASE DEFAULT
        PRINT '(A,I2)', 'Unrecognized map projection code: ', proj%code
        STOP 'LATLON_TO_IJ'
 
    END SELECT
    RETURN
  END SUBROUTINE latlon_to_ij
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ij_to_latlon(proj, i, j, lat, lon)
    ! Computes geographical latitude and longitude for a given (i,j) point
    ! in a grid with a projection of proj

    IMPLICIT NONE
    TYPE(proj_info),INTENT(IN)          :: proj
    REAL(r8), INTENT(IN)                :: i
    REAL(r8), INTENT(IN)                :: j
    REAL(r8), INTENT(OUT)               :: lat
    REAL(r8), INTENT(OUT)               :: lon

    IF (.NOT.proj%init) THEN
      PRINT '(A)', 'You have not called map_set for this projection!'
      STOP 'IJ_TO_LATLON'
    ENDIF
    SELECT CASE (proj%code)

      CASE (PROJ_LATLON)
        CALL ijll_latlon(i, j, proj, lat, lon)

      CASE (PROJ_MERC)
        CALL ijll_merc(i, j, proj, lat, lon)

      CASE (PROJ_PS)
        CALL ijll_ps(i, j, proj, lat, lon)

      CASE (PROJ_LC)
        CALL ijll_lc(i, j, proj, lat, lon)

      CASE DEFAULT
        PRINT '(A,I2)', 'Unrecognized map projection code: ', proj%code
        STOP 'IJ_TO_LATLON'

    END SELECT
    RETURN
  END SUBROUTINE ij_to_latlon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE set_ps(proj)
    ! Initializes a polar-stereographic map projection from the partially
    ! filled proj structure. This routine computes the radius to the
    ! southwest corner and computes the i/j location of the pole for use
    ! in llij_ps and ijll_ps.
    IMPLICIT NONE
 
    ! Declare args
    TYPE(proj_info), INTENT(INOUT)    :: proj

    ! Local vars
    REAL(r8)                          :: ala1
    REAL(r8)                          :: alo1
    REAL(r8)                          :: reflon
    REAL(r8)                          :: scale_top

    ! Executable code
    reflon = proj%stdlon + 90.0_r8

    ! Compute numerator term of map scale factor
    scale_top = 1.0_r8 + proj%hemi * SIN(proj%truelat1 * rad_per_deg)

    ! Compute radius to lower-left (SW) corner
    ala1 = proj%lat1 * rad_per_deg
    proj%rsw = proj%rebydx*COS(ala1)*scale_top/(1.0_r8+proj%hemi*SIN(ala1))

    ! Find the pole point
    alo1 = (proj%lon1 - reflon) * rad_per_deg
    proj%polei = proj%knowni - proj%rsw * COS(alo1)
    proj%polej = proj%knownj - proj%hemi * proj%rsw * SIN(alo1)
    PRINT '(A,2F10.1)', 'Computed (I,J) of pole point: ',proj%polei,proj%polej
    RETURN
  END SUBROUTINE set_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE llij_ps(lat,lon,proj,i,j)
    ! Given latitude (-90 to 90), longitude (-180 to 180), and the
    ! standard polar-stereographic projection information via the 
    ! proj structure, this routine returns the i/j indices which
    ! if within the domain range from 1->nx and 1->ny, respectively.

    IMPLICIT NONE

    ! Delcare input arguments
    REAL(r8), INTENT(IN)           :: lat
    REAL(r8), INTENT(IN)           :: lon
    TYPE(proj_info),INTENT(IN)     :: proj

    ! Declare output arguments     
    REAL(r8), INTENT(OUT)          :: i !(x-index)
    REAL(r8), INTENT(OUT)          :: j !(y-index)

    ! Declare local variables
    
    REAL(r8)                       :: reflon
    REAL(r8)                       :: scale_top
    REAL(r8)                       :: ala
    REAL(r8)                       :: alo
    REAL(r8)                       :: rm

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ijll_ps(i, j, proj, lat, lon)

    ! This is the inverse subroutine of llij_ps.  It returns the 
    ! latitude and longitude of an i/j point given the projection info 
    ! structure.  

    IMPLICIT NONE

    ! Declare input arguments
    REAL(r8), INTENT(IN)                :: i    ! Column
    REAL(r8), INTENT(IN)                :: j    ! Row
    TYPE (proj_info), INTENT(IN)        :: proj
    
    ! Declare output arguments
    REAL(r8), INTENT(OUT)               :: lat     ! -90 -> 90 North
    REAL(r8), INTENT(OUT)               :: lon     ! -180 -> 180 East

    ! Local variables
    REAL(r8)                            :: reflon
    REAL(r8)                            :: scale_top
    REAL(r8)                            :: xx,yy
    REAL(r8)                            :: gi2, r2
    REAL(r8)                            :: arccos

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
      gi2 = (proj%rebydx * scale_top)**2.0_r8
      lat = deg_per_rad * proj%hemi * ASIN((gi2-r2)/(gi2+r2))
      arccos = ACOS(xx/SQRT(r2))
      IF (yy .GT. 0) THEN
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE set_lc(proj)
    ! Initialize the remaining items in the proj structure for a
    ! lambert conformal grid.

    IMPLICIT NONE
    
    TYPE(proj_info), INTENT(INOUT)     :: proj

    REAL(r8)                           :: arg
    REAL(r8)                           :: deltalon1
    REAL(r8)                           :: tl1r
    REAL(r8)                           :: ctl1r

    ! Compute cone factor
    CALL lc_cone(proj%truelat1, proj%truelat2, proj%cone)
    PRINT '(A,F8.6)', 'Computed cone factor: ', proj%cone
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
    proj%polei = 1.0_r8 - proj%hemi * proj%rsw * SIN(arg)
    proj%polej = 1.0_r8 + proj%rsw * COS(arg)  
    PRINT '(A,2F10.3)', 'Computed pole i/j = ', proj%polei, proj%polej
    RETURN
  END SUBROUTINE set_lc                             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE lc_cone(truelat1, truelat2, cone)

  ! Subroutine to compute the cone factor of a Lambert Conformal projection

    IMPLICIT NONE
    
    ! Input Args
    REAL(r8), INTENT(IN)         :: truelat1  ! (-90 -> 90 degrees N)
    REAL(r8), INTENT(IN)         :: truelat2  !   "   "  "   "     "

    ! Output Args
    REAL(r8), INTENT(OUT)        :: cone

    ! Locals

    ! BEGIN CODE

    ! First, see if this is a secant or tangent projection.  For tangent
    ! projections, truelat1 = truelat2 and the cone is tangent to the 
    ! Earth's surface at this latitude.  For secant projections, the cone
    ! intersects the Earth's surface at each of the distinctly different
    ! latitudes
    IF (ABS(truelat1-truelat2) .GT. 0.1_r8) THEN
      cone = ALOG10(COS(truelat1*rad_per_deg)) - &
             ALOG10(COS(truelat2*rad_per_deg))
      cone = cone /(ALOG10(TAN((45.0_r8 - ABS(truelat1)/2.0_r8) * rad_per_deg)) - &
             ALOG10(TAN((45.0_r8 - ABS(truelat2)/2.0_r8) * rad_per_deg)))        
    ELSE
       cone = SIN(ABS(truelat1)*rad_per_deg )  
    ENDIF
  RETURN
  END SUBROUTINE lc_cone
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ijll_lc( i, j, proj, lat, lon)

  ! Subroutine to convert from the (i,j) cartesian coordinate to the 
  ! geographical latitude and longitude for a Lambert Conformal projection.

  ! History:
  ! 25 Jul 01: Corrected by B. Shaw, NOAA/FSL
  ! 
    IMPLICIT NONE

    ! Input Args
    REAL(r8), INTENT(IN)          :: i        ! Cartesian X coordinate
    REAL(r8), INTENT(IN)          :: j        ! Cartesian Y coordinate
    TYPE(proj_info),INTENT(IN)    :: proj     ! Projection info structure

    ! Output Args                 
    REAL(r8), INTENT(OUT)         :: lat      ! Latitude (-90->90 deg N)
    REAL(r8), INTENT(OUT)         :: lon      ! Longitude (-180->180 E)

    ! Locals 
    REAL(r8)                      :: inew
    REAL(r8)                      :: jnew
    REAL(r8)                      :: r
    REAL(r8)                      :: chi,chi1,chi2
    REAL(r8)                      :: r2
    REAL(r8)                      :: xx
    REAL(r8)                      :: yy

    ! BEGIN CODE

    chi1 = (90.0_r8 - proj%hemi*proj%truelat1)*rad_per_deg
    chi2 = (90.0_r8 - proj%hemi*proj%truelat2)*rad_per_deg

    ! See if we are in the southern hemispere and flip the indices
    ! if we are. 
    IF (proj%hemi .EQ. -1.0_r8) THEN 
      inew = -i + 2.0_r8
      jnew = -j + 2.0_r8
    ELSE
      inew = i
      jnew = j
    ENDIF

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE llij_lc( lat, lon, proj, i, j)

  ! Subroutine to compute the geographical latitude and longitude values
  ! to the cartesian x/y on a Lambert Conformal projection.
    
    IMPLICIT NONE

    ! Input Args
    REAL(r8), INTENT(IN)          :: lat      ! Latitude (-90->90 deg N)
    REAL(r8), INTENT(IN)          :: lon      ! Longitude (-180->180 E)
    TYPE(proj_info),INTENT(IN)    :: proj     ! Projection info structure

    ! Output Args                 
    REAL(r8), INTENT(OUT)         :: i        ! Cartesian X coordinate
    REAL(r8), INTENT(OUT)         :: j        ! Cartesian Y coordinate

    ! Locals 
    REAL(r8)                      :: arg
    REAL(r8)                      :: deltalon
    REAL(r8)                      :: tl1r
    REAL(r8)                      :: rm
    REAL(r8)                      :: ctl1r
    

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
    IF (proj%hemi .EQ. -1.0_r8) THEN
      i = 2.0_r8 - i  
      j = 2.0_r8 - j
    ENDIF
    RETURN
  END SUBROUTINE llij_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE set_merc(proj)
  
    ! Sets up the remaining basic elements for the mercator projection

    IMPLICIT NONE
    TYPE(proj_info), INTENT(INOUT)       :: proj
    REAL(r8)                             :: clain


    !  Preliminary variables

    clain = COS(rad_per_deg*proj%truelat1)
    proj%dlon = proj%dx / (earth_radius_m * clain)

    ! Compute distance from equator to origin, and store in the 
    ! proj%rsw tag.

    proj%rsw = 0.0_r8
    IF (proj%lat1 .NE. 0.0_r8) THEN
      proj%rsw = (ALOG(TAN(0.5_r8*((proj%lat1+90.0_r8)*rad_per_deg))))/proj%dlon
    ENDIF
    RETURN
  END SUBROUTINE set_merc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE llij_merc(lat, lon, proj, i, j)

    ! Compute i/j coordinate from lat lon for mercator projection
  
    IMPLICIT NONE
    REAL(r8), INTENT(IN)          :: lat
    REAL(r8), INTENT(IN)          :: lon
    TYPE(proj_info),INTENT(IN)    :: proj
    REAL(r8),INTENT(OUT)          :: i
    REAL(r8),INTENT(OUT)          :: j
    REAL(r8)                      :: deltalon

    deltalon = lon - proj%lon1
    IF (deltalon .LT. -180.0_r8) deltalon = deltalon + 360.0_r8
    IF (deltalon .GT. 180.0_r8) deltalon = deltalon - 360.0_r8
    i = 1.0_r8 + (deltalon/(proj%dlon*deg_per_rad))
    j = 1.0_r8 + (ALOG(TAN(0.5_r8*((lat + 90.0_r8) * rad_per_deg)))) / &
           proj%dlon - proj%rsw

    RETURN
  END SUBROUTINE llij_merc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ijll_merc(i, j, proj, lat, lon)

    ! Compute the lat/lon from i/j for mercator projection

    IMPLICIT NONE
    REAL(r8),INTENT(IN)           :: i
    REAL(r8),INTENT(IN)           :: j    
    TYPE(proj_info),INTENT(IN)    :: proj
    REAL(r8), INTENT(OUT)         :: lat
    REAL(r8), INTENT(OUT)         :: lon 


    lat = 2.0_r8*ATAN(EXP(proj%dlon*(proj%rsw + j-1.0_r8)))*deg_per_rad - 90.0_r8
    lon = (i-1.0_r8)*proj%dlon*deg_per_rad + proj%lon1
    IF (lon.GT.180.0_r8) lon = lon - 360.0_r8
    IF (lon.LT.-180.0_r8) lon = lon + 360.0_r8
    RETURN
  END SUBROUTINE ijll_merc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE llij_latlon(lat, lon, proj, i, j)
 
    ! Compute the i/j location of a lat/lon on a LATLON grid.
    IMPLICIT NONE
    REAL(r8), INTENT(IN)         :: lat
    REAL(r8), INTENT(IN)         :: lon
    TYPE(proj_info), INTENT(IN)  :: proj
    REAL(r8), INTENT(OUT)        :: i
    REAL(r8), INTENT(OUT)        :: j

    REAL(r8)                     :: deltalat
    REAL(r8)                     :: deltalon
    REAL(r8)                     :: lon360
    REAL(r8)                     :: latinc
    REAL(r8)                     :: loninc

    ! Extract the latitude and longitude increments for this grid
    ! (e.g., 2.5 deg for NCEP reanalysis) from the proj structure, where
    ! loninc is saved in the stdlon tag and latinc is saved in truelat1

    latinc = proj%truelat1
    loninc = proj%stdlon

    ! Compute deltalat and deltalon as the difference between the input 
    ! lat/lon and the origin lat/lon

    deltalat = lat - proj%lat1

    ! To account for issues around the dateline, convert the incoming
    ! longitudes to be 0->360.0_r8
    IF (lon .LT. 0) THEN 
      lon360 = lon + 360.0_r8
    ELSE 
      lon360 = lon
    ENDIF    
    deltalon = lon360 - proj%lon1      
    
    ! Compute i/j
    i = deltalon/loninc + 1.0_r8
    j = deltalat/latinc + 1.0_r8
    RETURN
    END SUBROUTINE llij_latlon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ijll_latlon(i, j, proj, lat, lon)
 
    ! Compute the lat/lon location of an i/j on a LATLON grid.
    IMPLICIT NONE
    REAL(r8), INTENT(IN)         :: i
    REAL(r8), INTENT(IN)         :: j
    TYPE(proj_info), INTENT(IN)  :: proj
    REAL(r8), INTENT(OUT)        :: lat
    REAL(r8), INTENT(OUT)        :: lon

    REAL(r8)                     :: deltalat
    REAL(r8)                     :: deltalon
    REAL(r8)                     :: latinc
    REAL(r8)                     :: loninc

    ! Extract the latitude and longitude increments for this grid
    ! (e.g., 2.5 deg for NCEP reanalysis) from the proj structure, where
    ! loninc is saved in the stdlon tag and latinc is saved in truelat1

    latinc = proj%truelat1
    loninc = proj%stdlon

    ! Compute deltalat and deltalon 

    deltalat = (j-1.0_r8)*latinc
    deltalon = (i-1.0_r8)*loninc
    lat = proj%lat1 + deltalat
    lon = proj%lon1 + deltalon

    IF ((ABS(lat) .GT. 90.0_r8).OR.(ABS(deltalon) .GT.360.0_r8)) THEN
      ! Off the earth for this grid
      lat = -999.0_r8
      lon = -999.0_r8
    ELSE
      lon = lon + 360.0_r8
      lon = MOD(lon,360.0_r8)
      IF (lon .GT. 180.0_r8) lon = lon -360.0_r8
    ENDIF

    RETURN
  END SUBROUTINE ijll_latlon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    ELSEIF ((proj%code .EQ. PROJ_MERC).OR.(proj%code .EQ. PROJ_LATLON))THEN
      utrue = ugrid
      vtrue = vgrid
    ELSE
      PRINT '(A)', 'Unrecognized projection.'
      STOP 'GRIDWIND_TO_TRUEWIND'
    ENDIF

    RETURN
  END SUBROUTINE gridwind_to_truewind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
END MODULE map_utils
