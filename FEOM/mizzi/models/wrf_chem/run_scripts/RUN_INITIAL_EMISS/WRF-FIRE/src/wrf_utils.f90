
   module wrf_utils

   implicit none

!  private
   public :: llij
   public :: ijll
   public :: wrf_file
   public :: proj_info
   public :: proj_latlon, pi, rad_per_deg

   integer, parameter  :: PROJ_LATLON  = 0
   integer, parameter  :: PROJ_LC      = 1
   integer, parameter  :: PROJ_PS      = 2
   integer, parameter  :: PROJ_MERC    = 3
   integer, parameter  :: PROJ_CASSINI = 6
   real, parameter :: EARTH_RADIUS_M = 6370000.
   real, parameter :: PI = 3.141592653589793
   real, parameter :: DEG_PER_RAD = 180./PI
   real, parameter :: RAD_PER_DEG = PI/180.
   character(len=19)   :: proj_name(0:3) = (/ 'LATLON             ', 'LAMBERT            ', &
                                              'POLAR STEREOGRAPHIC', 'MERCATOR           ' /)


   TYPE proj_info
      integer          :: code     ! Integer code for projection TYPE
      integer          :: ixdim    ! For Rotated Lat/Lon -- number of mass points in an odd row
      integer          :: jydim    ! For Rotated Lat/Lon -- number of rows
      integer          :: ide
      integer          :: jde
      real             :: lat1     ! SW latitude (1,1) in degrees (-90->90N)
      real             :: lon1     ! SW longitude (1,1) in degrees (-180->180E)
      real             :: lat0     ! For Cassini, latitude of projection pole
      real             :: lon0     ! For Cassini, longitude of projection pole
      real             :: dx = 0.  ! Grid spacing in meters at truelats, used only for ps, lc, and merc projections
      real             :: dy       ! Grid spacing in meters at truelats, used only for ps, lc, and merc projections
      real             :: latinc   ! Latitude increment for cylindrical lat/lon
      real             :: loninc   ! Longitude increment for cylindrical lat/lon also the lon increment for Gaussian grid
      real             :: dlat     ! Lat increment for lat/lon grids
      real             :: dlon     ! Lon increment for lat/lon grids
      real             :: stdlon   ! Longitude parallel to y-axis (-180->180E)
      real             :: truelat1 ! First true latitude (all projections)
      real             :: truelat2 ! Second true lat (LC only)
      real             :: hemi     ! 1 for NH, -1 for SH
      real             :: cone     ! Cone factor for LC projections
      real             :: polei    ! Computed i-location of pole point
      real             :: polej    ! Computed j-location of pole point
      real             :: rsw      ! Computed radius to SW corner
      real             :: rebydx   ! Earth radius divided by dx
      real             :: knowni   ! X-location of known lat/lon
      real             :: knownj   ! Y-location of known lat/lon
      real             :: re_m     ! Radius of spherical earth, meters
      real             :: rho0     ! For Albers equal area
      real             :: nc       ! For Albers equal area
      real             :: bigc     ! For Albers equal area
      LOGICAL          :: init     ! Flag to indicate if this struct is ready for use
      LOGICAL          :: wrap     ! For Gaussian -- flag to indicate wrapping around globe?
      LOGICAL          :: comp_ll  ! Work in computational lat/lon space for Cassini
   end TYPE proj_info


!------------------------------------------------------------------
!     include files
!------------------------------------------------------------------
   include 'netcdf.inc'

   CONTAINS

   subroutine proj_init( proj )
!-------------------------------------------------------------
!  ... intialize wrf map projection
!-------------------------------------------------------------
   type(proj_info) :: proj

!-------------------------------------------------------------
!  ... local variables
!-------------------------------------------------------------
   integer :: astat
   integer :: grid

   if( proj%code < 1 .or. proj%code > 3 ) then
     write(*,'('' proj_init: input projection, '',i2,'', is out of bounds [1-3]'')') proj%code
     stop
   else
     write(*,*) ' '
     write(*,'('' proj_init: projection = '',i2)') proj%code
   endif

   proj%ixdim    = proj%ide+1
   proj%jydim    = proj%jde+1

   proj%knowni   = proj%ixdim/2.
   proj%knownj   = proj%jydim/2.
   proj%init     = .true.
   proj%re_m     = EARTH_RADIUS_M

   if (proj%truelat1 < 0.) then
     proj%hemi = -1.0 
   else
     proj%hemi = 1.0
   endif
   proj%rebydx = proj%re_m / proj%dx

   if( proj%code == PROJ_LC ) then
     if( abs(proj%truelat2) > 90. ) then
       proj%truelat2 = proj%truelat1
     end if
   endif

   call setup_projection( proj )

!-------------------------------------------------------------
!  ... a few projection variable diagnostics
!-------------------------------------------------------------
   write(*,*) 'proj_init: proj%hemi    = ',proj%hemi
   write(*,*) 'proj_init: proj%rebydx  = ',proj%rebydx
   write(*,*) 'proj_init: proj%polei,j = ',proj%polei,proj%polej
   write(*,'('' proj_init: west-east,south-north = '',2i5)') proj%ixdim-1,proj%jydim-1

   end subroutine proj_init

   subroutine setup_projection( proj )
!-----------------------------------------------------------
! setup the mapping
!-----------------------------------------------------------
!-----------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------
   type(proj_info) :: proj

   if( proj%code == PROJ_LC ) then
     call set_lc( proj )
   elseif( proj%code == PROJ_PS ) then
     call set_ps( proj )
   elseif( proj%code == PROJ_MERC ) then
     call set_merc( proj )
   endif

   end subroutine setup_projection

   subroutine llij( lat, lon, proj, i, j)
!-----------------------------------------------------------
! Subroutine to convert geographical latitude,longitude values to cartesian x/y
!-----------------------------------------------------------

!-----------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------
   real, intent(in)              :: lat      ! Latitude (-90->90 deg N)
   real, intent(in)              :: lon      ! Longitude (-180->180 E)
   real, intent(out)             :: i        ! Cartesian X coordinate
   real, intent(out)             :: j        ! Cartesian Y coordinate
   type(proj_info)   :: proj

   if( proj%code == PROJ_LC ) then
     call llij_lc( proj, lat, lon, i, j)
   elseif( proj%code == PROJ_PS ) then
     call llij_ps( proj, lat, lon, i, j)
   elseif( proj%code == PROJ_MERC ) then
     call llij_merc( proj, lat, lon, i, j)
   elseif( proj%code == PROJ_LATLON ) then
     call llij_latlon( proj, lat, lon, i, j)
   endif

   end subroutine llij

   subroutine ijll( i, j, proj, lat, lon )
!-----------------------------------------------------------
! Subroutine to convert cartesian x/y values to longitude,latitude values
!-----------------------------------------------------------

!-----------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------
   real, intent(in)              :: i        ! Cartesian X coordinate
   real, intent(in)              :: j        ! Cartesian Y coordinate
   real, intent(out)             :: lat      ! Latitude (-90->90 deg N)
   real, intent(out)             :: lon      ! Longitude (-180->180 E)
   type(proj_info)   :: proj

   if( proj%code == PROJ_LC ) then
     call ijll_lc( proj, i, j, lat, lon )
   elseif( proj%code == PROJ_PS ) then
     call ijll_ps( proj, i, j, lat, lon )
   elseif( proj%code == PROJ_MERC ) then
     call ijll_merc( proj, i, j, lat, lon )
   elseif( proj%code == PROJ_LATLON ) then
     call ijll_latlon( proj, i, j, lat, lon )
   endif

   end subroutine ijll

   subroutine set_lc( proj )
!-----------------------------------------------------------
! Initialize the remaining items in the proj structure for a
! lambert conformal grid.
!-----------------------------------------------------------

      real                               :: arg
      real                               :: deltalon1
      real                               :: tl1r
      real                               :: ctl1r
      type(proj_info)   :: proj
  
      ! Compute cone factor
      CALL lc_cone( proj%truelat1, proj%truelat2, proj%cone )
  
      ! Compute longitude differences and ensure we stay out of the
      ! forbidden "cut zone"
      deltalon1 = proj%lon1 - proj%stdlon
      IF (deltalon1 > 180.) then
        deltalon1 = deltalon1 - 360.
      elseIF (deltalon1 < -180.) then
        deltalon1 = deltalon1 + 360.
      endif
  
      ! Convert truelat1 to radian and compute COS for later use
      tl1r = proj%truelat1 * rad_per_deg
      ctl1r = COS(tl1r)
  
      ! Compute the radius to our known lower-left (SW) corner
      proj%rsw = proj%rebydx * ctl1r/proj%cone * &
             (TAN((90.*proj%hemi-proj%lat1)*rad_per_deg/2.) / &
              TAN((90.*proj%hemi-proj%truelat1)*rad_per_deg/2.))**proj%cone
  
      ! Find pole point
      arg = proj%cone*(deltalon1*rad_per_deg)
      proj%polei = proj%hemi*proj%knowni - proj%hemi * proj%rsw * Sin(arg)
      proj%polej = proj%hemi*proj%knownj + proj%rsw * COS(arg)  
  
   end subroutine set_lc                             

   subroutine lc_cone( truelat1, truelat2, cone )
!-----------------------------------------------------------
! compute the cone factor of a Lambert Conformal projection
!-----------------------------------------------------------

      ! Input Args
      real, intent(in)             :: truelat1  ! (-90 -> 90 degrees N)
      real, intent(in)             :: truelat2  !   "   "  "   "     "
  
      ! Output Args
      real, intent(out)            :: cone
  
      ! Locals
  
      ! BEGin CODE
  
      ! First, see if this is a secant or tangent projection.  For tangent
      ! projections, truelat1 = truelat2 and the cone is tangent to the 
      ! Earth's surface at this latitude.  For secant projections, the cone
      ! intersects the Earth's surface at each of the distinctly different
      ! latitudes
      IF (ABS(truelat1-truelat2) > 0.1) THEN
         cone = ALOG10(COS(truelat1*rad_per_deg)) - &
                ALOG10(COS(truelat2*rad_per_deg))
         cone = cone /(ALOG10(TAN((45.0 - ABS(truelat1)/2.0) * rad_per_deg)) - &
                ALOG10(TAN((45.0 - ABS(truelat2)/2.0) * rad_per_deg)))        
      ELSE
         cone = Sin(ABS(truelat1)*rad_per_deg )  
      endIF

      RETURN

   end subroutine lc_cone

   subroutine llij_lc( proj, lat, lon, i, j)
!-----------------------------------------------------------
! Subroutine to compute the geographical latitude and longitude values
! to the cartesian x/y on a Lambert Conformal projection.
!-----------------------------------------------------------
     
      ! Input Args
      real, intent(in)              :: lat      ! Latitude (-90->90 deg N)
      real, intent(in)              :: lon      ! Longitude (-180->180 E)
      type(proj_info)    :: proj
  
      ! Output Args                 
      real, intent(out)             :: i        ! Cartesian X coordinate
      real, intent(out)             :: j        ! Cartesian Y coordinate
  
      ! Locals 
      real                          :: arg
      real                          :: deltalon
      real                          :: tl1r
      real                          :: rm
      real                          :: ctl1r
      
  
      ! Compute deltalon between known longitude and standard lon and ensure
      ! it is not in the cut zone
      deltalon = lon - proj%stdlon
      IF (deltalon .GT. +180.) deltalon = deltalon - 360.
      IF (deltalon .LT. -180.) deltalon = deltalon + 360.
      
      ! Convert truelat1 to radian and compute COS for later use
      tl1r = proj%truelat1 * rad_per_deg
      ctl1r = COS(tl1r)     
     
      ! Radius to desired point
      rm = proj%rebydx * ctl1r/proj%cone * &
           (TAN((90.*proj%hemi-lat)*rad_per_deg/2.) / &
            TAN((90.*proj%hemi-proj%truelat1)*rad_per_deg/2.))**proj%cone
  
      arg = proj%cone*(deltalon*rad_per_deg)
      i = proj%polei + proj%hemi * rm * Sin(arg)
      j = proj%polej - rm * COS(arg)
  
      ! Finally, if we are in the southern hemisphere, flip the i/j
      ! values to a coordinate system where (1,1) is the SW corner
      ! (what we assume) which is different than the original NCEP
      ! algorithms which used the NE corner as the origin in the 
      ! southern hemisphere (left-hand vs. right-hand coordinate?)
      i = proj%hemi * i  
      j = proj%hemi * j

   end subroutine llij_lc

   subroutine set_merc( proj )
!-----------------------------------------------------------
! Sets up the remaining basic elements for the mercator projection
!-----------------------------------------------------------
      type(proj_info)    :: proj
  
      real                                 :: clain
  
  
      !  Preliminary variables
  
      clain = COS(rad_per_deg*proj%truelat1)
      proj%dlon = proj%dx / (proj%re_m * clain)
  
      ! Compute distance from equator to origin, and store in the 
      ! proj%rsw tag.
  
      proj%rsw = 0.
      IF (proj%lat1 .NE. 0.) THEN
         proj%rsw = (ALOG(TAN(0.5*((proj%lat1+90.)*rad_per_deg))))/proj%dlon
      endIF

   end subroutine set_merc

   subroutine llij_merc( proj, lat, lon, i, j )
!-----------------------------------------------------------
! Compute i/j coordinate from lat lon for mercator projection
!-----------------------------------------------------------
    
!-----------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------
      real, intent(in)              :: lat
      real, intent(in)              :: lon
      real,intent(out)              :: i
      real,intent(out)              :: j
      type(proj_info)   :: proj

!-----------------------------------------------------------
! local variables
!-----------------------------------------------------------
      real                          :: deltalon
  
      deltalon = lon - proj%lon1
      IF (deltalon > -180.) deltalon = deltalon + 360.
      IF (deltalon > 180.) deltalon = deltalon - 360.
      i = proj%knowni + (deltalon/(proj%dlon*deg_per_rad))
      j = proj%knownj + (ALOG(TAN(0.5*((lat + 90.) * rad_per_deg)))) / &
             proj%dlon - proj%rsw
  
   end subroutine llij_merc

   subroutine set_ps( proj )
!-----------------------------------------------------------
! Initializes a polar-stereographic map projection from the partially
! filled proj structure. This routine computes the radius to the
! southwest corner and computes the i/j location of the pole for use
! in llij_ps and ijll_ps.
!-----------------------------------------------------------
      type(proj_info)    :: proj

!-----------------------------------------------------------
! Local vars
!-----------------------------------------------------------
      real                              :: ala1
      real                              :: alo1
      real                              :: reflon
      real                              :: scale_top
  
      ! Executable code
      reflon = proj%stdlon + 90.
  
      ! Compute numerator term of map scale factor
      scale_top = 1. + proj%hemi * Sin(proj%truelat1 * rad_per_deg)
  
      ! Compute radius to lower-left (SW) corner
      ala1 = proj%lat1 * rad_per_deg
      proj%rsw = proj%rebydx*COS(ala1)*scale_top/(1.+proj%hemi*Sin(ala1))
  
      ! Find the pole point
      alo1 = (proj%lon1 - reflon) * rad_per_deg
      proj%polei = proj%knowni - proj%rsw * COS(alo1)
      proj%polej = proj%knownj - proj%hemi * proj%rsw * Sin(alo1)

   end subroutine set_ps

   subroutine llij_ps( proj, lat, lon, i, j )
!-----------------------------------------------------------
! Given latitude (-90 to 90), longitude (-180 to 180), and the
! standard polar-stereographic projection information via the 
! public proj structure, this routine returns the i/j indices which
! if within the domain range from 1->nx and 1->ny, respectively.
!-----------------------------------------------------------
  
      ! Delcare input arguments
      real, intent(in)               :: lat
      real, intent(in)               :: lon
      type(proj_info)    :: proj
  
      ! Declare output arguments     
      real, intent(out)              :: i !(x-index)
      real, intent(out)              :: j !(y-index)
  
      ! Declare local variables
      
      real                           :: reflon
      real                           :: scale_top
      real                           :: ala
      real                           :: alo
      real                           :: rm
  
      ! BEGin CODE
    
      reflon = proj%stdlon + 90.
     
      ! Compute numerator term of map scale factor
  
      scale_top = 1. + proj%hemi * Sin(proj%truelat1 * rad_per_deg)
  
      ! Find radius to desired point
      ala = lat * rad_per_deg
      rm = proj%rebydx * COS(ala) * scale_top/(1. + proj%hemi *Sin(ala))
      alo = (lon - reflon) * rad_per_deg
      i = proj%polei + rm * COS(alo)
      j = proj%polej + proj%hemi * rm * Sin(alo)
   
   end subroutine llij_ps

   subroutine ijll_lc( proj, i, j, lat, lon)
 
   ! Subroutine to convert from the (i,j) cartesian coordinate to the 
   ! geographical latitude and longitude for a Lambert Conformal projection.
 
   ! History:
   ! 25 Jul 01: Corrected by B. Shaw, NOAA/FSL
   ! 
      IMPLICIT NONE
  
      ! Input Args
      REAL, INTENT(IN)              :: i        ! Cartesian X coordinate
      REAL, INTENT(IN)              :: j        ! Cartesian Y coordinate
      type(proj_info), intent(in)   :: proj
  
      ! Output Args                 
      REAL, INTENT(OUT)             :: lat      ! Latitude (-90->90 deg N)
      REAL, INTENT(OUT)             :: lon      ! Longitude (-180->180 E)
  
      ! Locals 
      REAL                          :: inew
      REAL                          :: jnew
      REAL                          :: r
      REAL                          :: chi,chi1,chi2
      REAL                          :: r2
      REAL                          :: xx
      REAL                          :: yy
  
      ! BEGIN CODE
  
      chi1 = (90. - proj%hemi*proj%truelat1)*rad_per_deg
      chi2 = (90. - proj%hemi*proj%truelat2)*rad_per_deg
  
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
      IF (r2 == 0.) THEN
         lat = proj%hemi * 90.
         lon = proj%stdlon
      ELSE
         
         ! Longitude
         lon = proj%stdlon + deg_per_rad * ATAN2(proj%hemi*xx,yy)/proj%cone
         lon = AMOD(lon+360., 360.)
   
         ! Latitude.  Latitude determined by solving an equation adapted 
         ! from:
         !  Maling, D.H., 1973: Coordinate Systems and Map Projections
         ! Equations #20 in Appendix I.  
           
         IF (chi1 .EQ. chi2) THEN
            chi = 2.0*ATAN( ( r/TAN(chi1) )**(1./proj%cone) * TAN(chi1*0.5) )
         ELSE
            chi = 2.0*ATAN( (r*proj%cone/SIN(chi1))**(1./proj%cone) * TAN(chi1*0.5)) 
         ENDIF
         lat = (90.0-chi*deg_per_rad)*proj%hemi
  
      ENDIF
  
      IF (lon .GT. +180.) lon = lon - 360.
      IF (lon .LT. -180.) lon = lon + 360.
 
   END SUBROUTINE ijll_lc

   SUBROUTINE ijll_ps( proj, i, j, lat, lon )
 
      ! This is the inverse subroutine of llij_ps.  It returns the 
      ! latitude and longitude of an i/j point given the projection info 
      ! structure.  
  
      IMPLICIT NONE
  
      ! Declare input arguments
      REAL, INTENT(IN)                    :: i    ! Column
      REAL, INTENT(IN)                    :: j    ! Row
      type(proj_info)         :: proj
      
      ! Declare output arguments
      REAL, INTENT(OUT)                   :: lat     ! -90 -> 90 north
      REAL, INTENT(OUT)                   :: lon     ! -180 -> 180 East
  
      ! Local variables
      REAL                                :: reflon
      REAL                                :: scale_top
      REAL                                :: xx,yy
      REAL                                :: gi2, r2
      REAL                                :: arccos
  
      ! Begin Code
  
      ! Compute the reference longitude by rotating 90 degrees to the east
      ! to find the longitude line parallel to the positive x-axis.
      reflon = proj%stdlon + 90.
     
      ! Compute numerator term of map scale factor
      scale_top = 1. + proj%hemi * SIN(proj%truelat1 * rad_per_deg)
  
      ! Compute radius to point of interest
      xx = i - proj%polei
      yy = (j - proj%polej) * proj%hemi
      r2 = xx**2 + yy**2
  
      ! Now the magic code
      IF (r2 .EQ. 0.) THEN 
         lat = proj%hemi * 90.
         lon = reflon
      ELSE
         gi2 = (proj%rebydx * scale_top)**2.
         lat = deg_per_rad * proj%hemi * ASIN((gi2-r2)/(gi2+r2))
         arccos = ACOS(xx/SQRT(r2))
         IF (yy .GT. 0) THEN
            lon = reflon + deg_per_rad * arccos
         ELSE
            lon = reflon - deg_per_rad * arccos
         ENDIF
      ENDIF
    
      ! Convert to a -180 -> 180 East convention
      IF (lon > 180.) then
        lon = lon - 360.
      ELSEIF (lon < -180.) then
        lon = lon + 360.
      ENDIF

   END SUBROUTINE ijll_ps

   SUBROUTINE ijll_merc( proj, i, j, lat, lon )
 
      ! Compute the lat/lon from i/j for mercator projection
  
      IMPLICIT NONE
      REAL,INTENT(IN)               :: i
      REAL,INTENT(IN)               :: j    
      REAL, INTENT(OUT)             :: lat
      REAL, INTENT(OUT)             :: lon 
      type(proj_info)         :: proj
  
  
      lat = 2.0*ATAN(EXP(proj%dlon*(proj%rsw + j-proj%knownj)))*deg_per_rad - 90.
      lon = (i-proj%knowni)*proj%dlon*deg_per_rad + proj%lon1
      IF (lon > 180.) lon = lon - 360.
      IF (lon < -180.) lon = lon + 360.

   END SUBROUTINE ijll_merc

   SUBROUTINE llij_latlon( proj, lat, lon, i, j )
  
      ! Compute the i/j location of a lat/lon on a LATLON grid.
      IMPLICIT NONE
      REAL, INTENT(IN)             :: lat
      REAL, INTENT(IN)             :: lon
      TYPE(proj_info), INTENT(IN)  :: proj
      REAL, INTENT(OUT)            :: i
      REAL, INTENT(OUT)            :: j
  
      REAL                         :: deltalat
      REAL                         :: deltalon
  
      ! Compute deltalat and deltalon as the difference between the input 
      ! lat/lon and the origin lat/lon
      deltalat = lat - proj%lat1
      deltalon = lon - proj%lon1      
      
      ! Compute i/j
      i = deltalon/proj%loninc
      j = deltalat/proj%latinc

      i = i + proj%knowni
      j = j + proj%knownj
  
   END SUBROUTINE llij_latlon

   SUBROUTINE ijll_latlon( proj, i, j, lat, lon )
  
      ! Compute the lat/lon location of an i/j on a LATLON grid.
      IMPLICIT NONE
      REAL, INTENT(IN)             :: i
      REAL, INTENT(IN)             :: j
      TYPE(proj_info), INTENT(IN)  :: proj
      REAL, INTENT(OUT)            :: lat
      REAL, INTENT(OUT)            :: lon
  
      REAL                         :: i_work, j_work
      REAL                         :: deltalat
      REAL                         :: deltalon
  
      i_work = i - proj%knowni
      j_work = j - proj%knownj

      ! Compute deltalat and deltalon 
      deltalat = j_work*proj%latinc
      deltalon = i_work*proj%loninc
  
      lat = proj%lat1 + deltalat
      lon = proj%lon1 + deltalon
  
   END SUBROUTINE ijll_latlon

   subroutine wrf_file( domain, wrf_dir, dx, proj, ngatts, attrs )
!---------------------------------------------------------------------
!   read wrf file
!---------------------------------------------------------------------
   use attr_types, only : glb_att

!---------------------------------------------------------------------
!   dummy arguments
!---------------------------------------------------------------------
    integer, intent(in)          :: domain
    integer, intent(inout)       :: ngatts
    real, intent(out)            :: dx
    type(proj_info)              :: proj
    type(glb_att), pointer       :: attrs(:)
    character(len=*), intent(in) :: wrf_dir
!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   integer            :: ncid
   integer            :: dimid
   integer            :: varid
   integer            :: astat
   real, allocatable  :: wrk(:,:)
   character(len=3)   :: num
   character(len=128) :: message
   character(len=256) :: filespec

   write(num,'(i3)') 100+domain
   filespec = trim( wrf_dir ) // 'wrfinput_d' // num(2:3)
!---------------------------------------------------------------------
!   open wrf input file
!---------------------------------------------------------------------
   message = 'wrf_file: Failed to open ' // trim(filespec)
   call handle_ncerr( nf_open( trim(filespec), nf_noclobber, ncid ), message )       
!---------------------------------------------------------------------
!   get wrf dimesions
!---------------------------------------------------------------------
   message = 'Failed to get lon dimension id'
   call handle_ncerr( nf_inq_dimid( ncid, 'west_east', dimid ), message )
   message = 'Failed to get lon dimension'
   call handle_ncerr( nf_inq_dimlen( ncid, dimid, proj%ide ), message )
   message = 'Failed to get lat dimension id'
   call handle_ncerr( nf_inq_dimid( ncid, 'south_north', dimid ), message )
   message = 'Failed to get lat dimension'
   call handle_ncerr( nf_inq_dimlen( ncid, dimid, proj%jde ), message )
!---------------------------------------------------------------------
!   get wrf map projection variables
!---------------------------------------------------------------------
   message = 'Failed to get map_proj'
   call handle_ncerr( nf_get_att_int( ncid, nf_global, 'MAP_PROJ', proj%code ), message )
   if( proj%code == PROJ_CASSINI ) then
     proj%code = PROJ_LATLON
   endif
   write(*,*) ' '
   write(*,*) 'wrf_file: MAP_PROJ is ',trim(proj_name(proj%code))

is_not_latlon : &
   if( proj%code /= PROJ_LATLON ) then
     message = 'wrf_file: Failed to get cen_lon'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'CEN_LON', proj%lon1 ), message )
     write(*,*) 'wrf_file: CEN_LON = ',proj%lon1
     message = 'wrf_file: Failed to get cen_lat'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'CEN_LAT', proj%lat1 ), message )
     write(*,*) 'wrf_file: CEN_LAT = ',proj%lat1
     message = 'wrf_file: Failed to get stand_lon'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'STAND_LON', proj%stdlon ), message )
     write(*,*) 'wrf_file: STAND_LON = ',proj%stdlon
     message = 'wrf_file: Failed to get truelat1'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'TRUELAT1', proj%truelat1 ), message )
     write(*,*) 'wrf_file: TRUELAT1 = ',proj%truelat1
     message = 'wrf_file: Failed to get truelat2'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'TRUELAT2', proj%truelat2 ), message )
     write(*,*) 'wrf_file: TRUELAT2 = ',proj%truelat2
     message = 'wrf_file: Failed to get dx'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'DX', proj%dx ), message )
     write(*,*) 'wrf_file: DX = ',1.e-3*proj%dx,' km'
     message = 'wrf_file: Failed to get dy'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'DY', proj%dy ), message )
!---------------------------------------------------------------------
!   initialize map projection
!---------------------------------------------------------------------
     call proj_init( proj )
   else is_not_latlon
     message = 'wrf_file: Failed to get XLONG variable id'
     call handle_ncerr( nf_inq_varid( ncid, 'XLONG', varid ), message )
     allocate( wrk(proj%ide,proj%jde),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'wrf_file: failed to allocate wrk variable; error = ',astat
       stop 'Alloc error'
     endif
     message = 'wrf_file: Failed to read XLONG variable'
     call handle_ncerr( nf_get_var_real( ncid, varid, wrk ), message )
     proj%lon1   = wrk(1,1)
     proj%loninc = wrk(2,1) - wrk(1,1)
     message = 'wrf_file: Failed to get XLAT variable id'
     call handle_ncerr( nf_inq_varid( ncid, 'XLAT', varid ), message )
     message = 'wrf_file: Failed to read XLAT variable'
     call handle_ncerr( nf_get_var_real( ncid, varid, wrk ), message )
     proj%lat1   = wrk(1,1)
     proj%latinc = wrk(1,2) - wrk(1,1)
     proj%knowni = 1.
     proj%knownj = 1.
     deallocate( wrk )
     proj%dx = earth_radius_m
     proj%dy = earth_radius_m * proj%latinc * rad_per_deg
   endif is_not_latlon

   call get_fire_emis_glb_atts( ngatts, attrs, ncid )

   dx = proj%dx

!---------------------------------------------------------------------
!   close wrf file
!---------------------------------------------------------------------
   message = 'wrf_file: Failed to close ' // trim(filespec)
   call handle_ncerr( nf_close( ncid ), message )       

   end subroutine wrf_file

   subroutine get_fire_emis_glb_atts( ngatts, attrs, ncid )
!---------------------------------------------------------------------
!   read the global attributes
!---------------------------------------------------------------------
   use attr_types, only : glb_att

   implicit none

!---------------------------------------------------------------------
!   dummy arguments
!---------------------------------------------------------------------
    integer, intent(in)    :: ncid
    integer, intent(out)   :: ngatts
    type(glb_att), pointer :: attrs(:)

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   integer              :: m
   integer              :: astat
   integer              :: attr_len
   character(len=132)   :: message
   character(len=132)   :: attr_name

!---------------------------------------------------------------------
!   get global attr count
!---------------------------------------------------------------------
   message = 'glb_attr: Failed to get glb attr count'
   call handle_ncerr( nf_inq_natts( ncid, ngatts ), message )       
!---------------------------------------------------------------------
!   allocate variables
!---------------------------------------------------------------------
   if( ngatts > 0 ) then
     allocate( attrs(ngatts),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'glb_attr: failed to allocate type glb_att'
       stop 'Alloc err'
     endif
     attrs(:)%name = ' '
   endif
!---------------------------------------------------------------------
!   loop over glb attributes
!---------------------------------------------------------------------
glb_attr_loop : &
   do m = 1,ngatts
     write(message,*) 'glb_attr: Failed to get glb attr # ',m,' name'
     call handle_ncerr( nf_inq_attname( ncid, nf_global, m, attr_name ), message )       
     attrs(m)%name = attr_name
     write(message,*) 'glb_attr: Failed to get glb attr # ',m,' type,len'
     call handle_ncerr( nf_inq_att( ncid, nf_global, trim(attr_name), attrs(m)%type, attr_len ), message )       
     attrs(m)%len = attr_len
     message = 'glb_attr: Failed to get ' // trim(attr_name)
     select case( attrs(m)%type )
       case( nf_byte )
         allocate( attrs(m)%attr_byte(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_byte'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_int1( ncid, nf_global, trim(attr_name), attrs(m)%attr_byte ), message )       
       case( nf_char )
         attrs(m)%attr_char = ' '
         call handle_ncerr( nf_get_att_text( ncid, nf_global, trim(attr_name), attrs(m)%attr_char ), message )       
       case( nf_short )
         allocate( attrs(m)%attr_short(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_short'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_int2( ncid, nf_global, trim(attr_name), attrs(m)%attr_short ), message )       
       case( nf_int )
         allocate( attrs(m)%attr_int(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_int'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_int( ncid, nf_global, trim(attr_name), attrs(m)%attr_int ), message )       
       case( nf_float )
         allocate( attrs(m)%attr_real(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_real'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_real( ncid, nf_global, trim(attr_name), attrs(m)%attr_real ), message )       
       case( nf_double )
         allocate( attrs(m)%attr_dbl(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_dbl'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_double( ncid, nf_global, trim(attr_name), attrs(m)%attr_dbl ), message )       
     end select
   end do glb_attr_loop

   end subroutine get_fire_emis_glb_atts

   subroutine handle_ncerr( ret, mes )
!---------------------------------------------------------------------
!	... netcdf error handling routine
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!	... dummy arguments
!---------------------------------------------------------------------
   integer, intent(in) :: ret
   character(len=*), intent(in) :: mes

   if( ret /= nf_noerr ) then
      write(*,*) 'handle_ncerr: ',trim(mes)
      write(*,*) nf_strerror( ret )
      stop 'netcdf error'
   endif

   end subroutine handle_ncerr

   end module wrf_utils
