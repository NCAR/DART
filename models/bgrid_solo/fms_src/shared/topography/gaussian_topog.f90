!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gaussian_topog_mod

! <CONTACT EMAIL="bw@gfdl.noaa.gov">
!   Bruce Wyman
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   Routines for creating Gaussian-shaped land surface topography
!   for latitude-longitude grids.
! </OVERVIEW>

! <DESCRIPTION>
!   Interfaces generate simple Gaussian-shaped mountains from
!   parameters specified by either argument list or namelist input.
!   The mountain shapes are controlled by the height, half-width,
!   and ridge-width parameters.
! </DESCRIPTION>

use types_mod, only : r8
use  fms_mod, only: file_exist, open_namelist_file,  &
                    close_file,     &
                    stdlog, write_version_number,    &
                    error_mesg, FATAL

use constants_mod, only: pi

use utilities_mod, only : find_namelist_in_file, check_namelist_read

implicit none
private

public :: gaussian_topog_init, get_gaussian_topog

!-----------------------------------------------------------------------
! <NAMELIST NAME="gaussian_topog_nml">
!   <DATA NAME="height" UNITS="meter" TYPE="real(r8)" DIM="(mxmtns)" DEFAULT="0.">
!     Height in meters of the Gaussian mountains.
!    </DATA>
!   <DATA NAME="olon, olat" UNITS="degree" TYPE="real(r8)" DIM="(mxmtns)" DEFAULT="0.">
!     The longitude and latitude of mountain origins (in degrees).
!    </DATA>
!   <DATA NAME="wlon, wlat" UNITS="degree" TYPE="real(r8)" DIM="(mxmtns)" DEFAULT="0.">
!     The longitude and latitude half-width of mountain tails (in degrees).
!    </DATA>
!   <DATA NAME="rlon, rlat" UNITS="degree" TYPE="real(r8)" DIM="(mxmtns)" DEFAULT="0.">
!     The longitude and latitude half-width of mountain ridges (in degrees).  For a
!     "standard" Gaussian mountain set rlon=rlat=0.
!    </DATA>
!
!    <DATA NAME="NOTE">
!     The variables in this namelist are only used when routine
!     <TT>gaussian_topog_init</TT> is called.  The namelist variables
!     are dimensioned (by 10), so that multiple mountains can be generated.
!
!     Internal parameter mxmtns = 10. By default no mountains are generated. 
!    </DATA>

   integer, parameter :: maxmts = 10

   real(r8), dimension(maxmts) :: height = 0.
   real(r8), dimension(maxmts) ::  olon  = 0.
   real(r8), dimension(maxmts) ::  olat  = 0.
   real(r8), dimension(maxmts) ::  wlon  = 0.
   real(r8), dimension(maxmts) ::  wlat  = 0.
   real(r8), dimension(maxmts) ::  rlon  = 0.
   real(r8), dimension(maxmts) ::  rlat  = 0.

   namelist /gaussian_topog_nml/ height, olon, olat, wlon, wlat, rlon, rlat
! </NAMELIST>

!-----------------------------------------------------------------------

character(len=128) :: version = '$Revision$'
character(len=128) :: tagname = '$Id$'

logical :: do_nml = .true.
logical :: module_is_initialized = .FALSE.

!-----------------------------------------------------------------------

contains

!#######################################################################

! <SUBROUTINE NAME="gaussian_topog_init">

!   <OVERVIEW>
!     Returns a surface height field that consists
!     of the sum of one or more Gaussian-shaped mountains.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Returns a land surface topography that consists of a "set" of
!     simple Gaussian-shaped mountains.  The height, position,
!     width, and elongation of the mountains can be controlled
!     by variables in namelist <LINK SRC="#NAMELIST">&#38;gaussian_topog_nml</LINK>.
!   </DESCRIPTION>
!   <TEMPLATE>
!     <B>call gaussian_topog_init</B> ( lon, lat, zsurf )
!   </TEMPLATE>

!   <IN NAME="lon" UNITS="radians" TYPE="real(r8)" DIM="(:)">
!     The mean grid box longitude in radians.
!   </IN>
!   <IN NAME="lat" UNITS="radians" TYPE="real(r8)" DIM="(:)">
!     The mean grid box latitude in radians.
!   </IN>
!   <OUT NAME="zsurf" UNITS="meter" TYPE="real(r8)" DIM="(:,:)">
!     The surface height (in meters).
!     The size of this field must be size(lon) by size(lat).
!   </OUT>

subroutine gaussian_topog_init ( lon, lat, zsurf )

real(r8), intent(in)  :: lon(:), lat(:)
real(r8), intent(out) :: zsurf(:,:)

integer :: n

  if (.not.module_is_initialized) then
     call write_version_number( version, tagname )
  endif

  if(any(shape(zsurf) /= (/size(lon),size(lat)/))) then
    call error_mesg ('get_gaussian_topog in topography_mod', &
     'shape(zsurf) is not equal to (/size(lon),size(lat)/)', FATAL)
  endif

  if (do_nml) call read_namelist

! compute sum of all non-zero mountains
  zsurf(:,:) = 0.
  do n = 1, maxmts
    if ( height(n) == 0. ) cycle
    zsurf = zsurf + get_gaussian_topog ( lon, lat, height(n), &
                olon(n), olat(n), wlon(n), wlat(n), rlon(n), rlat(n))
  enddo
 module_is_initialized = .TRUE.                    

end subroutine gaussian_topog_init
! </SUBROUTINE>

!#######################################################################

! <FUNCTION NAME="get_gaussian_topog">

!   <OVERVIEW>
!     Returns a simple surface height field that consists of a single
!     Gaussian-shaped mountain.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Returns a single Gaussian-shaped mountain.
!     The height, position, width, and elongation of the mountain
!     is controlled by optional arguments.
!   </DESCRIPTION>
!   <TEMPLATE>
!     zsurf = <B>get_gaussian_topog</B> ( lon, lat, height
!                    [, olond, olatd, wlond, wlatd, rlond, rlatd ] )
!   </TEMPLATE>

!   <IN NAME="lon" UNITS="radians" TYPE="real(r8)" DIM="(:)">
!     The mean grid box longitude in radians.
!   </IN>
!   <IN NAME="lat" UNITS="radians" TYPE="real(r8)" DIM="(:)">
!     The mean grid box latitude in radians.
!   </IN>
!   <IN NAME="height" UNITS="meter" TYPE="real(r8)" DIM="(scalar)">
!     Maximum surface height in meters.
!   </IN>
!   <IN NAME="olond, olatd" UNITS="degrees" TYPE="real(r8)" DIM="(scalar)">
!     Position/origin of mountain in degrees longitude and latitude.
!     This is the location of the maximum height.
!   </IN>
!   <IN NAME="wlond, wlatd" UNITS="degrees" TYPE="real(r8)" DIM="(scalar)" DEFAULT="15.">
!     Gaussian half-width of mountain in degrees longitude and latitude.
!   </IN>
!   <IN NAME="rlond, rlatd" UNITS="degrees" TYPE="real(r8)" DIM="(scalar)" DEFAULT="0.">
!     Ridge half-width of mountain in degrees longitude and latitude.
!                    This is the elongation of the maximum height.
!   </IN>
!   <OUT NAME="zsurf" UNITS="meter" TYPE="real(r8)" DIM="(:,:)">
!     The surface height (in meters).
!              The size of the returned field is size(lon) by size(lat).
!   </OUT>
!   <ERROR MSG="shape(zsurf) is not equal to (/size(lon),size(lat)/)" STATUS="FATAL">
!     Check the input grid size and output field size.
!     The input grid is defined at the midpoint of grid boxes.
!   </ERROR>
!   <NOTE>
!     Mountains do not wrap around the poles.
!   </NOTE>

function get_gaussian_topog ( lon, lat, height,                          &
                              olond, olatd, wlond, wlatd, rlond, rlatd ) &
                     result ( zsurf )

real(r8), intent(in)  :: lon(:), lat(:)
real(r8), intent(in)  :: height
real(r8), intent(in), optional :: olond, olatd, wlond, wlatd, rlond, rlatd
real(r8) :: zsurf(size(lon),size(lat))

integer :: i, j
real(r8)    :: olon, olat, wlon, wlat, rlon, rlat
real(r8)    :: tpi, dtr, dx, dy, xx, yy

  if (do_nml) call read_namelist

! no need to compute mountain if height=0
  if ( height == 0. ) then
       zsurf(:,:) = 0.
       return
  endif

  tpi = 2.0*pi
  dtr = tpi/360.

! defaults and convert degrees to radians (dtr)
  olon = 90.*dtr;  if (present(olond)) olon=olond*dtr
  olat = 45.*dtr;  if (present(olatd)) olat=olatd*dtr
  wlon = 15.*dtr;  if (present(wlond)) wlon=wlond*dtr
  wlat = 15.*dtr;  if (present(wlatd)) wlat=wlatd*dtr
  rlon =  0.    ;  if (present(rlond)) rlon=rlond*dtr
  rlat =  0.    ;  if (present(rlatd)) rlat=rlatd*dtr

! compute gaussian-shaped mountain
    do j=1,size(lat)
      dy = abs(lat(j) - olat)   ! dist from y origin
      yy = max(0.0_r8, dy-rlat)/wlat
      do i=1,size(lon)
        dx = abs(lon(i) - olon) ! dist from x origin
        dx = min(dx, abs(dx-tpi))  ! To ensure that: -pi <= dx <= pi
        xx = max(0.0_r8, dx-rlon)/wlon
        zsurf(i,j) = height*exp(-xx**2 - yy**2)
      enddo
    enddo

end function get_gaussian_topog
! </FUNCTION>

!#######################################################################

subroutine read_namelist

   integer :: unit, ierr, io, iunit
   real(r8)    :: dtr

!  read namelist

! Replace with new dart code, 8 June, 2006
call find_namelist_in_file("input.nml", "gaussian_topog_nml", iunit)
read(iunit, nml = gaussian_topog_nml, iostat = io)
call check_namelist_read(iunit, io, "gaussian_topog__nml")

!!!   if ( file_exist('input.nml')) then
!!!      unit = open_namelist_file ( )
!!!      ierr=1; do while (ierr /= 0)
!!!         read  (unit, nml=gaussian_topog_nml, iostat=io, end=10)
!!!         ierr = check_nml_error(io,'gaussian_topog_nml')
!!!      enddo
!!! 10   call close_file (unit)
!!!   endif

!  write version and namelist to log file

   write (stdlog(), nml=gaussian_topog_nml)

   do_nml = .false.

end subroutine read_namelist

!#######################################################################

end module gaussian_topog_mod

! <INFO>
!   <NOTE>
!     NAMELIST FOR GENERATING GAUSSIAN MOUNTAINS
!
!  * multiple mountains can be generated
!  * the final mountains are the sum of all
!
!       height = height in meters
!       olon, olat = longitude,latitude origin              (degrees)
!       rlon, rlat = longitude,latitude half-width of ridge (degrees)
!       wlon, wlat = longitude,latitude half-width of tail  (degrees)
!
!       Note: For the standard gaussian mountain
!             set rlon = rlat = 0 .
!
! <PRE>
!
!       height -->   ___________________________
!                   /                           \
!                  /              |              \
!    gaussian     /               |               \
!      sides --> /                |                \
!               /               olon                \
!         _____/                olat                 \______
!
!              |    |             |
!              |<-->|<----------->|
!              |wlon|    rlon     |
!               wlat     rlat
!
! </PRE>
!
!See the <LINK SRC="topography.html#TEST PROGRAM">topography </LINK>module documentation for a test program.
!   </NOTE>
! </INFO>
