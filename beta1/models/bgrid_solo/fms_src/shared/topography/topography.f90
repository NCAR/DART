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

module topography_mod

! <CONTACT EMAIL="bw@gfdl.noaa.gov">
!   Bruce Wyman
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   Routines for creating land surface topography fields and land-water masks
!   for latitude-longitude grids.
! </OVERVIEW>

! <DESCRIPTION>
!   This module generates real(r8)istic mountains and land-water masks
!   on a specified latitude-longitude grid by interpolating from the
!   1/6 degree Navy mean topography and percent water data sets.
!   The fields that can be generated are mean and standard deviation
!   of topography within the specified grid boxes; and land-ocean (or
!   water) mask and land-ocean (or water) fractional area.
!
!   The interpolation scheme conserves the area-weighted average
!   of the input data by using module horiz_interp.
!
!   The interfaces get_gaussian_topog and gaussian_topog_init are documented in <LINK SRC="gaussian_topog.html">gaussian_topog_mod</LINK>.
! </DESCRIPTION>

use types_mod, only : r8
use gaussian_topog_mod, only: gaussian_topog_init, get_gaussian_topog
use   horiz_interp_mod, only: horiz_interp

use utilities_mod, only : find_namelist_in_file, check_namelist_read
use            fms_mod, only: file_exist,               &
                              open_namelist_file, close_file, stdlog,    &
                              write_version_number, &
                              open_ieee32_file, error_mesg, FATAL

implicit none
private

public :: topography_init,                 &
          get_topog_mean, get_topog_stdev, &
          get_ocean_frac, get_ocean_mask,  &
          get_water_frac, get_water_mask,  &
          gaussian_topog_init, get_gaussian_topog

!-----------------------------------------------------------------------
! <NAMELIST NAME="topography_nml">
!   <DATA NAME="topog_file" TYPE="character" DEFAULT="DATA/navy_topography.data">
!       Name of topography file.
!   </DATA>
!   <DATA NAME="water_file" TYPE="character" DEFAULT="DATA/navy_pctwater.data">
!       Name of percent water file.
!   </DATA>

   character(len=128) :: topog_file = 'DATA/navy_topography.data', &
                         water_file = 'DATA/navy_pctwater.data'

   namelist /topography_nml/ topog_file, water_file
! </NAMELIST>

!-----------------------------------------------------------------------
! --- resolution of the topography data set ---
! <DATASET NAME="">
!   This module uses the 1/6 degree U.S. Navy mean topography
!   and percent water data sets.
!
!   These data sets have been re-formatted to separate 32-bit IEEE files.
!   The names of these files is specified by the <LINK SRC="#NAMELIST">namelist</LINK> input.
!
!The format for both files is as follows:
! <PRE>
!     record = 1    nlon, nlat
!     record = 2    blon, blat
!     record = 3    data
! </PRE>
!where:
! <PRE>
!     nlon, nlat = The number of longitude and latitude points
!                  in the horizontal grid.  For the 1/6 degree
!                  data sets this is 2160 x 1080. [integer]
!     blon, blat = The longitude and latitude grid box boundaries in degrees.
!                     [real(r8) :: blon(nlon+1), blat(nlat+1)]
!
!     data       = The topography or percent water data.
!                    [real(r8) :: data(nlon,nlat)]
! </PRE>
! </DATASET>
  integer :: unit
  integer :: ipts, jpts
  integer, parameter :: COMPUTE_STDEV = 123  ! use this flag to
                                             !   compute st dev

!-----------------------------------------------------------------------

 character(len=128) :: version = '$Revision$'
 character(len=128) :: tagname = '$Id$'

 logical :: module_is_initialized = .FALSE.

!-----------------------------------------------------------------------

 contains

!#######################################################################

   subroutine topography_init ()

     if ( module_is_initialized ) return

     call write_version_number (version,tagname)
     call read_namelist

     module_is_initialized = .TRUE.

   end subroutine topography_init

!#######################################################################

! <FUNCTION NAME="get_topog_mean">

!   <OVERVIEW>
!     Returns a "real(r8)istic" mean surface height field. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     Returns real(r8)istic mountains on a latitude-longtude grid.
!     The returned field is the mean topography for the given grid boxes.
!     Computed using a conserving area-weighted interpolation.
!     The current input data set is the 1/6 degree Navy mean topography.
!   </DESCRIPTION>
!   <TEMPLATE>
!     flag = <B>get_topog_mean</B> ( blon, blat, zmean )
!   </TEMPLATE>

!   <IN NAME="blon" TYPE="real(r8)" DIM="(:)">
!     The longitude (in radians) at grid box boundaries.
!   </IN>
!   <IN NAME="blat" TYPE="real(r8)" DIM="(:)">
!     The latitude (in radians) at grid box boundaries.
!   </IN>
!   <OUT NAME="zmean" UNIT=" meter" TYPE="real(r8)" DIM="(:,:)">
!     The mean surface height (meters).
!     The size of this field must be size(blon)-1 by size(blat)-1.
!   </OUT>
!   <OUT NAME="get_topog_mean" TYPE="logical">
!     A logical value of TRUE is returned if the surface height field
!     was successfully created. A value of FALSE may be returned if the
!     input topography data set was not readable.
!   </OUT>

!   <ERROR MSG="shape(zmean) is not equal to (/size(blon)-1,size(blat)-1/))" STATUS="FATAL">
!     Check the input grid size and output field size.
!   </ERROR>

 function get_topog_mean (blon, blat, zmean)

   real(r8), intent(in),  dimension(:)   :: blon, blat
   real(r8), intent(out), dimension(:,:) :: zmean
   logical :: get_topog_mean

!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(zmean) /= (/size(blon)-1,size(blat)-1/)) ) &
        call error_mesg('get_topog_mean','shape(zmean) is not&
            & equal to (/size(blon)-1,size(blat)-1/))', FATAL)

   if ( open_topog_file(topog_file) ) then
       call interp_topog ( blon, blat, zmean )
       get_topog_mean = .true.
   else
       get_topog_mean = .false.
   endif

!-----------------------------------------------------------------------

 end function get_topog_mean
! </FUNCTION>

!#######################################################################

! <FUNCTION NAME="get_topog_stdev">

!   <OVERVIEW>
!     Returns a standard deviation of higher resolution topography with 
!     the given model grid boxes. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     Returns the standard deviation of the "finer" input topography data set,
!     currently the Navy 1/6 degree mean topography data, within the
!     boundaries of the given input grid.
!   </DESCRIPTION>
!   <TEMPLATE>
!     flag = <B>get_topog_stdev</B> ( blon, blat, stdev )
!   </TEMPLATE>
!   <IN NAME="blon" TYPE="real(r8)" DIM="(:)">
!     The longitude (in radians) at grid box boundaries.
!   </IN>
!   <IN NAME="blat" TYPE="real(r8)" DIM="(:)">
!     The latitude (in radians) at grid box boundaries.
!   </IN>
!   <OUT NAME="stdev" UNITS="meter" TYPE="real(r8)" DIM="(:,:)">
!     The standard deviation of surface height (in meters) within
!     given input model grid boxes.
!     The size of this field must be size(blon)-1 by size(blat)-1.
!   </OUT>
!   <OUT NAME="get_topog_stdev" TYPE="logical">
!     A logical value of TRUE is returned if the output field was
!     successfully created. A value of FALSE may be returned if the
!     input topography data set was not readable.
!   </OUT>

 function get_topog_stdev (blon, blat, stdev)

   real(r8), intent(in),  dimension(:)   :: blon, blat
   real(r8), intent(out), dimension(:,:) :: stdev
   logical :: get_topog_stdev

   real(r8) :: zsurf (size(stdev,1),size(stdev,2))
!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(stdev) /= (/size(blon)-1,size(blat)-1/)) ) &
       call error_mesg('get_topog_stdev','shape(stdev) is not&
            & equal to (/size(blon)-1,size(blat)-1/))', FATAL)

   if ( open_topog_file(topog_file) ) then
       call interp_topog ( blon, blat, stdev, flag=COMPUTE_STDEV )
       get_topog_stdev = .true.
   else
       get_topog_stdev = .false.
   endif

!-----------------------------------------------------------------------

 end function get_topog_stdev
! </FUNCTION>

!#######################################################################

! <FUNCTION NAME="get_ocean_frac">

!   <OVERVIEW>
!      Returns fractional area covered by ocean in a grid box.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Returns fractional area covered by ocean in the given model grid boxes.
!   </DESCRIPTION>
!   <TEMPLATE>
!     flag = <B>get_ocean_frac</B> ( blon, blat, ocean_frac )
!   </TEMPLATE>

!   <IN NAME="blon" UNITS="radians" TYPE="real(r8)" DIM="(:)">
!     The longitude (in radians) at grid box boundaries.
!   </IN>
!   <IN NAME="blat" UNITS="radians" TYPE="real(r8)" DIM="(:)">
!     The latitude (in radians) at grid box boundaries.
!   </IN>
!   <OUT NAME="ocean_frac" TYPE="real(r8)" DIM="(:,:)">
!     The fractional amount (0 to 1) of ocean in a grid box.
!     The size of this field must be size(blon)-1 by size(blat)-1.
!   </OUT>
!   <OUT NAME="get_ocean_frac" TYPE="logical">
!     A logical value of TRUE is returned if the output field
!     was successfully created. A value of FALSE may be returned
!     if the Navy 1/6 degree percent water data set was not readable.
!   </OUT>

 function get_ocean_frac (blon, blat, ocean_frac)

 real(r8), intent(in),  dimension(:)   :: blon, blat
 real(r8), intent(out), dimension(:,:) :: ocean_frac
 logical :: get_ocean_frac

!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(ocean_frac) /= (/size(blon)-1,size(blat)-1/)) ) &
        call error_mesg('get_ocean_frac','shape(ocean_frac) is not&
                 & equal to (/size(blon)-1,size(blat)-1/))', FATAL)

   if ( open_topog_file(water_file) ) then
       call interp_water ( blon, blat, ocean_frac, do_ocean=.true. )
       get_ocean_frac = .true.
   else
       get_ocean_frac = .false.
   endif

!-----------------------------------------------------------------------

 end function get_ocean_frac
! </FUNCTION>

!#######################################################################
! <FUNCTION NAME="get_ocean_mask">

!   <OVERVIEW>
!     Returns a land-ocean mask in a grid box.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Returns a land-ocean mask in the given model grid boxes.
!   </DESCRIPTION>
!   <TEMPLATE>
!     flag = <B>get_ocean_mask</B> ( blon, blat, ocean_mask )
!   </TEMPLATE>

!   <IN NAME="blon" UNITS="radians" TYPE="real(r8)" DIM="(:)">
!     The longitude (in radians) at grid box boundaries.
!   </IN>
!   <IN NAME="blat" UNITS="radians" TYPE="real(r8)" DIM="(:)">
!     The latitude (in radians) at grid box boundaries.
!   </IN>
!   <OUT NAME="ocean_frac" TYPE="real(r8)" DIM="(:,:)">
!     The fractional amount (0 to 1) of ocean in a grid box.
!     The size of this field must be size(blon)-1 by size(blat)-1.
!   </OUT>
!   <OUT NAME="get_ocean_mask" TYPE="logical">
!     A logical value of TRUE is returned if the output field
!     was successfully created. A value of FALSE may be returned
!     if the Navy 1/6 degree percent water data set was not readable.
!   </OUT>

 function get_ocean_mask (blon, blat, ocean_mask)

 real(r8)   , intent(in),  dimension(:)   :: blon, blat
 logical, intent(out), dimension(:,:) :: ocean_mask
 logical :: get_ocean_mask

 real(r8), dimension(size(ocean_mask,1),size(ocean_mask,2)) :: ocean_frac
!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

 if ( get_ocean_frac(blon, blat, ocean_frac) ) then
   where (ocean_frac > 0.50)
     ocean_mask = .true.
   elsewhere
     ocean_mask = .false.
   end where
   get_ocean_mask = .true.
 else
   get_ocean_mask = .false.
 endif

!-----------------------------------------------------------------------

 end function get_ocean_mask
! </FUNCTION>

!#######################################################################
! <FUNCTION NAME="get_water_frac">

!   <OVERVIEW>
!     Returns fractional area covered by water.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Returns the percent of water in a grid box.
!   </DESCRIPTION>
!   <TEMPLATE>
!     flag = <B>get_water_frac</B> ( blon, blat, water_frac )
!   </TEMPLATE>

!   <IN NAME="blon" UNITS="radians" TYPE="real(r8)" DIM="(:)">
!     The longitude (in radians) at grid box boundaries.
!   </IN>
!   <IN NAME="blat" UNITS="radians" TYPE="real(r8)" DIM="(:)">
!     The latitude (in radians) at grid box boundaries.
!   </IN>
!   <OUT NAME="water_frac" TYPE="real(r8)" DIM="(:,:)">
!     The fractional amount (0 to 1) of water in a grid box.
!     The size of this field must be size(blon)-1 by size(blat)-1.
!   </OUT>
!   <OUT NAME="get_water_frac" TYPE="logical">
!     A logical value of TRUE is returned if the output field
!     was successfully created. A value of FALSE may be returned
!     if the Navy 1/6 degree percent water data set was not readable.
!   </OUT>
!   <ERROR MSG="shape(water_frac) is not equal to (/size(blon)-1,size(blat)-1/))" STATUS="FATAL">
!      Check the input grid size and output field size.
!   </ERROR>

 function get_water_frac (blon, blat, water_frac)

 real(r8), intent(in),  dimension(:)   :: blon, blat
 real(r8), intent(out), dimension(:,:) :: water_frac
 logical :: get_water_frac

!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

   if ( any(shape(water_frac) /= (/size(blon)-1,size(blat)-1/)) ) &
        call error_mesg('get_water_frac','shape(water_frac) is not&
                 & equal to (/size(blon)-1,size(blat)-1/))', FATAL)

   if ( open_topog_file(water_file) ) then
       call interp_water ( blon, blat, water_frac )
       get_water_frac = .true.
   else
       get_water_frac = .false.
   endif

!-----------------------------------------------------------------------

 end function get_water_frac
! </FUNCTION>

!#######################################################################
! <FUNCTION NAME="get_water_mask">

!   <OVERVIEW>
!     Returns a land-water mask in a grid box.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Returns a land-water mask in the given model grid boxes.
!   </DESCRIPTION>
!   <TEMPLATE>
!     flag = <B>get_water_mask</B> ( blon, blat, water_mask )
!   </TEMPLATE>

!   <IN NAME="blon" UNITS="radians" TYPE="real(r8)" DIM="(:)">
!     The longitude (in radians) at grid box boundaries.
!   </IN>
!   <IN NAME="blat" UNITS="radians" TYPE="real(r8)" DIM="(:)">
!     The latitude (in radians) at grid box boundaries.
!   </IN>
!   <OUT NAME="water_mask" TYPE="real(r8)" DIM="(:,:)">
!     A binary mask for water (true) or land (false).
!     The size of this field must be size(blon)-1 by size(blat)-1.
!   </OUT>
!   <OUT NAME="get_water_mask" TYPE="logical">
!     A logical value of TRUE is returned if the output field
!     was successfully created. A value of FALSE may be returned
!     if the Navy 1/6 degree percent water data set was not readable.
!   </OUT>

 function get_water_mask (blon, blat, water_mask)

 real(r8)   , intent(in),  dimension(:)   :: blon, blat
 logical, intent(out), dimension(:,:) :: water_mask
 logical :: get_water_mask

 real(r8), dimension(size(water_mask,1),size(water_mask,2)) :: water_frac
!-----------------------------------------------------------------------
   if (.not. module_is_initialized) call topography_init()

 if ( get_water_frac(blon, blat, water_frac) ) then
   where (water_frac > 0.50)
     water_mask = .true.
   elsewhere
     water_mask = .false.
   end where
   get_water_mask = .true.
 else
   get_water_mask = .false.
 endif

!-----------------------------------------------------------------------

 end function get_water_mask
! </FUNCTION>

!#######################################################################
!##################   private interfaces below here   ##################
!#######################################################################

 function open_topog_file ( filename )
 character(len=*), intent(in) :: filename
 logical :: open_topog_file

  if ( file_exist(filename) ) then
       unit = open_ieee32_file (trim(filename), 'read')
       read (unit) ipts, jpts
       open_topog_file = .true.
  else
       open_topog_file = .false.
  endif

 end function open_topog_file

!#######################################################################

 subroutine interp_topog ( blon, blat, zout, flag )
 real(r8)   , intent(in)  :: blon(:), blat(:)
 real(r8)   , intent(out) :: zout(:,:)
 integer, intent(in), optional :: flag

 real(r8) :: xdat(ipts+1), ydat(jpts+1)
 real(r8) :: zdat(ipts,jpts)
 real(r8) :: zout2(size(zout,1),size(zout,2))

! note: ipts,jpts,unit are global

    read (unit) xdat, ydat    ! read lon/lat edges in radians
    read (unit) zdat          ! read land surface height in meters
    call close_file (unit)

    call horiz_interp ( zdat, xdat, ydat, blon, blat, zout )

! compute standard deviation if necessary
    if (present(flag)) then
       if (flag == COMPUTE_STDEV) then
           zdat = zdat*zdat
           call horiz_interp ( zdat, xdat, ydat, blon, blat, zout2 )
           zout = zout2 - zout*zout
           where (zout > 0.0)
             zout = sqrt ( zout )
           elsewhere
             zout = 0.0
           endwhere
       endif
    endif

 end subroutine interp_topog

!#######################################################################

 subroutine interp_water ( blon, blat, zout, do_ocean )
 real(r8)   , intent(in)  :: blon(:), blat(:)
 real(r8)   , intent(out) :: zout(:,:)
 logical, intent(in), optional :: do_ocean

 real(r8) :: xdat(ipts+1), ydat(jpts+1), zdat(ipts,jpts)

! note: ipts,jpts,unit are global

    read (unit) xdat, ydat    ! read lon/lat edges in radians
    read (unit) zdat          ! read fractional water
    call close_file (unit)

! only use designated ocean points
    if (present(do_ocean)) then
        if (do_ocean) call determine_ocean_points (zdat)
    endif

! interpolate onto output grid
    call horiz_interp ( zdat, xdat, ydat, blon, blat, zout )

 end subroutine interp_water

!#######################################################################

 subroutine determine_ocean_points ( pctwater )
 real(r8), intent(inout) :: pctwater(:,:)
 logical :: ocean(size(pctwater,1),size(pctwater,2))
 integer :: i, j, m, n, im, ip, jm, jp, new 

 real(r8) :: ocean_pct_crit = .500

  ! resolution of the grid
    m = size(pctwater,1)
    n = size(pctwater,2)

  ! the 1/6 degree navy percent water data set
  ! designates ocean grid boxes as 100 percent water
  ! all other grid boxes have <= 99 percent water

  ! set a mask for ocean grid boxes
    ocean = (pctwater > .999) 
    new = count(ocean)

  ! set land grid boxes that have sufficient amount of water
  ! to ocean grid boxes when they are adjacent to ocean points
  ! iterate until there are no new ocean points
    do 
    if (new == 0) exit 
    new = 0 

       do j = 1, n
       do i = 1, m
          if (.not.ocean(i,j) .and. pctwater(i,j) > ocean_pct_crit) then
             im = i-1; ip = i+1; jm = j-1; jp = j+1
             if (im == 0)   im = m  
             if (ip == m+1) ip = 1
             if (jm == 0)   jm = 1
             if (jp == n+1) jp = n
           ! check the 8 grid boxes that surround this grid box
             if (ocean(im,j ) .or. ocean(ip,j ) .or. ocean(i ,jm) .or. ocean(i ,jp) .or. &
                 ocean(im,jm) .or. ocean(ip,jm) .or. ocean(ip,jp) .or. ocean(im,jp)) then
                 ocean(i,j) = .true.
                 new = new + 1
             endif
          endif
       enddo
       enddo
      !print *, 'new=',new

    enddo

  ! final step is to elimate water percentage if land
    where (.not.ocean) pctwater = 0.

 end subroutine determine_ocean_points

!#######################################################################
! reads the namelist file, write namelist to log file, 
! and initializes constants

subroutine read_namelist

   integer :: unit, ierr, io, iunit
   real(r8)    :: dtr, ght

!  read namelist

call find_namelist_in_file("input.nml", "topography_nml", iunit)
read(iunit, nml = topography_nml, iostat = io)
call check_namelist_read(iunit, io, "topography_nml")

! Old fms namelist read replace 8 June, 2006
!!!   if ( file_exist('input.nml')) then
!!!      unit = open_namelist_file ( )
!!!      ierr=1; do while (ierr /= 0)
!!!         read  (unit, nml=topography_nml, iostat=io, end=10)
!!!         ierr = check_nml_error(io,'topography_nml')
!!!      enddo
!!! 10   call close_file (unit)
!!!   endif

!  write version and namelist to log file

    write (stdlog(), nml=topography_nml)

end subroutine read_namelist

!#######################################################################

end module topography_mod

! <INFO>

!   <TESTPROGRAM NAME="">  
!  
!  To run this program you will need the topography and percent water
!  data sets and use the following namelist (in file input.nml).
!  
!   &amp;gaussian_topog_nml
!     height = 5000., 3000., 3000., 3000.,
!     olon   =   90.,  255.,  285.,    0.,
!     olat   =   45.,   45.,  -15.,  -90.,
!     wlon   =   15.,   10.,    5.,  180.,
!     wlat   =   15.,   25.,   25.,   20., /
!  
!  program test
!  
!  ! test program for topography and gaussian_topog modules
!  <PRE>  
!  use topography_mod
!  implicit none
!  
!  integer, parameter :: nlon=24, nlat=18
!  real(r8) :: x(nlon), y(nlat), xb(nlon+1), yb(nlat+1), z(nlon,nlat)
!  real(r8) :: hpi, rtd
!  integer :: i,j
!  logical :: a
!  
!  ! gaussian mountain parameters
!  real(r8), parameter :: ht=4000.
!  real(r8), parameter :: x0=90., y0=45. ! origin in degrees
!  real(r8), parameter :: xw=15., yw=15. ! half-width in degees
!  real(r8), parameter :: xr=30., yr= 0. ! ridge-width in degrees
!  
!  ! create lat/lon grid in radians
!    hpi = acos(0.0)
!    rtd = 90./hpi ! rad to deg
!    do i=1,nlon
!      xb(i) = 4.*hpi*real(r8)(i-1)/real(nlon)
!    enddo
!      xb(nlon+1) = xb(1)+4.*hpi
!      yb(1) = -hpi
!    do j=2,nlat
!      yb(j) = yb(j-1) + 2.*hpi/real(r8)(nlat)
!    enddo
!      yb(nlat+1) = hpi
!  ! mid-point of grid boxes
!    x(1:nlon) = 0.5*(xb(1:nlon)+xb(2:nlon+1))
!    y(1:nlat) = 0.5*(yb(1:nlat)+yb(2:nlat+1))
!  ! test topography_mod routines
!    a = get_topog_mean(xb,yb,z)
!    call printz ('get_topog_mean')
!  
!    a = get_water_frac(xb,yb,z)
!    z = z*100. ! in percent
!    call printz ('get_water_frac')
!  
!    a = get_ocean_frac(xb,yb,z)
!    z = z*100. ! in percent
!    call printz ('get_ocean_frac')
!  
!  ! test gaussian_topog_mod routines
!    a = .true.
!    z = get_gaussian_topog(x,y,ht,x0,y0,xw,yw,xr,yr)
!    call printz ('get_gaussian_topog')
!  
!    call gaussian_topog_init (x,y,z)
!    call printz ('gaussian_topog_init')
!  
!  contains
!  
!  ! simple printout of topog/water array
!    subroutine printz (lab)
!    character(len=*), intent(in) :: lab
!     if (a) then
!        print '(/a)', trim(lab)
!     else
!        print '(/a)', 'no data available: '//trim(lab)
!        return
!     endif
!      ! print full grid
!        print '(3x,25i5)', (nint(x(i)*rtd),i=1,nlon)
!      do j=nlat,1,-1
!        print '(i3,25i5)', nint(y(j)*rtd), (nint(z(i,j)),i=1,nlon)
!      enddo
!    end subroutine printz
!  
!  end program test
!   </PRE>
!   </TESTPROGRAM>

!   <BUG>                  
!      Water mask produces some possible erroneous water points along
!      the coast of Antarctic (at about 90W).
!   </BUG>

!   <FUTURE>Use of netcdf data sets. </FUTURE>
!   <FUTURE>Incorporate other topography and ocean data sets. </FUTURE>
! 
! </INFO>
