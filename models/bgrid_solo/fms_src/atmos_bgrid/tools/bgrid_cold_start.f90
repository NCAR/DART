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

module bgrid_cold_start_mod

use types_mod, only : r8
use  bgrid_horiz_mod, only: horiz_grid_type, get_horiz_grid_bound, TGRID
use   bgrid_halo_mod, only: update_halo, TEMP
use   topography_mod, only: gaussian_topog_init, get_topog_mean
use utilities_mod, only : find_namelist_in_file, check_namelist_read
use          fms_mod, only: file_exist, open_namelist_file, open_restart_file, &
                            open_ieee32_file, close_file,     &
                            stdlog, error_mesg, FATAL, write_version_number
use    constants_mod, only: GRAV, RDGAS

implicit none
private

public cold_start_resol, cold_start

!-----------------------------------------------------------------------
!----- namelist ------
!
!  nlon = number of grid points along the longitude axis (1st dimension)
!  nlat = number of grid points along the latitude axis (2nd dimension)
!  nlev = number of vertical levels (equally spaced in sigma)
!  pref = initial surface pressure in pascals
!  tref = initial temperature in deg kelvin
!
!  NOTE: nlon and nlat are for the global compute grid
!

integer :: nlon = 0
integer :: nlat = 0
integer :: nlev = 0
   real(r8) :: pref = 1000.e2
   real(r8) :: tref = 255.
logical :: equal_vert_spacing = .true.

namelist /bgrid_cold_start_nml/ nlon, nlat, nlev, pref, tref, equal_vert_spacing

!-----------------------------------------------------------------------

character(len=128) :: version = '$Revision$'
character(len=128) :: tag = '$Id$'

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine cold_start_resol ( nx, ny, nz )

      integer, intent(out) :: nx, ny, nz

      integer :: unit, ierr, io, iunit

!-------------- read namelist --------------
! Old fms version replaced with dart 8 June, 2006
call find_namelist_in_file("input.nml", "bgrid_cold_start_nml", iunit)
read(iunit, nml = bgrid_cold_start_nml, iostat = io)
call check_namelist_read(iunit, io, "bgrid_cold_start_nml")


!!!   if ( file_exist('input.nml')) then
!!!      unit = open_namelist_file ( )
!!!      ierr=1; do while (ierr /= 0)
!!!         read  (unit, nml=bgrid_cold_start_nml, iostat=io, end=10)
!!!         ierr = check_nml_error(io,'bgrid_cold_start_nml')
!!!      enddo
!!! 10   call close_file (unit)
!!!   endif

!-------- write version and namelist to log file --------

   call write_version_number (version,tag)
   write (stdlog(), nml=bgrid_cold_start_nml)

!------- must specify a resolution -----

   if (nlon == 0 .or. nlat == 0 .or. nlev == 0)  &
   call error_mesg ('bgrid_cold_start_mod', 'resolution not specified', FATAL)

!------- otherwise, return resolution to calling program ------

   nx  = nlon
   ny  = nlat
   nz  = nlev

end subroutine cold_start_resol

!#######################################################################

subroutine cold_start ( Hgrid, eta, peta, fis, res, ps, pssl, u, v, t )

   type(horiz_grid_type), intent(inout) :: Hgrid
   real(r8), intent(out), dimension(:)      :: eta, peta
   real(r8), intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:)     :: fis, res, &
                                                              ps, pssl
   real(r8), intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:)   :: u, v, t
  !real(r8), intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:,:) :: r
   integer :: m
   real(r8)    :: rt

!   very simple initial condition
!   -----------------------------
!       equally spaced levels
!       no topography
!       uniform pressure
!       isothermal
!       no wind


!--- no hybrid option ---
   peta = 0.0

   eta = compute_sigma (size(eta)-1)
   fis = compute_topog (Hgrid)

!--- no eta option ---
   res = 1.0

   rt = RDGAS * tref

     ps   = pref*exp(-fis/rt)
     pssl = ps                ! no eta option
     u    = 0.0
     v    = 0.0
     t    = tref

end subroutine cold_start

!#######################################################################

 function compute_topog (Hgrid) result (fis)

   type(horiz_grid_type), intent(inout) :: Hgrid
   real(r8), dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: fis

   real(r8) :: hpi, dtr
   integer :: unit, ires, jres
!-----------------------------------------------------------------------

   hpi = acos(0.0)
   dtr = hpi/90.

   if ( file_exist('INPUT/topography.res') ) then

   !  OPTION 1:
   !  read topography from restart file if this file exists

      unit = open_restart_file ('INPUT/topography.res', action='read')
      read (unit) ires,jres
      if ( ires /= Hgrid%nlon .or. jres /= Hgrid%nlat )  &
                      call error_mesg ('bgrid_cold_start_mod',  &
                  'incorrect resolution in file topography.res', FATAL)
      !!!call set_domain (Hgrid%Tmp%Domain)
      !!!call read_data  (unit, fis)
      write(*, *) 'In compute_topog in bgrid_cold_start'
      write(*, *) 'Need to read in topogrphy only onto compute domain'
      write(*, *) 'Should be simple?'
      if(1 == 1) stop
      !!!read(unit) fis(:, :)
      call close_file (unit)

   else

   !  OPTION 2: interp real topography from data set, otherwise
   !  OPTION 3: generate gaussian-shaped mountains

      if ( .not.compute_real_topog(Hgrid,fis) ) then
          call gaussian_topog_init ( Hgrid%Tmp%alm(:,Hgrid%Tmp%js), &
                                     Hgrid%Tmp%aph(Hgrid%Tmp%is,:), &
                                     fis(:,:) )
        ! convert topog from meters to geop ht 
          fis = fis * GRAV
        ! only in case of roundoff error
          where (fis(:,:) < 0.0) fis(:,:) = 0.0
      endif

   endif

   call update_halo ( Hgrid, TEMP, fis )

!-----------------------------------------------------------------------

 end function compute_topog

!#######################################################################

 function compute_real_topog ( Hgrid, fis )

   type(horiz_grid_type), intent(in)  :: Hgrid
   real(r8),                  intent(out) :: fis(Hgrid%ilb:,Hgrid%jlb:)
   logical :: compute_real_topog

   real(r8), dimension(Hgrid%Tmp%is:Hgrid%Tmp%ie+1) :: blon
   real(r8), dimension(Hgrid%Tmp%js:Hgrid%Tmp%je+1) :: blat
!-----------------------------------------------------------------------

   call get_horiz_grid_bound (Hgrid, TGRID, blon, blat, global=.false.)

   compute_real_topog = get_topog_mean ( blon, blat, &
           fis(Hgrid%Tmp%is:Hgrid%Tmp%ie,Hgrid%Tmp%js:Hgrid%Tmp%je) )

!  convert topog from meters to geop ht 
   if (compute_real_topog) fis = fis * GRAV

!-----------------------------------------------------------------------

 end function compute_real_topog

!#######################################################################

 function compute_sigma (nlev) result (eta)

   integer, intent(in) :: nlev
   real(r8), dimension(nlev+1) :: eta
   real(r8) :: dz, qk
   integer :: k

     eta(1) = 0.0
     eta(nlev+1) = 1.0

   if ( equal_vert_spacing ) then

       dz = 1./nlev
       do k = 2, nlev
         eta(k) = eta(k-1) + dz
       enddo

   else

       do k = 1, nlev-1
         qk = (2.0*k)/(2.0*nlev)
         eta(k+1) = qk*qk*(3.-2.*qk)
       enddo

   endif

 end function compute_sigma

!#######################################################################

end module bgrid_cold_start_mod

