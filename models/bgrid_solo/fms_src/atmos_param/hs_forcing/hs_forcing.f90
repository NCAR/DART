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

module hs_forcing_mod

!-----------------------------------------------------------------------

use types_mod, only : r8
use     constants_mod, only: KAPPA, CP, GRAV

use utilities_mod, only : find_namelist_in_file, check_namelist_read
use           fms_mod, only: error_mesg, FATAL, file_exist,       &
                             open_namelist_file, &
                             close_file,     &
                             write_version_number, stdlog,        &
                             uppercase

use  time_manager_mod, only: time_type

!use  diag_manager_mod, only: register_diag_field, send_data

use  field_manager_mod, only: MODEL_ATMOS, parse
use tracer_manager_mod, only: query_method, get_tracer_index

implicit none
private

!-----------------------------------------------------------------------
!---------- interfaces ------------

   public :: hs_forcing, hs_forcing_init

!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------

   real(r8) :: t_zero=315., t_strat=200., delh=60., delv=10., eps=0., sigma_b=0.7
   real(r8) :: P00 = 1.e5

   real(r8) :: ka = -40. !  negative values are damping time in days
   real(r8) :: ks =  -4., kf = -1.

   logical :: do_conserve_energy = .true.

   real(r8) :: trflux = 1.e-5   !  surface flux for optional tracer
   real(r8) :: trsink = -4.     !  damping time for tracer

!-----------------------------------------------------------------------

   namelist /hs_forcing_nml/  t_zero, t_strat, delh, delv, eps, &
                              sigma_b, ka, ks, kf, do_conserve_energy, &
                              trflux, trsink

!-----------------------------------------------------------------------

   character(len=128) :: version='$Revision$'
   character(len=128) :: tag='$Id$'

   real(r8) :: tka, tks, vkf
   real(r8) :: trdamp

   integer :: id_teq, id_tdt, id_udt, id_vdt,  &
              id_tdt_diss, id_diss_heat
   real(r8)    :: missing_value = -1.e10
   character(len=14) :: mod_name = 'hs_forcing'

   logical :: do_init = .true.

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine hs_forcing ( is, ie, js, je, dt, Time, lat, p_half, p_full, &
                         u, v, t, r, um, vm, tm, rm, udt, vdt, tdt, rdt,&
                         mask, kbot )

!-----------------------------------------------------------------------
   integer, intent(in)                        :: is, ie, js, je
      real(r8), intent(in)                        :: dt
 type(time_type), intent(in)                  :: Time
      real(r8), intent(in),    dimension(:,:)     :: lat
      real(r8), intent(in),    dimension(:,:,:)   :: p_half, p_full
      real(r8), intent(in),    dimension(:,:,:)   :: u, v, t, um, vm, tm
      real(r8), intent(in),    dimension(:,:,:,:) :: r, rm
      real(r8), intent(inout), dimension(:,:,:)   :: udt, vdt, tdt
      real(r8), intent(inout), dimension(:,:,:,:) :: rdt

      real(r8), intent(in),    dimension(:,:,:), optional :: mask
   integer, intent(in),    dimension(:,:)  , optional :: kbot
!-----------------------------------------------------------------------
   real(r8), dimension(size(t,1),size(t,2))           :: ps, diss_heat
   real(r8), dimension(size(t,1),size(t,2),size(t,3)) :: ttnd, utnd, vtnd, teq, pmass
   real(r8), dimension(size(r,1),size(r,2),size(r,3)) :: rst, rtnd
   integer :: i, j, k, kb, n
   logical :: used
   real(r8)    :: flux, sink, value
   character(len=128) :: scheme, params

!-----------------------------------------------------------------------

     if (do_init) call error_mesg ('hs_forcing','hs_forcing_init has not been called', FATAL)

!-----------------------------------------------------------------------
!     surface pressure

     if (present(kbot)) then
         do j=1,size(p_half,2)
         do i=1,size(p_half,1)
            kb = kbot(i,j)
            ps(i,j) = p_half(i,j,kb+1)
         enddo
         enddo
     else
            ps(:,:) = p_half(:,:,size(p_half,3))
     endif

!-----------------------------------------------------------------------
!     rayleigh damping of wind components near the surface

      ! unnecessary extra test on the optional mask argument; this code
      ! is a work-around for a gfortran compiler bug; it was incorrectly
      ! passing down garbage as the optional mask argument if not specified.
      if (present(mask)) then
        call rayleigh_damping ( ps, p_full, u, v, utnd, vtnd, mask )
      else
        call rayleigh_damping ( ps, p_full, u, v, utnd, vtnd )
      endif

      if (do_conserve_energy) then
         ttnd = -((um+.5*utnd*dt)*utnd + (vm+.5*vtnd*dt)*vtnd)/CP
         tdt = tdt + ttnd
!         if (id_tdt_diss > 0) used = send_data ( id_tdt_diss, ttnd, Time, is, js)
       ! vertical integral of ke dissipation
         if ( id_diss_heat > 0 ) then
          do k = 1, size(t,3)
            pmass(:,:,k) = p_half(:,:,k+1)-p_half(:,:,k)
          enddo
          diss_heat = CP/GRAV * sum( ttnd*pmass, 3)
!          used = send_data ( id_diss_heat, diss_heat, Time, is, js)
         endif
      endif

      udt = udt + utnd
      vdt = vdt + vtnd

!      if (id_udt > 0) used = send_data ( id_udt, utnd, Time, is, js)
!      if (id_vdt > 0) used = send_data ( id_vdt, vtnd, Time, is, js)

!-----------------------------------------------------------------------
!     thermal forcing for held & suarez (1994) benchmark calculation

      ! unnecessary extra test on the optional mask argument; this code
      ! is a work-around for a gfortran compiler bug; it was incorrectly
      ! passing down garbage as the optional mask argument if not specified.
      if (present(mask)) then
        call newtonian_damping ( lat, ps, p_full, t, ttnd, teq, mask )
      else
        call newtonian_damping ( lat, ps, p_full, t, ttnd, teq )
      endif

      tdt = tdt + ttnd

!      if (id_tdt > 0) used = send_data ( id_tdt, ttnd, Time, is, js)
!      if (id_teq > 0) used = send_data ( id_teq, teq,  Time, is, js)

!-----------------------------------------------------------------------
!     -------- tracers -------

      do n = 1, size(rdt,4)
         flux = trflux
         sink = trsink
         if (query_method('tracer_sms', MODEL_ATMOS, n, scheme, params)) then
             if (uppercase(trim(scheme)) == 'NONE') cycle
             if (parse(params,'flux',value) == 1) flux = value
             if (parse(params,'sink',value) == 1) sink = value
         endif
         rst = rm(:,:,:,n) + dt*rdt(:,:,:,n)
         call tracer_source_sink ( flux, sink, p_half, rst, rtnd, kbot )
         rdt(:,:,:,n) = rdt(:,:,:,n) + rtnd
      enddo

!-----------------------------------------------------------------------

 end subroutine hs_forcing

!#######################################################################

 subroutine hs_forcing_init ( axes, Time )

!-----------------------------------------------------------------------
!
!           routine for initializing the model with an
!              initial condition at rest (u & v = 0)
!
!-----------------------------------------------------------------------

           integer, intent(in) :: axes(4)
   type(time_type), intent(in) :: Time

!-----------------------------------------------------------------------
   integer  unit, io, ierr, iunit

!     ----- read namelist -----
! fms replaced with dart, 8 June, 2006
call find_namelist_in_file("input.nml", "hs_forcing_nml", iunit)
read(iunit, nml = hs_forcing_nml, iostat = io)
call check_namelist_read(iunit, io, "hs_forcing_nml")

!!!      if (file_exist('input.nml')) then
!!!         unit = open_namelist_file ( )
!!!         ierr=1; do while (ierr /= 0)
!!!            read  (unit, nml=hs_forcing_nml, iostat=io, end=10)
!!!            ierr = check_nml_error (io, 'hs_forcing_nml')
!!!         enddo
!!!  10     call close_file (unit)
!!!      endif

!     ----- write version info and namelist to log file -----

      call write_version_number (version,tag)
      write (stdlog(),nml=hs_forcing_nml)


!     ----- compute coefficients -----

      if (ka < 0.) ka = -86400.*ka
      if (ks < 0.) ks = -86400.*ks
      if (kf < 0.) kf = -86400.*kf

      tka = 0.; if (ka > 0.) tka = 1./ka
      tks = 0.; if (ks > 0.) tks = 1./ks
      vkf = 0.; if (kf > 0.) vkf = 1./kf

!     ----- for tracers -----

      if (trsink < 0.) trsink = -86400.*trsink
      trdamp = 0.; if (trsink > 0.) trdamp = 1./trsink

!     ----- register diagnostic fields -----

!      id_teq = register_diag_field ( mod_name, 'teq', axes(1:3), Time, &
!                      'equilibrium temperature', 'deg_K'   , &
!                      missing_value=missing_value, range=(/100.,400./) )

!      id_tdt = register_diag_field ( mod_name, 'tdt_ndamp', axes(1:3), Time, &
!                      'newtonian damping', 'deg_K/sec' ,    &
!                       missing_value=missing_value     )

!      id_udt = register_diag_field ( mod_name, 'udt_rdamp', axes(1:3), Time, &
!                      'rayleigh damping (zonal wind)', 'm/s2',       &
!                       missing_value=missing_value     )

!      id_vdt = register_diag_field ( mod_name, 'vdt_rdamp', axes(1:3), Time, &
!                      'rayleigh damping (meridional wind)', 'm/s2',  &
!                       missing_value=missing_value     )

!      if (do_conserve_energy) then
!         id_tdt_diss = register_diag_field ( mod_name, 'tdt_diss_rdamp', axes(1:3), &
!                   Time, 'Dissipative heating from Rayleigh damping', 'deg_K/sec',&
!                   missing_value=missing_value     )

!         id_diss_heat = register_diag_field ( mod_name, 'diss_heat_rdamp', axes(1:2), &
!                   Time, 'Integrated dissipative heating for Rayleigh damping', 'W/m2')
!      endif

      do_init  = .false.

!-----------------------------------------------------------------------

 end subroutine hs_forcing_init

!#######################################################################

 subroutine hs_forcing_end 

!-----------------------------------------------------------------------
!
!       routine for terminating held-suarez benchmark module
!             (this routine currently does nothing)
!
!-----------------------------------------------------------------------


 end subroutine hs_forcing_end

!#######################################################################

 subroutine newtonian_damping ( lat, ps, p_full, t, tdt, teq, mask )

!-----------------------------------------------------------------------
!
!   routine to compute thermal forcing for held & suarez (1994)
!   benchmark calculation.
!
!-----------------------------------------------------------------------

real(r8), intent(in),  dimension(:,:)   :: lat, ps
real(r8), intent(in),  dimension(:,:,:) :: p_full, t
real(r8), intent(out), dimension(:,:,:) :: tdt, teq
real(r8), intent(in),  dimension(:,:,:), optional :: mask

!-----------------------------------------------------------------------

          real(r8), dimension(size(t,1),size(t,2)) :: &
     sin_lat, sin_lat_2, cos_lat_2, t_star, cos_lat_4, &
     tstr, sigma, the, tfactr, rps, p_norm

       real(r8), dimension(size(t,1),size(t,2),size(t,3)) :: tdamp

       integer :: k
       real(r8)    :: tcoeff, pref

!-----------------------------------------------------------------------
!------------latitudinal constants--------------------------------------

      sin_lat  (:,:) = sin(lat(:,:))
      sin_lat_2(:,:) = sin_lat(:,:)*sin_lat(:,:)
      cos_lat_2(:,:) = 1.0-sin_lat_2(:,:)
      cos_lat_4(:,:) = cos_lat_2(:,:)*cos_lat_2(:,:)

      t_star(:,:) = t_zero - delh*sin_lat_2(:,:) - eps*sin_lat(:,:)
      tstr  (:,:) = t_strat - eps*sin_lat(:,:)

!-----------------------------------------------------------------------

      tcoeff = (tks-tka)/(1.0-sigma_b)
      pref = P00
      rps  = 1./ps

      do k = 1, size(t,3)

!  ----- compute equilibrium temperature (teq) -----
         p_norm(:,:) = p_full(:,:,k)/pref
         the   (:,:) = t_star(:,:) - delv*cos_lat_2(:,:)*log(p_norm(:,:))
         teq(:,:,k) = the(:,:)*(p_norm(:,:))**KAPPA
         teq(:,:,k) = max( teq(:,:,k), tstr(:,:) )

!  ----- compute damping -----
         sigma(:,:) = p_full(:,:,k)*rps(:,:)
         where (sigma(:,:) <= 1.0 .and. sigma(:,:) > sigma_b)
           tfactr(:,:) = tcoeff*(sigma(:,:)-sigma_b)
           tdamp(:,:,k) = tka + cos_lat_4(:,:)*tfactr(:,:)
         elsewhere
           tdamp(:,:,k) = tka
         endwhere

      enddo

!*** note: if the following loop uses vector notation for all indices
!          then the code will not run ??????

      do k=1,size(t,3)
         tdt(:,:,k) = -tdamp(:,:,k)*(t(:,:,k)-teq(:,:,k))
      enddo

      if (present(mask)) then
         tdt = tdt * mask
         teq = teq * mask
      endif

!-----------------------------------------------------------------------

 end subroutine newtonian_damping

!#######################################################################

 subroutine rayleigh_damping ( ps, p_full, u, v, udt, vdt, mask )

!-----------------------------------------------------------------------
!
!           rayleigh damping of wind components near surface
!
!-----------------------------------------------------------------------

real(r8), intent(in),  dimension(:,:)   :: ps
real(r8), intent(in),  dimension(:,:,:) :: p_full, u, v
real(r8), intent(out), dimension(:,:,:) :: udt, vdt
real(r8), intent(in),  dimension(:,:,:), optional :: mask

!-----------------------------------------------------------------------

real(r8), dimension(size(u,1),size(u,2)) :: sigma, vfactr, rps

integer :: i,j,k
real(r8)    :: vcoeff

!-----------------------------------------------------------------------
!----------------compute damping----------------------------------------

      vcoeff = -vkf/(1.0-sigma_b)
      rps = 1./ps

      do k = 1, size(u,3)

         sigma(:,:) = p_full(:,:,k)*rps(:,:)

         where (sigma(:,:) <= 1.0 .and. sigma(:,:) > sigma_b)
            vfactr(:,:) = vcoeff*(sigma(:,:)-sigma_b)
            udt(:,:,k)  = vfactr(:,:)*u(:,:,k)
            vdt(:,:,k)  = vfactr(:,:)*v(:,:,k)
         elsewhere
            udt(:,:,k) = 0.0
            vdt(:,:,k) = 0.0
         endwhere
      enddo

      if (present(mask)) then
          udt = udt * mask
          vdt = vdt * mask
      endif

!-----------------------------------------------------------------------

 end subroutine rayleigh_damping

!#######################################################################

 subroutine tracer_source_sink ( flux, damp, p_half, r, rdt, kbot )

!-----------------------------------------------------------------------
      real(r8), intent(in)  :: flux, damp, p_half(:,:,:), r(:,:,:)
      real(r8), intent(out) :: rdt(:,:,:)
   integer, intent(in), optional :: kbot(:,:)
!-----------------------------------------------------------------------
      real(r8), dimension(size(r,1),size(r,2),size(r,3)) :: source, sink
      real(r8), dimension(size(r,1),size(r,2))           :: pmass

      integer :: i, j, kb
      real(r8)    :: rdamp
!-----------------------------------------------------------------------

      rdamp = damp
      if (rdamp < 0.) rdamp = -86400.*rdamp   ! convert days to seconds
      if (rdamp > 0.) rdamp = 1./rdamp

!------------ simple surface source and global sink --------------------

      source(:,:,:)=0.0

   if (present(kbot)) then
      do j=1,size(r,2)
      do i=1,size(r,1)
         kb = kbot(i,j)
         pmass (i,j)    = p_half(i,j,kb+1) - p_half(i,j,kb)
         source(i,j,kb) = flux/pmass(i,j)
      enddo
      enddo
   else
         kb = size(r,3)
         pmass (:,:)    = p_half(:,:,kb+1) - p_half(:,:,kb)
         source(:,:,kb) = flux/pmass(:,:)
   endif

     sink(:,:,:) = rdamp*r(:,:,:)
     rdt(:,:,:) = source(:,:,:)-sink(:,:,:)

!-----------------------------------------------------------------------

 end subroutine tracer_source_sink

!#######################################################################

end module hs_forcing_mod
