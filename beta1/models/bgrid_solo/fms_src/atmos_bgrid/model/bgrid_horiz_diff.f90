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

module bgrid_horiz_diff_mod

!-----------------------------------------------------------------------

use types_mod, only : r8
use bgrid_horiz_mod      , only: horiz_grid_type
use bgrid_vert_mod       , only: vert_grid_type, compute_pres_full
use bgrid_masks_mod      , only: grid_mask_type
use bgrid_prog_var_mod   , only: prog_var_type, var_init
use bgrid_halo_mod       , only: update_halo, vel_flux_boundary, &
                                 EAST, NORTH, NOPOLE, POLEONLY,  &
                                 TEMP, UWND, VWND
use bgrid_change_grid_mod, only: mass_to_vel

use         fms_mod, only:  error_mesg, FATAL, write_version_number, &
                            uppercase, stdlog
use   constants_mod, only:  RADIUS

use  field_manager_mod, only: MODEL_ATMOS, parse
use tracer_manager_mod, only: query_method, get_tracer_names, get_number_tracers

implicit none
private

!=======================================================================
!
!                    linear horizontal mixing
!               with option for any order accuracy
!
!=======================================================================

public hdiff_control_type

type hdiff_control_type
     type(horiz_grid_type), pointer :: Hgrid
     integer, pointer :: order(:)
     real(r8)   , pointer :: coeff(:), slope(:,:)
     logical       :: do_damping, do_slope_adj_temp, do_slope_adj_wind
     integer       :: nplevs
     integer       :: num_adjust_dt
     integer       :: damping_scheme
     real(r8), dimension(:,:),   pointer :: areahx, areahy, &
                                        areavx, areavy, &
                                        wth,  wtv
end type hdiff_control_type

!-----------------------------------------------------------------------
!--------- public interfaces ----------

public   horiz_diff, horiz_diff_init

!-----------------------------------------------------------------------
!--------- private data ----------

   integer  :: nlev

!  ---- slope correction weights ----
!  real(r8), dimension(3) :: slope_weights = (/ 1.00, 0.75, 0.25 /)
!!!!!!!!!!!!!!!!!!!!!   calgary values = (/ 1.00, 0.50, 0.10 /)
!!!!!!!!!!!!!!!!!!!!!    bombay values = (/ 1.00, 1.00, 0.50 /)

 character(len=128) :: version='$Revision$'
 character(len=128) :: tagname='$Id$'
 logical :: do_log = .true.

!  timing data

 integer, parameter :: time_level = 5
 integer :: id_total
 logical :: do_clock_init = .true.

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine horiz_diff ( Control, Masks, dt, dpde, pres, Var, Var_dt )

!-----------------------------------------------------------------------
!
!   dt     = adjustment time step
!   dpde   = pressure weight for model layers
!   pres   = pressure at full model levels
!   Var    = prognostic variables at the last updated time level
!   Var_dt = tendency of prognostic variables since the last
!            updated time level
!
!-----------------------------------------------------------------------

type(hdiff_control_type), intent(inout) :: Control
type (grid_mask_type),    intent(in)    :: Masks

  real(r8),                 intent(in)    :: dt
  real(r8),                 intent(in)    :: dpde(Control%Hgrid%ilb:,Control%Hgrid%jlb:,:),&
                                         pres(Control%Hgrid%ilb:,Control%Hgrid%jlb:,:)
  type (prog_var_type), intent(in)    :: Var
  type (prog_var_type), intent(inout) :: Var_dt

!-----------------------------------------------------------------------

  real(r8), dimension (Control%Hgrid%ilb:Control%Hgrid%iub,    &
                   Control%Hgrid%jlb:Control%Hgrid%jub) :: &
                                         hkew2, hkns2, coeff

  real(r8), dimension(Control%Hgrid%ilb:Control%Hgrid%iub,       &
                  Control%Hgrid%jlb:Control%Hgrid%jub,       &
                  size(dpde,3)) :: hkew3, hkns3, hkew, hkns, &
                                   hcew, hcns, dat, vdat, hdac

  real(r8), dimension(Control%Hgrid%ilb:Control%Hgrid%iub, &
                  Control%Hgrid%jlb:Control%Hgrid%jub, &
                     size(dpde,3),3) :: pterms

  real(r8)    :: dt_inv, hsign
  integer :: i, j, k, n, is, ie, js, je, nplev, ntp

!=======================================================================
!  --- should damping be done ??? ----

   if ( .not. Control%do_damping ) return

!-----------------------------------------------------------------------
!       --- check the horizontal dimensions of the input array ---

      nlev = size(dpde,3)

    if (size(dpde,1) /= Control%Hgrid%isize .or.               &
        size(dpde,2) /= Control%Hgrid%jsize ) call error_mesg  &
       ('bgrid_horiz_diff', 'input array has the wrong dimensions.', FATAL)

!   ---- time step related values ----

    dt_inv  = 1./dt

!   ---- do not use pressure weight at all levels ----

    nplev = Control%nplevs

!   ---- initialize flux weights ----

    hkew2 = 0.0;  hkns2 = 0.0
    hkew3 = 0.0;  hkns3 = 0.0

!-----------------------------------------------------------------------
!-------------setup temperature and tracer diffusion -------------------
!-----------------------------------------------------------------------
   ntp = count( Control%order(1:Var_dt%ntrace) > 0 )
   if (Control%order(0) > 0 .or. ntp > 0) then

      is = Control%Hgrid%Tmp%is-(Control%Hgrid%ihalo-1)
      ie = Control%Hgrid%Tmp%ie+(Control%Hgrid%ihalo-1)
      js = Control%Hgrid%Tmp%js-(Control%Hgrid%jhalo-1)
      je = Control%Hgrid%Tmp%je+(Control%Hgrid%jhalo-1)

      hdac = spread(Control%Hgrid%Tmp%rarea,3,nlev)
      do k = nplev+1, nlev
         hdac(:,:,k) = hdac(:,:,k) / (2.0*dpde(:,:,k))
      enddo

!     ----- mass-weighted fluxes -----

      hkew3(is-1:ie,js-1:je,1:nplev) = &
                       spread(Control%areahx(is-1:ie,js-1:je),3,nplev)
      hkns3(is-1:ie,js-1:je,1:nplev) = &
                       spread(Control%areahy(is-1:ie,js-1:je),3,nplev)
!dir$ IVDEP
      do k = nplev+1, nlev
      do j = js-1, je
      do i = is-1, ie
         hkew3(i,j,k) = Control%areahx(i,j) * (dpde(i,j,k)+dpde(i+1,j,k))
         hkns3(i,j,k) = Control%areahy(i,j) * (dpde(i,j,k)+dpde(i,j+1,k))
      enddo
      enddo
      enddo

!     ----- mask out fluxes below step-mountain ----

    if (.not.Masks%sigma) then
!dir$ IVDEP
       do k = 1, nlev
       do j = js-1, je
       do i = is-1, ie
          hkew3(i,j,k) = hkew3(i,j,k)*Masks%Tmp%mask(i,j,k)*Masks%Tmp%mask(i+1,j,k)
          hkns3(i,j,k) = hkns3(i,j,k)*Masks%Tmp%mask(i,j,k)*Masks%Tmp%mask(i,j+1,k)
       enddo
       enddo
       enddo
    endif

!-----------------------------------------------------------------------
!----- slope adjustment ------------------------------------------------

      if ( Control % do_slope_adj_temp ) then
         call slope_correction_init ( Control%Hgrid,      &
                                      Control%nplevs,     &
                                      pres, pterms        )
      else
         hcew = 0.0
         hcns = 0.0
      endif

!-----------------------------------------------------------------------
!------------------temperature diffusion--------------------------------

      if (Control%order(0) > 0) then

         hkew2 = Control%coeff(0)*Control%wth
         hkns2 = Control%coeff(0)*Control%wth
!        -- stability condition --
         if ( Control%damping_scheme > 1) then
            hkew2 = min(0.125_r8,hkew2)
            hkns2 = min(0.125_r8,hkew2)
         endif
         hkew = hkew3 * spread(hkew2,3,nlev)
         hkns = hkns3 * spread(hkns2,3,nlev)

         dat  = Var%t + dt*Var_dt%t
!        --- slope adjustent ---
         if ( Control%do_slope_adj_temp ) &
         call slope_correction ( Control%Hgrid, Control%nplevs, &
                                 Control%slope(:,0), pterms,    &
                                 dat, hcew, hcns )
         call diff_mass (Control%Hgrid, dat, hcew, hcns, hkew, hkns, hdac,  &
                         Control%order(0))
         hsign = 1.; if (mod(Control%order(0),4) == 0) hsign = -1.
         Var_dt%t = Var_dt%t + hsign * dt_inv * dat

      endif

!-----------------------------------------------------------------------
!------------ tracer diffusion (prognostic tracers only) ---------------

      if (ntp > 0) then

         do n = 1, Var_dt%ntrace

           !--- set up diff coeff for tracers ---
            hkew2 = Control%coeff(n)*Control%wth
            hkns2 = Control%coeff(n)*Control%wth
           !-- stability condition --
            if ( Control%damping_scheme > 1) then
               hkew2 = min(0.125_r8,hkew2)
               hkns2 = min(0.125_r8,hkns2)
            endif
            hkew = hkew3 * spread(hkew2,3,nlev)
            hkns = hkns3 * spread(hkns2,3,nlev)

            if (Control%order(n) == 0) cycle
            dat  = Var%r(:,:,:,n) + dt*Var_dt%r(:,:,:,n)
           !--- slope adjustent ---
            if ( Control%do_slope_adj_temp ) &
            call slope_correction ( Control%Hgrid, Control%nplevs, &
                                    Control%slope(:,n), pterms,    &
                                    dat, hcew, hcns )
            call diff_mass (Control%Hgrid, dat, hcew, hcns, hkew, hkns, &
                            hdac, Control%order(n))
            hsign = 1.; if (mod(Control%order(n),4) == 0) hsign = -1.
            Var_dt%r(:,:,:,n) = Var_dt%r(:,:,:,n) + hsign * dt_inv * dat
         enddo

      endif

   endif
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------setup momentum diffusion -----------------------------

   if (Control%order(-1) > 0) then

      is = Control%Hgrid%Vel%is-(Control%Hgrid%ihalo-1)
      ie = Control%Hgrid%Vel%ie+(Control%Hgrid%ihalo-1)
      js = Control%Hgrid%Vel%js-(Control%Hgrid%jhalo-1)
      je = Control%Hgrid%Vel%je+(Control%Hgrid%jhalo-1)

      call mass_to_vel (Control%Hgrid, dpde(:,:,nplev+1:nlev), vdat(:,:,nplev+1:nlev))
      call update_halo (Control%Hgrid, UWND, vdat(:,:,nplev+1:nlev), EAST+NORTH+NOPOLE)

      hdac(:,js:je,:) = spread(Control%Hgrid%Vel%rarea(:,js:je),3,nlev)
      do k = nplev+1, nlev
         hdac(:,js:je,k) = hdac(:,js:je,k) / (2.*vdat(:,js:je,k))
      enddo

!     ----- mass-weighted fluxes -----

      hkew3(is:ie+1,js:je+1,1:nplev) = &
                       spread (Control%areavx(is:ie+1,js:je+1),3,nplev)
      hkns3(is:ie+1,js:je+1,1:nplev) = &
                       spread (Control%areavy(is:ie+1,js:je+1),3,nplev)

!dir$ IVDEP
      do k = nplev+1, nlev
      do j = js, je+1
      do i = is, ie+1
         hkew3(i,j,k) = Control%areavx(i,j)*(vdat(i,j,k)+vdat(i-1,j,k))
         hkns3(i,j,k) = Control%areavy(i,j)*(vdat(i,j,k)+vdat(i,j-1,k))
      enddo
      enddo
      enddo

!     ----- mask out fluxes below step-mountain ----

    if (.not.Masks%sigma) then
!dir$ IVDEP
       do k = 1, nlev
       do j = js, je+1
       do i = is, ie+1
          hkew3(i,j,k) = hkew3(i,j,k)*Masks%Vel%mask(i,j,k)*Masks%Vel%mask(i-1,j,k)
          hkns3(i,j,k) = hkns3(i,j,k)*Masks%Vel%mask(i,j,k)*Masks%Vel%mask(i,j-1,k)
       enddo
       enddo
       enddo
    endif

!-----------------------------------------------------------------------
!----- slope adjustment setup ------

    if ( Control % do_slope_adj_wind ) then
!      ---- pressure at velocity points ----
       vdat(:,:,1:nplev) = pres(:,:,1:nplev)
       call mass_to_vel ( Control%Hgrid, pres(:,:,nplev+1:nlev), &
                                         vdat(:,:,nplev+1:nlev)  )
       call update_halo (Control%Hgrid, UWND, vdat(:,:,nplev+1:nlev), &
                                              EAST+NORTH+NOPOLE)

       call vel_slope_correction_init ( Control%Hgrid,       &
                                        Control%nplevs,      &
                                        Control%slope(:,-1), &
                                        vdat, pterms         )
    endif

!-----------------------------------------------------------------------

      hkew2(is:ie+1,js:je+1) = Control%wtv(is:ie+1,js:je+1) * &
                               Control%coeff(-1)
      hkns2(is:ie+1,js:je+1) = Control%wtv(is:ie+1,js:je+1) * &
                               Control%coeff(-1)

!     ---- zero-out cross-polar fluxes ----

      call vel_flux_boundary (Control%Hgrid, hkns2)

!     -- stability condition --
      if ( Control%damping_scheme > 1) then
           hkew2 = min(0.125_r8,hkew2)
           hkns2 = min(0.125_r8,hkns2)
      endif
      hkew = hkew3 * spread(hkew2,3,nlev)
      hkns = hkns3 * spread(hkns2,3,nlev)

!-----------------------------------------------------------------------
!-------------------------momentum diffusion----------------------------

       dat = Var%u + dt*Var_dt%u
      vdat = Var%v + dt*Var_dt%v
      if ( Control % do_slope_adj_wind ) then
         call diff_vel (Control%Hgrid, dat,vdat, hkew,hkns, hdac,  &
                        Control%order(-1), pterms)
      else
         call diff_vel (Control%Hgrid, dat,vdat, hkew,hkns, hdac,  &
                        Control%order(-1))
      endif
      hsign = 1.; if (mod(Control%order(-1),4) == 0) hsign = -1.
      Var_dt%u = Var_dt%u + hsign * dt_inv *  dat
      Var_dt%v = Var_dt%v + hsign * dt_inv * vdat

   endif


!-----------------------------------------------------------------------

 end subroutine horiz_diff

!#######################################################################
!#######################################################################

 subroutine diff_mass (Hgrid, rdat, hcew, hcns, hkew, hkns, hdac, &
                       order)

!----------------------------------------------------------------------
!
!        diff_mass is a private interface that performs multiple
!        2nd order lapacians.
!
!----------------------------------------------------------------------

   type(horiz_grid_type), intent(inout) :: Hgrid
   integer, intent(in)    :: order
   real(r8)   , intent(inout), dimension(Hgrid%ilb:Hgrid%iub, &
                                     Hgrid%jlb:Hgrid%jub, &
                                     nlev) :: rdat
   real(r8)   , intent(in),    dimension(Hgrid%ilb:Hgrid%iub, &
                                     Hgrid%jlb:Hgrid%jub, &
                                     nlev) :: hcew, hcns, &
                                              hkew, hkns, hdac

!----------------------------------------------------------------------
!  hcew, hcns are weighted corrections to the diffusive fluxes for
!  sloping sigma surfaces
!----------------------------------------------------------------------

   real(r8), dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub) :: &
              rew, rns
   integer :: i, j, k, n, is, ie, js, je, nordr
   integer :: is0, ie0, js0, je0
   logical :: do_halos

!-----------------------------------------------------------------------

   is0 = Hgrid%Tmp%is-(Hgrid%ihalo-1);  ie0 = Hgrid%Tmp%ie+(Hgrid%ihalo-1)
   js0 = Hgrid%Tmp%js-(Hgrid%jhalo-1);  je0 = Hgrid%Tmp%je+(Hgrid%jhalo-1)

   is = is0;  ie = ie0;  js = js0;  je = je0

!-----------------------------------------------------------------------
!------------loop for order of diffusion scheme-------------------------

   do n = 1, order/2

        do k = 1, nlev

!------------------contributions (fluxes) ------------------------------

     if ( n == 1 ) then
!dir$ IVDEP
        do j = js-1, je
        do i = is-1, ie
           rew(i,j) = (rdat(i+1,j,k)-rdat(i,j,k)+hcew(i,j,k))*hkew(i,j,k)
           rns(i,j) = (rdat(i,j+1,k)-rdat(i,j,k)+hcns(i,j,k))*hkns(i,j,k)
        enddo
        enddo
     else
!dir$ IVDEP
        do j = js-1, je
        do i = is-1, ie
           rew(i,j) = (rdat(i+1,j,k)-rdat(i,j,k))*hkew(i,j,k)
           rns(i,j) = (rdat(i,j+1,k)-rdat(i,j,k))*hkns(i,j,k)
        enddo
        enddo
     endif

!-----------------------------------------------------------------------

!dir$ IVDEP
     do j = js, je
     do i = is, ie
        rdat(i,j,k)=(rew(i,j)-rew(i-1,j)+rns(i,j)-rns(i,j-1)) &
                    *hdac(i,j,k)
     enddo
     enddo

     enddo

!-----------------------------------------------------------------------
!---- update all halo rows ? ----
!   do not update on last pass, halos will updated in the main program

                                                       do_halos = .false.
     if ( is == Hgrid%Tmp%is .or. js == Hgrid%Tmp%js ) do_halos = .true.
     if ( n == order/2 )                               do_halos = .false.

     if ( do_halos ) then
         call update_halo (Hgrid, TEMP, rdat)
         is = is0;  ie = ie0;  js = js0;  je = je0
     else
         call update_halo (Hgrid, TEMP, rdat, POLEONLY)
         is = is+1;  ie = ie-1
         js = js+1;  je = je-1
     endif

!-----------------------------------------------------------------------

   enddo

!-----------------------------------------------------------------------

 end subroutine diff_mass

!#######################################################################

 subroutine diff_vel (Hgrid, udat, vdat, vkew, vkns, hdac, order, terms)

!----------------------------------------------------------------------
!
!        diff_vel is a private interface that performs multiple
!        2nd order lapacians for the momentum components.
!
!----------------------------------------------------------------------

   type(horiz_grid_type), intent(inout) :: Hgrid
   integer, intent(in)    :: order
   real(r8)   , intent(inout), dimension(Hgrid%ilb:Hgrid%iub, &
                                     Hgrid%jlb:Hgrid%jub, &
                                     nlev) :: udat, vdat
   real(r8)   , intent(in),    dimension(Hgrid%ilb:Hgrid%iub, &
                                     Hgrid%jlb:Hgrid%jub, &
                                     nlev) :: vkew, vkns, hdac
   real(r8)   , intent(in),    dimension(Hgrid%ilb:Hgrid%iub, &
                                     Hgrid%jlb:Hgrid%jub, &
                                     nlev, 3), optional :: terms

!----------------------------------------------------------------------

   real(r8), dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub) ::  &
              uew, uns, vew, vns, dudp, dvdp, ucew, ucns, vcew, vcns
   integer :: i, j, k, n, is, ie, js, je, k1, k2
   integer :: is0, ie0, js0, je0
   logical :: do_halos

!-----------------------------------------------------------------------

   is0 = Hgrid%Vel%is-(Hgrid%ihalo-1);  ie0 = Hgrid%Vel%ie+(Hgrid%ihalo-1)
   js0 = Hgrid%Vel%js-(Hgrid%jhalo-1);  je0 = Hgrid%Vel%je+(Hgrid%jhalo-1)

   is = is0;  ie = ie0;  js = js0;  je = je0

!------------loop for order of diffusion scheme-------------------------

   do n = 1, order/2

      do k = 1, nlev

!------------------contributions (fluxes) ------------------------------

      if ( n == 1 .and. present(terms) ) then

        k1 = max(k-1,1)
        k2 = min(k+1,nlev)
        dudp(:,:) = (udat(:,:,k2)-udat(:,:,k1))*terms(:,:,k,1)
        dvdp(:,:) = (vdat(:,:,k2)-vdat(:,:,k1))*terms(:,:,k,1)

!dir$ IVDEP
        do j = js, je+1
        do i = is, ie+1
!          ---- slope correction terms ----
           ucew(i,j) = (dudp(i,j)+dudp(i-1,j))*terms(i,j,k,2)
           ucns(i,j) = (dudp(i,j)+dudp(i,j-1))*terms(i,j,k,3)
           vcew(i,j) = (dvdp(i,j)+dvdp(i-1,j))*terms(i,j,k,2)
           vcns(i,j) = (dvdp(i,j)+dvdp(i,j-1))*terms(i,j,k,3)

           uew(i,j) = (udat(i,j,k)-udat(i-1,j  ,k)+ucew(i,j))*vkew(i,j,k)
           uns(i,j) = (udat(i,j,k)-udat(i  ,j-1,k)+ucns(i,j))*vkns(i,j,k)
           vew(i,j) = (vdat(i,j,k)-vdat(i-1,j  ,k)+vcew(i,j))*vkew(i,j,k)
           vns(i,j) = (vdat(i,j,k)-vdat(i  ,j-1,k)+vcns(i,j))*vkns(i,j,k)
        enddo
        enddo

      else

!dir$ IVDEP
        do j = js, je+1
        do i = is, ie+1
           uew(i,j) = (udat(i,j,k)-udat(i-1,j  ,k))*vkew(i,j,k)
           uns(i,j) = (udat(i,j,k)-udat(i  ,j-1,k))*vkns(i,j,k)
           vew(i,j) = (vdat(i,j,k)-vdat(i-1,j  ,k))*vkew(i,j,k)
           vns(i,j) = (vdat(i,j,k)-vdat(i  ,j-1,k))*vkns(i,j,k)
        enddo
        enddo

      endif

!-----------------------------------------------------------------------
!dir$ IVDEP
      do j = js, je
      do i = is, ie
        udat(i,j,k) = (uew(i+1,j  )-uew(i,j)+ &
                       uns(i  ,j+1)-uns(i,j))*hdac(i,j,k)
        vdat(i,j,k) = (vew(i+1,j  )-vew(i,j)+ &
                       vns(i  ,j+1)-vns(i,j))*hdac(i,j,k)
      enddo
      enddo

      enddo
!-----------------------------------------------------------------------
!---- update all halo rows ? ----
!   do not update on last pass, halos will updated in the main program

                                                       do_halos = .false.
     if ( is == Hgrid%Vel%is .or. js == Hgrid%Vel%js ) do_halos = .true.
     if ( n == order/2 )                               do_halos = .false.

      if ( do_halos ) then
         call update_halo (Hgrid, UWND, udat)
         call update_halo (Hgrid, VWND, vdat)
         is = is0;  ie = ie0;  js = js0;  je = je0
      else
         call update_halo (Hgrid, UWND, udat, POLEONLY)
         call update_halo (Hgrid, VWND, vdat, POLEONLY)
         is = is+1;  ie = ie-1
         js = js+1;  je = je-1
      endif

!-----------------------------------------------------------------------

   enddo

!-----------------------------------------------------------------------

 end subroutine diff_vel

!#######################################################################

 function horiz_diff_init (Hgrid, nt_adj, nplevs,           &
                           order_vel, order_tmp, order_trs, &
                           coeff_vel, coeff_tmp, coeff_trs, &
                           slope_vel, slope_tmp,            &
                           fis, Vgrid, damping_scheme)      &
                   result (Control)

!-----------------------------------------------------------------------
!       initialization of horizontal diffusion coefficients
!-----------------------------------------------------------------------
!   nplevs = number of "pure" pressure levels at the top of the model

   type(horiz_grid_type), intent(in), target :: Hgrid
   integer, intent(in)               :: nt_adj, nplevs
   integer, intent(in), optional  :: order_vel, order_tmp, order_trs
   real(r8),    intent(in), optional  :: coeff_vel, coeff_tmp, coeff_trs
   real(r8),    intent(in), optional  :: slope_vel(4), slope_tmp(4)

   real(r8),    intent(in), optional  :: fis ( Hgrid%ilb:Hgrid%iub, &
                                           Hgrid%jlb:Hgrid%jub  )
   type(vert_grid_type), intent(in), optional :: Vgrid
   integer,              intent(in), optional :: damping_scheme

type(hdiff_control_type) :: Control

!-----------------------------------------------------------------------

real(r8), dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: dxh2, dxv2
real(r8)    :: eps = 1.e-6
real(r8)    :: dxdy2eq, dx2eq, dyh2, dyv2
integer :: i, j

 integer            :: n, nv, order, ntrace
 real(r8)               :: slope_trs(4)  ! might make this an argument?
 real(r8)               :: coeff, slope(4)
 character(len=128) :: scheme, params, tname

!-----------------------------------------------------------------------

   if (do_log) then
      call write_version_number (version,tagname)
      do_log = .false.
   endif

   Control%Hgrid         => Hgrid
   Control%nplevs        = nplevs
   Control%num_adjust_dt = nt_adj


   call get_number_tracers ( MODEL_ATMOS, num_prog=ntrace )

   allocate ( Control%order   (-1:ntrace), &
              Control%coeff   (-1:ntrace), &
              Control%slope (4,-1:ntrace)  )

 ! defaults
   Control%order = 4
   Control%coeff = .35
   Control%slope =  0.

 ! check optional arguments
   if (present(order_vel)) Control%order  (-1) = order_vel
   if (present(coeff_vel)) Control%coeff  (-1) = coeff_vel
   if (present(slope_vel)) Control%slope(:,-1) = slope_vel

   if (present(order_tmp)) Control%order  (0) = order_tmp
   if (present(coeff_tmp)) Control%coeff  (0) = coeff_tmp
   if (present(slope_tmp)) Control%slope(:,0) = slope_tmp

   if (present(order_trs)) Control%order  (1:ntrace) = order_trs
   if (present(coeff_trs)) Control%coeff  (1:ntrace) = coeff_trs
   ! use default slope values from tmp
   if (present(slope_tmp)) Control%slope(:,1:ntrace) = spread(slope_tmp,2,ntrace)
!!!if (present(slope_trs)) Control%slope(:,1:ntrace) = spread(slope_trs,2,ntrace)


 ! process tracer table information for horizontal diffusion methods
   do n = 1, ntrace
      if (query_method('diff_horiz', MODEL_ATMOS, n, scheme, params)) then
          if (uppercase(trim(scheme)) /= 'LINEAR' .and. uppercase(trim(scheme)) /= 'NONE') &
              call error_mesg  ('bgrid_horiz_diff_mod',  &
                            'invalid diffusion method, '//uppercase(trim(scheme)), FATAL)
          if (parse(params,'order', order) == 1) Control%order(n) = order
          if (uppercase(trim(scheme)) == 'NONE') Control%order(n) = 0
          if (parse(params,'coeff', coeff) == 1) Control%coeff(n) = coeff
          nv = parse(params,'slope', slope)
          Control%slope (1:nv,n) = slope(1:nv)
      endif
         call get_tracer_names (MODEL_ATMOS, n, tname)
         write (stdlog(),10) n, trim(tname), Control%order(n), Control%coeff(n), Control%slope(:,n)
      10 format (i3, a24, ', Order=',i2, ', Coeff=',f10.5, ', Slope=',4f10.5)
   enddo

 ! error checking
   do n = -1, ntrace
     !-- order
      if (Control%order(n) < 0 .or. mod(Control%order(n),2) /= 0) &
          call error_mesg ('bgrid_horiz_diff_mod', 'invalid diffusion order', FATAL)
     !-- non-dimension, normalized coefficient
      if (Control%coeff(n) < 0. .or. Control%coeff(n) > 1.) call error_mesg &
                   ('bgrid_horiz_diff_mod', 'invalid diffusion coeff', FATAL)
     !-- slope correction weights
      if (minval(Control%slope(:,n)) < 0. .or. maxval(Control%slope(:,n)) > 1.) &
            call error_mesg ('bgrid_horiz_diff_mod', 'invalid diffusion coeff', FATAL)
   enddo

 ! set flags
   Control%do_damping        = maxval(Control%order)             > 0
   Control%do_slope_adj_wind = maxval(Control%slope(:,-1))       > 1.e-6
   Control%do_slope_adj_temp = maxval(Control%slope(:,0:ntrace)) > 1.e-6

!-----------------------------------------------------------------------
!----- pre-compute metric terms ------

    Control % damping_scheme = 1
    if (present(damping_scheme)) Control%damping_scheme = damping_scheme

     Control % areahx => var_init(Hgrid)
     Control % areahy => var_init(Hgrid)
     Control % areavx => var_init(Hgrid)
     Control % areavy => var_init(Hgrid)

     Control % wth   => var_init(Hgrid)
     Control % wtv   => var_init(Hgrid)


!    ---- areas averaged along axes ----

     do j = Hgrid%jlb, Hgrid%jub-1
     do i = Hgrid%ilb, Hgrid%iub-1
        Control%areahx(i,j) = 0.5*(Hgrid%Tmp%area(i,j)+Hgrid%Tmp%area(i+1,j))
        Control%areahy(i,j) = 0.5*(Hgrid%Tmp%area(i,j)+Hgrid%Tmp%area(i,j+1))
     enddo
     enddo

     do j = Hgrid%jlb+1, Hgrid%jub
     do i = Hgrid%ilb+1, Hgrid%iub
        Control%areavx(i,j) = 0.5*(Hgrid%Vel%area(i,j)+Hgrid%Vel%area(i-1,j))
        Control%areavy(i,j) = 0.5*(Hgrid%Vel%area(i,j)+Hgrid%Vel%area(i,j-1))
     enddo
     enddo


!    ---- damping weight of x and y axis varies ----
!    ---- depending which damping scheme is used ----
!
!      scheme 1:   equal/constant
!      scheme 2:   function of diagonal grid distance
!      scheme 3:   function of x-axis grid distance
!

  select case (Control%damping_scheme)

!   ---- uniform diffusion ----
    case (1)
       Control%wth = 0.125
       Control%wtv = 0.125

    case (2:3)
       dxh2=0.0; dxv2=0.0
       do j = Hgrid%jlb, Hgrid%jub-1
       do i = Hgrid%ilb, Hgrid%iub-1
          dxh2(i,j) = (0.5*(Hgrid%Tmp%dx(i,j)+Hgrid%Tmp%dx(i+1,j)))**2
       enddo
       enddo
       do j = Hgrid%jlb+1, Hgrid%jub
       do i = Hgrid%ilb+1, Hgrid%iub
          dxv2(i,j) = (0.5*(Hgrid%Vel%dx(i,j)+Hgrid%Vel%dx(i-1,j)))**2
       enddo
       enddo
       dyh2 = Hgrid%Tmp%dy*Hgrid%Tmp%dy
       dyv2 = Hgrid%Vel%dy*Hgrid%Vel%dy


!      ---- function of diagonal distance ----
       if (Control%damping_scheme == 2) then
           dxdy2eq = RADIUS**2*(Hgrid%dlm**2+Hgrid%dph**2)
           Control%wth = 0.125*dxdy2eq/(dxh2+dyh2)
           Control%wtv = 0.125*dxdy2eq/(dxv2+dyv2)
       endif

!      ---- function of x-distance ----
       if (Control%damping_scheme == 3) then
           dx2eq = (RADIUS*Hgrid%dlm)**2
           where (dxh2 /= 0.0) Control%wth = 0.125*dx2eq/dxh2
           where (dxv2 /= 0.0) Control%wtv = 0.125*dx2eq/dxv2
       endif


    case default

       call error_mesg ('bgrid_horiz_diff_mod', 'invalid damping scheme', &
                         FATAL)

  end select


! initialize code sections for performance timing 

!-----------------------------------------------------------------------

 end function horiz_diff_init

!#######################################################################

 subroutine slope_correction_init ( Hgrid, nplevs, pres, terms )

  type(horiz_grid_type), intent(in)  :: Hgrid
  integer,               intent(in)  :: nplevs
  real(r8),                  intent(in)  :: pres (Hgrid%ilb:,Hgrid%jlb:,:)
  real(r8),                  intent(out) :: terms(Hgrid%ilb:,Hgrid%jlb:,:,:)
  integer :: i, j, k, k1, k2, nlev

!  initialization of pressure terms for the sigma slope correction
!  these pressure terms do not change between mass variables
!  may want make weight a function of variable and/or level
!      USE ONE-HALF OF SPECIFIED WEIGHT AT LOWEST LEVEL

   nlev = size(pres,3)

   do k = nplevs+1, nlev
      k1 = max(k-1,1)
      k2 = min(k+1,nlev)
!     --- reciprocal of vert gradient ---
      terms(:,:,k,1) = 1.0/(pres(:,:,k2)-pres(:,:,k1))

!     --- horiz gradients (flip sign) ---
      do j = Hgrid%jlb, Hgrid%jub-1
      do i = Hgrid%ilb, Hgrid%iub-1
           terms(i,j,k,2) = (pres(i,j,k)-pres(i+1,j,k))
           terms(i,j,k,3) = (pres(i,j,k)-pres(i,j+1,k))
      enddo
      enddo
   enddo

 end subroutine slope_correction_init

!#######################################################################

 subroutine vel_slope_correction_init ( Hgrid, nplevs, weights, pres, terms )

  type(horiz_grid_type), intent(in)  :: Hgrid
  integer,               intent(in)  :: nplevs
  real(r8),                  intent(in)  :: weights(4)
  real(r8),                  intent(in)  :: pres (Hgrid%ilb:,Hgrid%jlb:,:)
  real(r8),                  intent(out) :: terms(Hgrid%ilb:,Hgrid%jlb:,:,:)
  real(r8)    :: wt2
  integer :: i, j, k, k1, k2, ks, nlev, isd, ied, jsd, jed
  real(r8) :: dp(size(pres,1),size(pres,2))

!  initialization of pressure terms for the sigma slope correction
!  these pressure terms do not change between mass variables
!  may want make weight a function of variable and/or level
!      USE ONE-HALF OF SPECIFIED WEIGHT AT LOWEST LEVEL

   nlev = size(pres,3)
   isd = Hgrid%Vel%isd; ied = Hgrid%Vel%ied
   jsd = Hgrid%Vel%jsd; jed = Hgrid%Vel%jed

   do k = 1, nplevs
      terms(:,:,k,:) = 0.0
   enddo

   do k = nplevs+1, nlev
      k1 = max(k-1,1)
      k2 = min(k+1,nlev)
!     --- reciprocal of vert gradient ---
      dp = (pres(:,:,k2)-pres(:,:,k1))
      where (dp > 0.)
        terms(:,:,k,1) = 1.0/dp
      elsewhere
        terms(:,:,k,1) = 1.e30 ! these values should not be used where it counts
      endwhere

!     --- horiz gradients (flip sign) ---
      ks = max(1,k-nlev+4)
      wt2  = 0.5*weights(ks)
      do j = Hgrid%jlb+1, Hgrid%jub
      do i = Hgrid%ilb+1, Hgrid%iub
           terms(i,j,k,2) = wt2*(pres(i-1,j,k)-pres(i,j,k))
           terms(i,j,k,3) = wt2*(pres(i,j-1,k)-pres(i,j,k))
      enddo
      enddo
   enddo

 end subroutine vel_slope_correction_init

!#######################################################################

 subroutine slope_correction ( Hgrid, nplevs, weights, terms, temp, cew, cns )

  type(horiz_grid_type), intent(in)  :: Hgrid
  integer,               intent(in)  :: nplevs
  real(r8),                  intent(in)  :: weights(4)
  real(r8),                  intent(in)  :: terms(Hgrid%ilb:,Hgrid%jlb:,:,:)
  real(r8),                  intent(in)  :: temp (Hgrid%ilb:,Hgrid%jlb:,:)
  real(r8),                  intent(out) :: cew  (Hgrid%ilb:,Hgrid%jlb:,:),&
                                        cns  (Hgrid%ilb:,Hgrid%jlb:,:)
  integer :: i, j, k, k1, k2, ks, nlev
  real(r8) :: wt2
  real(r8), dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: dtdp

!  computes weighted-corrections for the slope of sigma surfaces
!  to east-west and north-south fluxes of field temp

 ! check for no correction
   do k = 1, 4
      if (weights(k) > 1.e-6) go to 10
   enddo
   cew = 0.0; cns = 0.0
   return

10 nlev = size(temp,3)

   do k = 1, nplevs
      cew(:,:,k) = 0.0
      cns(:,:,k) = 0.0
   enddo

   do k = nplevs+1, nlev
      k1 = max(k-1,1)
      k2 = min(k+1,nlev)
      dtdp(:,:) = (temp(:,:,k2)-temp(:,:,k1))*terms(:,:,k,1)

      ks = max(1,k-nlev+4)
      wt2  = 0.5*weights(ks)
      do j = Hgrid%jlb, Hgrid%jub-1
      do i = Hgrid%ilb, Hgrid%iub-1
           ! note: terms are grouped for reproducibility with previous version
           cew(i,j,k) = (dtdp(i,j)+dtdp(i+1,j))*(terms(i,j,k,2)*wt2)
           cns(i,j,k) = (dtdp(i,j)+dtdp(i,j+1))*(terms(i,j,k,3)*wt2)
      enddo
      enddo
   enddo

 end subroutine slope_correction

!#######################################################################

end module bgrid_horiz_diff_mod

