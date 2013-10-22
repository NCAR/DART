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

                   module bgrid_advection_mod

!-----------------------------------------------------------------------
!
!            performs vertical and horizontal advection
!
!-----------------------------------------------------------------------
!--------------------------- modules -----------------------------------
!-----------------------------------------------------------------------
use types_mod, only : r8
use bgrid_horiz_mod       , only: horiz_grid_type
use bgrid_vert_mod        , only: vert_grid_type, compute_pres_depth
use bgrid_masks_mod       , only: grid_mask_type
use bgrid_prog_var_mod    , only: prog_var_type
use bgrid_polar_filter_mod, only: pfilt_control_type, polar_filter, &
                                  polar_filter_wind
use bgrid_halo_mod        , only: update_halo, vel_flux_boundary,   &
                                  EAST, WEST, SOUTH, NORTH, NOPOLE, &
                                  TEMP, UWND, VWND, POLEONLY
use bgrid_change_grid_mod , only: mass_to_vel
use vert_advection_mod    , only: vert_advection, &
                                  SECOND_CENTERED, FOURTH_CENTERED, &
                                  FINITE_VOLUME_LINEAR => VAN_LEER_LINEAR
!!!use horiz_advection_mod   , only: horiz_advection

use fms_mod, only: error_mesg, FATAL, write_version_number, &
                   stdlog, uppercase
use  field_manager_mod, only: MODEL_ATMOS, parse
use tracer_manager_mod, only: query_method, get_tracer_names, get_number_tracers
!-----------------------------------------------------------------------

implicit none
private

 public :: advection_init, advection
 public :: advec_control_type

 type advec_control_type
    integer, pointer  :: scheme(:,:)
    real(r8)   , pointer  :: weight(:)
    integer, pointer :: fill_scheme(:), npass_horiz(:), npass_vert(:)
    logical  ::  do_mask4_tmp, do_mask4_vel, do_finite_volume
 end type advec_control_type

 integer, parameter :: NONE = 7000, & 
       FINITE_VOLUME_PARABOLIC = 7004, & 
       EQUAL_BORROW            = 8005, & 
       BOTTOM_BORROW           = 8006, & 
       GLOBAL_BORROW           = 8007  

!-----------------------------------------------------------------------
! private data

 real(r8) :: c4  = -1./6.

 character(len=128) :: version='$Revision$'
 character(len=128) :: tagname='$Id$'

 logical :: stability_check = .false.
 logical :: do_log = .true.

!  timing data

 integer, parameter :: time_level = 5
 integer :: id_advection
 logical :: do_clock_init = .true. 

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine advection ( Control, Pfilt, Hgrid, Vgrid, Masks,   &
                       pfilt_opt, dt, dpde_old, dpde,         &
                       few, fns, etadot, u, v, t, Var, Var_dt )

type(advec_control_type), intent(in) :: Control
type(pfilt_control_type), intent(in) :: Pfilt
type(horiz_grid_type), intent(inout) :: Hgrid
type (vert_grid_type), intent(in)    :: Vgrid
type (grid_mask_type), intent(in)    :: Masks
              integer, intent(in)    :: pfilt_opt
                 real(r8), intent(in)    :: dt
real(r8), intent(in), dimension(Hgrid%ilb:,Hgrid%jlb:,:) ::  &
                                            dpde_old, dpde, u, v, t
real(r8), intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: &
                                            few, fns, etadot
type (prog_var_type), intent(in)    :: Var
type (prog_var_type), intent(inout) :: Var_dt

!-----------------------------------------------------------------------
!
!  Var       prognostic variables at the end of the last dynamics time
!            step, current values would be, uc = Var%u + dt*Var_dt%u
!
!  Var_dt    prognostic variable tendencies, accumulated since the
!            variable were updated in Var
!
!  u, v, t   prognostic variables at the end of the last advective time
!            step, note that tracers have not been updated since the last
!            advective time step, therefore, r = Var%r + dt*Var_dt%r
!
!-----------------------------------------------------------------------
  real(r8), dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub) ::  &
                                 few_2d, fns_2d
  real(r8), dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub,  &
                  Vgrid%nlev) :: dpdt, dpde_xy, r, uc, vc
  real(r8), dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub,  &
                  Vgrid%nlev,2) :: mask4
  integer :: i, j, k, n, is, ie, js, je, kx, order
!-----------------------------------------------------------------------


!  ---- update halos for fluxes ----
   call update_halo (Hgrid, TEMP, few)
   call update_halo (Hgrid, VWND, fns)

!  ---- compute pressure tendency in each layer ----

   dpdt = (dpde - dpde_old) / dt

!---- stability diagnostic (aka CFL for finite volume scheme) ----
   if (stability_check) call stability_diag ( Hgrid, dt, few, dpde_old, dpde )

!-----------------------------------------------------------------------
!------------------ advection of mass fields ---------------------------
!-----------------------------------------------------------------------
!------ initialize fourth-order mask ----

   if (Control%do_mask4_tmp) &
         call mask4_init (Hgrid, 1, Masks%sigma, Masks%Tmp%mask, mask4)

 ! compute wind components for finite-volume advection scheme
   if (Control%do_finite_volume) then
        call compute_advecting_wind (Hgrid, dpde_old, dpde, few, fns, uc, vc)
   else
        uc = 0.
        vc = 0.
   endif
!------ temperature advection ------
   call advect_mass (Hgrid, Pfilt, Control % scheme(:,0),   &
                     Control % weight(0), pfilt_opt,        &
                     dt, dpdt, dpde, Masks%Tmp%mask, mask4, &
                     few, fns, etadot, uc, vc, Var%t, t, Var_dt%t   )
!------ tracer advection -----------
!  (need to pass all tracers to update pressure tendency part)

   do n = 1, Var_dt%ntrace

      r = Var%r(:,:,:,n) + dt*Var_dt%r(:,:,:,n)

      call advect_mass (Hgrid, Pfilt, Control%scheme(:,n), &
                        Control % weight(n),               &
                        pfilt_opt, dt, dpdt, dpde,         &
                        Masks%Tmp%mask, mask4, few, fns, etadot, uc, vc, &
                        Var%r(:,:,:,n), r, Var_dt%r(:,:,:,n)     )

   enddo
!------ tracer hole filling -------
!  (remove negative tracer values with vert/horiz borrowing)
     if(Var_dt%ntrace /= 0) then
        call vert_borrow ( dt, dpde, Var%r(:,:,:,1:Var_dt%ntrace), &
                               Var_dt%r(:,:,:,1:Var_dt%ntrace), &
                                  iters = Control % npass_vert  )

        call horiz_borrow ( Hgrid, dt, dpde, Masks%Tmp%mask,    &
                               Var%r(:,:,:,1:Var_dt%ntrace), &
                            Var_dt%r(:,:,:,1:Var_dt%ntrace), &
                              iters = Control % npass_horiz  )
      endif

!-----------------------------------------------------------------------
!------------------ advection of momentum fields -----------------------
!-----------------------------------------------------------------------
!------------- mass fluxes between velocity points ---------------------

   is = Hgrid % Vel % is;  ie = Hgrid % Vel % ie
   js = Hgrid % Vel % js;  je = Hgrid % Vel % je
   kx = size(few,3)

! Need to set unused values of few_2d to zero
few_2d = 0.0; fns_2d = 0.0

   do k = 1, kx
   do j = js, je
   do i = is, ie+1
      few_2d(i,j) = (few(i,j  ,k)+few(i-1,j  ,k)+  &
                     few(i,j+1,k)+few(i-1,j+1,k))*0.25
   enddo
   enddo
   do j = js  , je+1
   do i = is, ie
      fns_2d(i,j) = (fns(i  ,j,k)+fns(i  ,j-1,k)+  &
                     fns(i+1,j,k)+fns(i+1,j-1,k))*0.25
   enddo
   enddo
!     ---- eliminate meridional fluxes at sub pole row ----
!     ---- for mass conservation ----
      call vel_flux_boundary (Hgrid, fns_2d)

!     ---- store back into 3d arrays ----
      few(:,:,k) = few_2d
      fns(:,:,k) = fns_2d
   enddo

!  ---- interpolate mass fields to momentum grid ----

   call mass_to_vel (Hgrid,   dpde, dpde_xy)
   call mass_to_vel (Hgrid, etadot,  etadot)
   call mass_to_vel (Hgrid,   dpdt,    dpdt)

   if (Control%do_mask4_vel) then
       call mask4_init (Hgrid, 2, Masks%sigma, Masks%Vel%mask, mask4)
   endif

   call advect_vel (Hgrid, Pfilt, Control % scheme(:,-1),    &
                    Control % weight(-1),                    &
                    pfilt_opt, dt, dpdt, dpde_xy,            &
                    Masks%Vel%mask, mask4, few, fns, etadot, &
                    Var%u, Var%v, u, v, Var_dt%u, Var_dt%v   )

!-----------------------------------------------------------------------
!------ done with fluxes - zero-out ??? -------

   few    = 0.0
   fns    = 0.0
   etadot = 0.0

!-----------------------------------------------------------------------


!-----------------------------------------------------------------------

end subroutine advection

!#######################################################################

subroutine advect_mass ( Hgrid, Pfilt, scheme, coeff, fopt, dt, &
                         dpdt, dpn, mask, mask4, few, fns, fud, &
                         uc, vc, ro, r, r_dt )

  type(horiz_grid_type),    intent(inout) :: Hgrid
  type(pfilt_control_type), intent(in) :: Pfilt
  integer,                  intent(in) :: scheme(2)
  integer,                  intent(in) :: fopt
  real(r8), intent(in)                     :: coeff, dt
  real(r8), intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:)   :: dpdt, dpn,&
                                     mask, few, fns, fud, ro, r, uc, vc
  real(r8), intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:,:) :: mask4
  real(r8), intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:)   :: r_dt

!-----------------------------------------------------------------------

  real(r8), dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub,  &
                  size(r,3)) :: rst, rst_dt, rst_dt_v, rdpdt
  real(r8), dimension(Hgrid%jlb:Hgrid%jub) :: dx, dy

  real(r8)    :: rcoef
  integer :: pass, npass, j, k, order
  integer :: vert_scheme, horiz_scheme
  integer, parameter :: FINITE_VOLUME = 1234567 !can not match other numbers?

!-----------------------------------------------------------------------
 ! set the vertical and horizontal differencing scheme
    horiz_scheme = scheme(1)
     vert_scheme = scheme(2)
    if (horiz_scheme == SECOND_CENTERED) order = 2
    if (horiz_scheme == FOURTH_CENTERED) order = 4

!-----------------------------------------------------------------------

   rst = r
   rdpdt = r*dpdt
   rcoef = coeff

                      npass = 1
   if (rcoef > 1.e-4) npass = 2
 ! no advection, psdt correction only
   if (vert_scheme ==  NONE .and. horiz_scheme == NONE) then 
       npass = 1
       rst_dt_v = 0.0
   endif
 ! both vert and horiz are using finite volume
   if ( vert_scheme == FINITE_VOLUME_LINEAR .and. &
       horiz_scheme == FINITE_VOLUME_LINEAR  ) rcoef = 1.

   do pass = 1, npass

     !---- vertical tendency ----
      if (vert_scheme == SECOND_CENTERED .or. &
          vert_scheme == FOURTH_CENTERED .or. &
         (vert_scheme == FINITE_VOLUME_LINEAR .and. pass == 1)) then
! need to add parabolic check
         call vert_advection ( dt, fud, dpn, rst, rst_dt_v,  &
                               mask=mask, scheme=vert_scheme )
      endif

     !---- horizontal tendency ----
      if (horiz_scheme == SECOND_CENTERED .or. &
          horiz_scheme == FOURTH_CENTERED) then
         call advect_mass_horiz (Hgrid, order, few, fns,    &
                                 mask4(:,:,:,1), mask4(:,:,:,2), &
                                 rst, rst_dt                     )
       ! polar filter horizontal tendency
         if ( fopt >= 1 ) call polar_filter (Pfilt, rst_dt, 1, mask)
         call update_halo (Hgrid, TEMP, rst_dt)

    ! compute horizontal finite volume scheme on 2nd pass
      else if ((horiz_scheme == FINITE_VOLUME_LINEAR .or. &
                horiz_scheme == FINITE_VOLUME_PARABOLIC) .and. pass == 2) then
        ! temporary error check
         call error_mesg ('bgrid_advection_mod', &
                          'horizontal finite volume schemes &
                          &are not implemented with this release', FATAL)
         dx = Hgrid%Tmp%dx(Hgrid%Tmp%is,:)
         dy = Hgrid%Tmp%dy
        !call horiz_advection ( Hgrid%Tmp%Domain, dt, dx, dy, &
        !                       uc, vc, few, fns, rst, rst_dt )
         do k = 1, size(r,3)
           rst_dt(:,:,k) = rst_dt(:,:,k)*Hgrid%Tmp%rarea
         enddo
       ! halos updated by finite volume scheme
       ! call update_halo (Hgrid, TEMP, rst_dt)

      else
            rst_dt = 0.0
      endif
      

     !---- combine vertical and horizontal tendencies ----
     !---- adjust for pressure tendency ----
        rst_dt = (rst_dt + rst_dt_v - rdpdt) / dpn * mask

     !---- compute new value or return tendency ----
        if (pass < npass) then
            rst = ro + dt * (r_dt + rcoef*rst_dt)
        else
            r_dt = r_dt + rst_dt
        endif

   enddo

!-----------------------------------------------------------------------

end subroutine advect_mass

!#######################################################################

subroutine advect_vel ( Hgrid, Pfilt, scheme, coeff, fopt, dt,  &
                        dpdt, dpn, mask, mask4, few, fns, fud, &
                        uo, vo, u, v, u_dt, v_dt               )

  type(horiz_grid_type),    intent(inout) :: Hgrid
  type(pfilt_control_type), intent(in) :: Pfilt
  integer,                  intent(in) :: scheme(2)
  integer,                  intent(in) :: fopt
  real(r8), intent(in)                     :: coeff, dt
  real(r8), intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:)   :: dpdt, dpn,&
                                     mask, few, fns, fud, uo, vo, u, v
  real(r8), intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:,:) :: mask4
  real(r8), intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: u_dt, v_dt

!-----------------------------------------------------------------------

  real(r8), dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub,  &
                  size(u,3)) :: ust, ust_dt, ust_dt_v, &
                                vst, vst_dt, vst_dt_v
  integer :: order, vert_scheme, horiz_scheme
  integer :: pass, npass, i, j, k, is, ie, js, je, nlev

!-----------------------------------------------------------------------

 ! set the vertical differencing scheme
   horiz_scheme = scheme(2)
    vert_scheme = scheme(2)

    if (horiz_scheme == SECOND_CENTERED) order = 2
    if (horiz_scheme == FOURTH_CENTERED) order = 4

!-----------------------------------------------------------------------

   is = Hgrid%Vel%is; ie = Hgrid%Vel%ie
   js = Hgrid%Vel%js; je = Hgrid%Vel%je
   nlev = size(u,3)

   ust = u
   vst = v

                      npass = 1
   if (coeff > 1.e-4) npass = 2

   do pass = 1, npass

     !---- vertical tendency ----
      if (vert_scheme == SECOND_CENTERED .or. &
          vert_scheme == FOURTH_CENTERED .or. &
         (vert_scheme == FINITE_VOLUME_LINEAR .and. pass == 1)) then
         call vert_advection ( dt, fud, dpn, ust, ust_dt_v,  &
                               mask=mask, scheme=vert_scheme )
         call vert_advection ( dt, fud, dpn, vst, vst_dt_v,  &
                               mask=mask, scheme=vert_scheme )
      endif

     !---- horizontal tendency ----
        call advect_vel_horiz (Hgrid, order, few, fns,         &
                               mask4(:,:,:,1), mask4(:,:,:,2), &
                               ust, vst, ust_dt, vst_dt        )
     ! polar filter horizontal tendency
        if (fopt == 2) call polar_filter_wind (Pfilt, ust_dt, vst_dt, mask)

     !---- combine vertical and horizontal tendencies ----
     !---- adjust for pressure tendency ----
      do k = 1, nlev
      do j = js, je
      do i = is, ie
         ust_dt(i,j,k) = (ust_dt(i,j,k) + ust_dt_v(i,j,k) - &
                       u(i,j,k)*dpdt(i,j,k)) / dpn(i,j,k) * mask(i,j,k)
         vst_dt(i,j,k) = (vst_dt(i,j,k) + vst_dt_v(i,j,k) - &
                       v(i,j,k)*dpdt(i,j,k)) / dpn(i,j,k) * mask(i,j,k)
      enddo
      enddo
      enddo

     !---- compute new value or return tendency ----
      if (pass < npass) then
         do k = 1, nlev
         do j = js, je
         do i = is, ie
           ust(i,j,k) = uo(i,j,k) + dt * (u_dt(i,j,k) + coeff*ust_dt(i,j,k))
           vst(i,j,k) = vo(i,j,k) + dt * (v_dt(i,j,k) + coeff*vst_dt(i,j,k))
         enddo
         enddo
         enddo
         call update_halo (Hgrid, UWND, ust)
         call update_halo (Hgrid, VWND, vst)
      else
         do k = 1, nlev
         do j = js, je
         do i = is, ie
           u_dt(i,j,k) = u_dt(i,j,k) + ust_dt(i,j,k)
           v_dt(i,j,k) = v_dt(i,j,k) + vst_dt(i,j,k)
         enddo
         enddo
         enddo
      endif
   
   enddo

!-----------------------------------------------------------------------

end subroutine advect_vel

!#######################################################################

 subroutine advect_mass_horiz (Hgrid, order, few, fns, maskx, masky, &
                               rst, rst_dt)

  type(horiz_grid_type),    intent(inout) :: Hgrid
  integer,                  intent(in)    :: order
  real(r8), intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: few, fns, &
                                                       maskx, masky, rst
  real(r8), intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: rst_dt

  real(r8), dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub, &
                  size(rst,3)) :: rstx, rsty, rew, rns
  real(r8), dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub) :: &
                               x2, x4, y2, y4
  integer :: i, j, k, is, ie, js, je, kx

   if (order == 0) then
      rst_dt = 0.0
      return
   endif

   is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie
   js = Hgrid % Tmp % js;  je = Hgrid % Tmp % je
   kx = size(rst,3)

   if (order == 4) then
      rsty(:,Hgrid%jub,:) = 0.0
      do k = 1, kx
        do j = Hgrid%jlb,Hgrid%jub
        do i = Hgrid%ilb,Hgrid%iub-1
           rstx(i,j,k) = rst(i+1,j,k)+rst(i,j,k)
        enddo
        enddo
        do j = Hgrid%jlb,Hgrid%jub-1
        do i = Hgrid%ilb,Hgrid%iub
           rsty(i,j,k) = rst(i,j+1,k)+rst(i,j,k)
        enddo
        enddo
      enddo
      if (Hgrid%ihalo == 1) call update_halo (Hgrid,TEMP,rstx,flags=EAST)
      if (Hgrid%jhalo == 1) call update_halo (Hgrid,UWND,rsty,flags=NORTH+NOPOLE)
   endif


!     ---- horizontal mass fluxes ----

      if ( order == 4 ) then
         rew=0.0; rns=0.0
         do k = 1, kx
           do j = js,          je
           do i = Hgrid%ilb+1, ie
             x4(i,j) = maskx(i,j,k)
             x2(i,j) = 1.0 - 2.*x4(i,j)
             rew(i,j,k) = few(i,j,k) * ( x2(i,j) * rstx(i,j,k) &
                         + x4(i,j) * (rstx(i-1,j,k) + rstx(i+1,j,k)) )
           enddo
           enddo
           do j = Hgrid%jlb+1, je
           do i = is,          ie
             y4(i,j) = masky(i,j,k)
             y2(i,j) = 1.0 - 2.*y4(i,j)
             rns(i,j,k) = fns(i,j,k) * ( y2(i,j) * rsty(i,j,k) &
                         + y4(i,j) * (rsty(i,j-1,k) + rsty(i,j+1,k)) )
           enddo
           enddo
         enddo
         if (Hgrid%ihalo == 1) call update_halo (Hgrid,TEMP,rew,flags=WEST)
         if (Hgrid%jhalo == 1) call update_halo (Hgrid,VWND,rns,flags=SOUTH)
      else if ( order == 2 ) then
         do k = 1, kx
           do j = js-1, je
           do i = is-1, ie
             rew(i,j,k) = few(i,j,k) * (rst(i+1,j,k) + rst(i,j,k))
             rns(i,j,k) = fns(i,j,k) * (rst(i,j+1,k) + rst(i,j,k))
           enddo
           enddo
         enddo
      else
         ! actually invalid
           rew = 0.0
           rns = 0.0
      endif

!     ---- horizontal advective tendency ----

      do k = 1, kx
      do j = js, je
      do i = is, ie
         rst_dt(i,j,k) = -0.5*Hgrid%Tmp%rarea(i,j)*            &
                                 ((rew(i  ,j,k)+rns(i,j  ,k))  &
                                 -(rew(i-1,j,k)+rns(i,j-1,k)))
      enddo
      enddo
      enddo


 end subroutine advect_mass_horiz

!#######################################################################

 subroutine advect_vel_horiz (Hgrid, order, few, fns, maskx, masky, &
                              ust, vst, ust_dt, vst_dt)

  type(horiz_grid_type),    intent(inout) :: Hgrid
  integer,                  intent(in)    :: order
  real(r8), intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: few, fns, &
                                                 maskx, masky, ust, vst
  real(r8), intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: ust_dt, vst_dt

  real(r8), dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub, &
                  size(ust,3)) :: ustx, usty, uew, uns, &
                                  vstx, vsty, vew, vns
  real(r8), dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub) :: &
                               x2, x4, y2, y4
  integer :: i, j, k, is, ie, js, je, kx

   if (order == 0) then
      ust_dt = 0.0
      vst_dt = 0.0
      return
   endif

   is = Hgrid % Vel % is;  ie = Hgrid % Vel % ie
   js = Hgrid % Vel % js;  je = Hgrid % Vel % je
   kx = size(ust,3)

   if (order == 4) then
      do k = 1, kx
        do j = Hgrid%jlb  ,Hgrid%jub
        do i = Hgrid%ilb+1,Hgrid%iub
           ustx(i,j,k) = ust(i-1,j,k)+ust(i,j,k)
           vstx(i,j,k) = vst(i-1,j,k)+vst(i,j,k)
        enddo
        enddo
        do j = Hgrid%jlb+1,Hgrid%jub
        do i = Hgrid%ilb  ,Hgrid%iub
           usty(i,j,k) = ust(i,j-1,k)+ust(i,j,k)
           vsty(i,j,k) = vst(i,j-1,k)+vst(i,j,k)
        enddo
        enddo
      enddo
      if (Hgrid%ihalo == 1) call update_halo (Hgrid,UWND,ustx,flags=WEST)
      if (Hgrid%ihalo == 1) call update_halo (Hgrid,VWND,vstx,flags=WEST)
      if (Hgrid%jhalo == 1) call update_halo (Hgrid,TEMP,usty,flags=SOUTH)
      if (Hgrid%jhalo == 1) call update_halo (Hgrid,TEMP,vsty,flags=SOUTH)
   endif


!     ---- horizontal mass fluxes ----

      if ( order == 4 ) then
         do k = 1, kx
           do j = js, je
           do i = is, Hgrid%iub-1
             x4(i,j) = maskx(i,j,k)
             x2(i,j) = 1.0 - 2.*x4(i,j)
             uew(i,j,k) = few(i,j,k) * ( x2(i,j) * ustx(i,j,k) &
                         + x4(i,j) * (ustx(i-1,j,k) + ustx(i+1,j,k)) )
             vew(i,j,k) = few(i,j,k) * ( x2(i,j) * vstx(i,j,k) &
                         + x4(i,j) * (vstx(i-1,j,k) + vstx(i+1,j,k)) )
           enddo
           enddo
           do j = js, Hgrid%jub-1
           do i = is, ie
             y4(i,j) = masky(i,j,k)
             y2(i,j) = 1.0 - 2.*y4(i,j)
             uns(i,j,k) = fns(i,j,k) * ( y2(i,j) * usty(i,j,k) &
                         + y4(i,j) * (usty(i,j-1,k) + usty(i,j+1,k)) )
             vns(i,j,k) = fns(i,j,k) * ( y2(i,j) * vsty(i,j,k) &
                         + y4(i,j) * (vsty(i,j-1,k) + vsty(i,j+1,k)) )
           enddo
           enddo
         enddo
         if (Hgrid%ihalo == 1) call update_halo (Hgrid,UWND,uew,flags=EAST)
         if (Hgrid%ihalo == 1) call update_halo (Hgrid,VWND,vew,flags=EAST)
         if (Hgrid%jhalo == 1) call update_halo (Hgrid,TEMP,uns,flags=NORTH)
         if (Hgrid%jhalo == 1) call update_halo (Hgrid,TEMP,vns,flags=NORTH)
      else if ( order == 2 ) then
         do k = 1, kx
           do j = js, je
           do i = is, ie+1
             uew(i,j,k) = few(i,j,k) * (ust(i-1,j,k) + ust(i,j,k))
             vew(i,j,k) = few(i,j,k) * (vst(i-1,j,k) + vst(i,j,k))
           enddo
           enddo
           do j = js, je+1
           do i = is, ie
             uns(i,j,k) = fns(i,j,k) * (ust(i,j-1,k) + ust(i,j,k))
             vns(i,j,k) = fns(i,j,k) * (vst(i,j-1,k) + vst(i,j,k))
           enddo
           enddo
         enddo
      else
         ! actually invalid
           uew = 0.0
           uns = 0.0
      endif

!     ---- horizontal advective tendency ----

      do k = 1, kx
      do j = js, je
      do i = is, ie
         ust_dt(i,j,k) = -0.5*Hgrid%Vel%rarea(i,j)*            &
                                 ((uew(i+1,j,k)+uns(i,j+1,k))  &
                                 -(uew(i  ,j,k)+uns(i,j  ,k)))
         vst_dt(i,j,k) = -0.5*Hgrid%Vel%rarea(i,j)*            &
                                 ((vew(i+1,j,k)+vns(i,j+1,k))  &
                                 -(vew(i  ,j,k)+vns(i,j  ,k)))
      enddo
      enddo
      enddo


 end subroutine advect_vel_horiz

!#######################################################################
! computes advecting wind for finite-volume advection

 subroutine compute_advecting_wind ( Hgrid, dpo, dpn, few, fns, uc, vc )
 type(horiz_grid_type), intent(inout) :: Hgrid
 real(r8), intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: dpo, dpn, few, fns
 real(r8), intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: uc, vc

 real(r8) :: dp
 integer :: i, j, k

   do k = 1, size(uc,3)
   do j = Hgrid%jlb, Hgrid%jub
   do i = Hgrid%ilb, Hgrid%iub-1
     dp = (dpo(i,j,k)+dpn(i,j,k)+dpo(i+1,j,k)+dpn(i+1,j,k))*0.25
     uc(i,j,k) = few(i,j,k)/(dp*Hgrid%Tmp%dy)
   enddo
   enddo
   enddo

   do k = 1, size(uc,3)
   do j = Hgrid%jlb, Hgrid%jub-1
   do i = Hgrid%ilb, Hgrid%iub
     dp = (dpo(i,j,k)+dpn(i,j,k)+dpo(i,j+1,k)+dpn(i,j+1,k))*0.25
     vc(i,j,k) = fns(i,j,k)/(dp*Hgrid%Vel%dx(i,j))
   enddo
   enddo
   enddo

   call update_halo ( Hgrid, TEMP, uc, flags=EAST )
   call update_halo ( Hgrid, VWND, vc, flags=NORTH )

 end subroutine compute_advecting_wind

!#######################################################################
!---- stability diagnostic (aka CFL for finite volume scheme) ----

 subroutine stability_diag ( Hgrid, dt, few, dpo, dpn )
 type(horiz_grid_type), intent(inout) :: Hgrid
 real(r8), intent(in)                     :: dt
 real(r8), intent(in), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: few, dpo, dpn

 real(r8), dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: uc
 real(r8)    :: cflmax, cfl, dpa
 integer :: i, j, k

    cflmax = 0.
    do k = 1, size(few,3)
     ! compute zonal wind
       do j = Hgrid%Tmp%js  , Hgrid%Tmp%je
       do i = Hgrid%Tmp%is-1, Hgrid%Tmp%ie
         dpa = (dpo(i,j,k)+dpn(i,j,k)+dpo(i+1,j,k)+dpn(i+1,j,k))*0.25
         uc(i,j) = few(i,j,k)/(dpa*Hgrid%Tmp%dy)
       enddo     
       enddo   
       do j = Hgrid%Tmp%js, Hgrid%Tmp%je
       do i = Hgrid%Tmp%is, Hgrid%Tmp%ie
          cfl = abs(uc(i,j)-uc(i-1,j))*dt/Hgrid%Tmp%dx(i,j)
          cflmax = max(cflmax,cfl)
       enddo      
       enddo   
    enddo   

! Running on single PE, don't need to find maxes from other pe's
        if (cflmax > 1.0) then
            print *, 'x-axis stability violated, cfl = ', cflmax
        endif       

 end subroutine stability_diag

!#######################################################################

 subroutine mask4_init (Hgrid, grid, sigma, mask, mask4)

 type(horiz_grid_type), intent(inout) :: Hgrid
 integer, intent(in)                                      :: grid
 logical, intent(in)                                      :: sigma
 real(r8), intent(in),   dimension(Hgrid%ilb:,Hgrid%jlb:,:)   :: mask
 real(r8), intent(out),  dimension(Hgrid%ilb:,Hgrid%jlb:,:,:) :: mask4
 real(r8), dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub, &
                 size(mask,3)) :: mew, mns

 integer :: i,j,k,is,ie,js,je,n


! horizontal masks

  select case (grid)

  case (1)
     if (sigma) then
        mask4(:,:,:,1:2) = c4
     else

      ! initialize
        mask4 = 0.0
        mns(:,Hgrid%jub,:) = 0.0

        do k = 1, size(mask,3)
          do j = Hgrid%jlb,Hgrid%jub
          do i = Hgrid%ilb,Hgrid%iub-1
             mew(i,j,k) = mask(i+1,j,k)*mask(i,j,k)
          enddo
          enddo
          do j = Hgrid%jlb,Hgrid%jub-1
          do i = Hgrid%ilb,Hgrid%iub
             mns(i,j,k) = mask(i,j+1,k)*mask(i,j,k)
          enddo
          enddo
        enddo
        if (Hgrid%ihalo == 1) call update_halo (Hgrid,TEMP,mew,flags=EAST)
        if (Hgrid%jhalo == 1) call update_halo (Hgrid,UWND,mns,flags=NORTH+NOPOLE)

         do k = 1, size(mask,3)
           do j = Hgrid%Tmp%js, Hgrid%Tmp%je
           do i = Hgrid%ilb+1,  Hgrid%Tmp%ie
              mask4(i,j,k,1) = c4 * mew(i-1,j,k) * mew(i+1,j,k)
           enddo
           enddo
           do j = Hgrid%jlb+1,  Hgrid%Tmp%je
           do i = Hgrid%Tmp%is, Hgrid%Tmp%ie
              mask4(i,j,k,2) = c4 * mns(i,j-1,k) * mns(i,j+1,k)
           enddo
           enddo
         enddo
     endif

  case (2)
     if (sigma) then
        mask4(:,:,:,1:2) = c4
        mns = 1.0
     else
      ! initialize
        mask4 = 0.0
        do k = 1, size(mask,3)
          do j = Hgrid%jlb  ,Hgrid%jub
          do i = Hgrid%ilb+1,Hgrid%iub
             mew(i,j,k) = mask(i-1,j,k)*mask(i,j,k)
          enddo
          enddo
          do j = Hgrid%jlb+1,Hgrid%jub
          do i = Hgrid%ilb  ,Hgrid%iub
             mns(i,j,k) = mask(i,j-1,k)*mask(i,j,k)
          enddo
          enddo
        enddo
        if (Hgrid%ihalo == 1) call update_halo (Hgrid,UWND,mew,flags=WEST)
        if (Hgrid%jhalo == 1) call update_halo (Hgrid,TEMP,mns,flags=SOUTH+NOPOLE)
      endif
    ! use second order in meridional direction near poles
      call vel_flux_boundary (Hgrid,mns)

      if (.not.sigma) then
        do k = 1, size(mask,3)
          do j = Hgrid%Vel%js, Hgrid%Vel%je
          do i = Hgrid%Vel%is, Hgrid%iub-1
             mask4(i,j,k,1) = c4 * mew(i-1,j,k) * mew(i+1,j,k)
          enddo
          enddo
        enddo
      endif

      do k = 1, size(mask,3)
        do j = Hgrid%Vel%js, Hgrid%jub-1
        do i = Hgrid%Vel%is, Hgrid%Vel%ie
           mask4(i,j,k,2) = c4 * mns(i,j-1,k) * mns(i,j+1,k)
        enddo
        enddo
      enddo

  end select


 end subroutine mask4_init

!#######################################################################

  function advection_init (Hgrid,  &
                            order_vel,  order_tmp,  order_trs, &
                           weight_vel, weight_tmp, weight_trs, &
                           horiz_fill, vert_fill) result (Control)

  type(horiz_grid_type), intent(in) :: Hgrid
  integer, intent(in), optional     ::  order_vel,  order_tmp,  &
                                        order_trs
  real(r8),    intent(in), optional     :: weight_vel, weight_tmp,  &
                                       weight_trs
  integer, intent(in), optional     :: horiz_fill, vert_fill
  type(advec_control_type)          :: Control

 real(r8)    :: wt
 integer :: n, m, np, ntrace
 character(len=32)  :: advec_methods(2) = (/ 'advec_horiz', 'advec_vert ' /)
 character(len=128) :: scheme, params, name
 integer            ::  advec_default(2)
!-----------------------------------------------------------------------

   if (do_log) then
      call write_version_number (version,tagname)
      do_log = .false.
   endif

   call get_number_tracers ( MODEL_ATMOS, num_prog=ntrace )

   allocate ( Control%scheme (2,-1:ntrace), &
              Control%weight   (-1:ntrace)  )
   allocate ( Control%fill_scheme (ntrace), &
              Control%npass_horiz (ntrace), &
              Control%npass_vert  (ntrace)  )

 ! defaults
   Control%scheme = SECOND_CENTERED
   Control%weight = 0.70
   Control%fill_scheme = NONE
   Control%npass_horiz = 0
   Control%npass_vert  = 0

 ! check optional arguments
   if (present( order_vel)) Control%scheme(:,-1) = set_default_advec (order_vel)
   if (present(weight_vel)) Control%weight  (-1) = weight_vel

   if (present( order_tmp)) Control%scheme(:,0) = set_default_advec (order_tmp)
   if (present(weight_tmp)) Control%weight  (0) = weight_tmp

   if (present( order_trs)) then
                    advec_default = set_default_advec (order_trs)
                    Control%scheme(:,1:ntrace) = spread(advec_default,2,ntrace)
   endif
   if (present(weight_trs)) Control%weight  (1:ntrace) = weight_trs

 ! optional arguments for filling
   if (present(horiz_fill)) Control%npass_horiz = max(0,horiz_fill)
   if (present( vert_fill)) Control%npass_vert  = max(0, vert_fill)
   if (maxval(Control%npass_horiz) > 0 .or. maxval(Control%npass_vert) > 0) &
                            Control%fill_scheme = EQUAL_BORROW

 ! process tracer table information for advection methods
   do n = 1, ntrace
      do m = 1, 2
         if (query_method(trim(advec_methods(m)), MODEL_ATMOS, n, scheme, params)) then
             Control%scheme(m,n) = set_advec_scheme( scheme )
             ! parse e-b weight
             if (Control%scheme(m,n) == SECOND_CENTERED .or. &
                 Control%scheme(m,n) == FOURTH_CENTERED) then
                   if (parse(params,'wt', wt) == 1) Control%weight(n) = wt
             endif
         endif
      enddo
   enddo

 ! error check on Euler weight
   do n = -1, ntrace
        if (Control%weight(n) < 0.0 .or. Control%weight(n) > 1.0) &
                          call error_mesg ('bgrid_advection_mod', &
                            'E-B weight out of range [0,1]', FATAL)
   enddo

 ! process tracer table information for filling method
   do n = 1, ntrace
      if (query_method('filling', MODEL_ATMOS, n, scheme, params)) then
          if (uppercase(trim(scheme)) == 'LOCAL') then
              Control%fill_scheme(n) = EQUAL_BORROW
              Control%npass_horiz(n) = 1
              Control%npass_vert (n) = 1
         !else if (uppercase(trim(scheme)) == 'GLOBAL') then
         !    Control%fill_scheme(n) = GLOBAL_BORROW
          else if (uppercase(trim(scheme)) == 'NONE') then
              Control%fill_scheme(n) = NONE
          else
              call error_mesg ('bgrid_advection_mod',  &
                   'invalid filling scheme, '//uppercase(trim(scheme)), FATAL)
          endif
         ! parse number of filling passes
           if (Control%fill_scheme(n) == EQUAL_BORROW) then
              if (parse(params,'hp', np) == 1) Control%npass_horiz(n) = np
              if (parse(params,'vp', np) == 1) Control%npass_vert (n) = np
           endif
      endif
   enddo

 ! print results
      do n = 1, ntrace
           call get_tracer_names (MODEL_ATMOS, n, name)
           write (stdlog(),10) n, trim(name), &
                    ( trim(echo_advec_scheme(Control%scheme(m,n))),m=1,2), &
                                             Control%weight(n),            &
                                             Control%npass_horiz(n),       &
                                             Control%npass_vert(n)
        10 format (i3, a24, ', HORIZ=',a24, ', VERT=',a24, ', wt=',f7.3, &
                   ', hp=',i2, ', vp=',i2)
      enddo

 ! check if fourth order centered horizontal scheme
   Control%do_mask4_vel = Control%scheme(1,-1) == FOURTH_CENTERED
   Control%do_mask4_tmp = .false.
   do n = 0, ntrace
     if (Control%scheme(1,n) == FOURTH_CENTERED) then
        Control%do_mask4_tmp = .true.
        exit
     endif
   enddo
 ! check if finite_volume horizontal scheme
   Control%do_finite_volume = .false.
   do n = 0, ntrace
     if (Control%scheme(1,n) == FINITE_VOLUME_LINEAR .or. &
         Control%scheme(1,n) == FINITE_VOLUME_PARABOLIC) then
         Control%do_finite_volume = .true.
         exit
     endif
   enddo
 ! error checks
 ! many finite volume schemes are not support with this release
   do n = -1, ntrace
     if (Control%scheme(1,n) == FINITE_VOLUME_LINEAR    .or. &
         Control%scheme(1,n) == FINITE_VOLUME_PARABOLIC .or. &
         Control%scheme(2,n) == FINITE_VOLUME_PARABOLIC) then
         call error_mesg ('bgrid_advection_mod',  &
                   'advection scheme not supported', FATAL)
     endif
   enddo

! initialize code sections for performance timing 

!-----------------------------------------------------------------------

  end function advection_init

!#######################################################################

 function set_default_advec ( order ) result ( advec_default )
 integer, intent(in) :: order
 integer :: advec_default(2)

 ! set horiz and vert scheme depending on input order
       select case (order)
           case(-4)
              advec_default = (/ FINITE_VOLUME_LINEAR, FINITE_VOLUME_LINEAR /)
           case(-2)
              advec_default = (/ SECOND_CENTERED, FINITE_VOLUME_LINEAR /)
           case(2)
              advec_default = (/ SECOND_CENTERED, SECOND_CENTERED /)
           case(4)
              advec_default = (/ FOURTH_CENTERED, FOURTH_CENTERED /)
           case(0)
              advec_default = (/ NONE, NONE /)
           case default
              call error_mesg  ('bgrid_advection_mod', &
                           'invalid default advection scheme', FATAL)
       end select

 end function set_default_advec
!-----------------------------------------------

 function set_advec_scheme ( scheme ) result ( advec_scheme )
 character(len=*), intent(in)  :: scheme
 integer                       :: advec_scheme

      if (uppercase(trim(scheme)) == 'SECOND_CENTERED') then
          advec_scheme = SECOND_CENTERED
      else if (uppercase(trim(scheme)) == 'FOURTH_CENTERED') then
          advec_scheme = FOURTH_CENTERED
      else if (uppercase(trim(scheme)) == 'FINITE_VOLUME_LINEAR') then
          advec_scheme = FINITE_VOLUME_LINEAR
      else if (uppercase(trim(scheme)) == 'FINITE_VOLUME_PARABOLIC') then
          advec_scheme = FINITE_VOLUME_PARABOLIC
      else if (uppercase(trim(scheme)) == 'NONE') then
          advec_scheme = NONE
      else
          call error_mesg ('bgrid_advection_mod',  &
               'invalid advection scheme, '//uppercase(trim(scheme)), FATAL)
      endif

 end function set_advec_scheme
!-----------------------------------------------

 function echo_advec_scheme ( scheme_number ) result ( scheme_name )
 integer, intent(in)  :: scheme_number
 character(len=24)    :: scheme_name

   select case (scheme_number)
       case(NONE)
            scheme_name = 'none                    '
       case(SECOND_CENTERED)
            scheme_name = 'second_centered         '
       case(FOURTH_CENTERED)
            scheme_name = 'fourth_centered         '
       case(FINITE_VOLUME_LINEAR)
            scheme_name = 'finite_volume_linear    '
       case(FINITE_VOLUME_PARABOLIC)
            scheme_name = 'finite_volume_parabolic '
       case default
            call error_mesg  ('bgrid_advection_mod', &
                     'invalid advection scheme number', FATAL)
   end select

 end function echo_advec_scheme

!#######################################################################
!         routines to minimize negative tracer values 
!        conservative horizontal and vertical borrowing
!#######################################################################

subroutine horiz_borrow (Hgrid, dt, dpde, mask, var, var_dt,  &
                         iters, eps)

!-----------------------------------------------------------------------
!
!       removes negative values by borrowing from horizontal
!       neighbors.  (usually for specific humidity or tke )
!
!-----------------------------------------------------------------------
   type(horiz_grid_type), intent(inout)                      :: Hgrid
   real(r8), intent(in)                                          :: dt
   real(r8), intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:)   :: dpde, &
                                                                mask
   real(r8), intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:,:) :: var
   real(r8), intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:,:) :: var_dt
integer, intent(in), optional                                :: iters(:)
   real(r8), intent(in), optional                                :: eps
!-----------------------------------------------------------------------
!  ( note:   0.0 < flt <= 0.125 )
   real(r8), parameter :: flt  = 0.125
   real(r8), parameter :: flt4 = flt*0.25
!-----------------------------------------------------------------------

   real(r8), dimension(lbound(var,1):ubound(var,1),                 &
                   lbound(var,2):ubound(var,2),size(var,3)) ::  &
                                hew, hns, rap, few, fns, rcur, rdif

   real(r8), dimension(lbound(var,1):ubound(var,1), &
                   lbound(var,2):ubound(var,2)) :: rew, rns,  &
                                                   areax, areay

integer :: i, j, k, n, is, ie, js, je, nlev, ntrace, knt, it, mxknt(size(var,4))
integer :: is0, ie0, js0, je0, norder
logical :: last, do_halos
   real(r8) :: rmin, hsign
!-----------------------------------------------------------------------

   mxknt = 1;  if (present(iters)) mxknt = iters
   rmin = 0.0; if (present(eps))   rmin  = eps
   if (maxval(mxknt) == 0) return

   is0 = Hgrid%Tmp%is-(Hgrid%ihalo-1);  ie0 = Hgrid%Tmp%ie+(Hgrid%ihalo-1)
   js0 = Hgrid%Tmp%js-(Hgrid%jhalo-1);  je0 = Hgrid%Tmp%je+(Hgrid%jhalo-1)

   is = is0;  ie = ie0;  js = js0;  je = je0

   nlev   = size(var,3)
   ntrace = size(var,4)

!-----------------------------------------------------------------------
!------ compute flux coeff common to all variables -------

   do j = js-1, je
   do i = is-1, ie
     areax(i,j) = (Hgrid%Tmp%area(i,j)+Hgrid%Tmp%area(i+1,j))
     areay(i,j) = (Hgrid%Tmp%area(i,j)+Hgrid%Tmp%area(i,j+1))
   enddo
   enddo

   do k = 1, nlev
   do j = js-1, je
   do i = is-1, ie
     hew(i,j,k) = areax(i,j)*(dpde(i,j,k)+dpde(i+1,j,k))* &
                              mask(i,j,k)*mask(i+1,j,k)
     hns(i,j,k) = areay(i,j)*(dpde(i,j,k)+dpde(i,j+1,k))* &
                              mask(i,j,k)*mask(i,j+1,k)
   enddo
   enddo
   enddo

   do k = 1, nlev
     rap(:,js:je,k) = flt4/(Hgrid%Tmp%area(:,js:je)*dpde(:,js:je,k))
   enddo

!-----------------------------------------------------------------------
!---- variable loop -----
!---- store current value of tracer in rdif ----

   do n = 1, ntrace

      rcur(:,:,:) = var(:,:,:,n) + var_dt(:,:,:,n)*dt

!---- iteration loop -----

      few = 0.0;  fns = 0.0
      is = is0;  ie = ie0;  js = js0;  je = je0

   do knt = 1, mxknt(n)

      rdif(:,:,:) = rcur(:,:,:)

!--------------masking lat/lon fluxes-----------------------------------

   do k = 1, nlev

!     --- do borrowing where adjacent values have opposite sign ---
!            but do not turn off fluxes previously turned on

      do j = js-1, je
      do i = is-1, ie
         if ((rdif(i+1,j,k) <  rmin .and. rdif(i,j,k) >= rmin) .or. &
             (rdif(i+1,j,k) >= rmin .and. rdif(i,j,k) <  rmin))     &
         few(i,j,k) = hew(i,j,k)
      enddo
      enddo

      do j = js-1, je
      do i = is-1, ie
         if ((rdif(i,j+1,k) <  rmin .and. rdif(i,j,k) >= rmin) .or. &
             (rdif(i,j+1,k) >= rmin .and. rdif(i,j,k) <  rmin))     &
         fns(i,j,k) = hns(i,j,k)
      enddo
      enddo

   enddo

!---- loop for fourth order iteration ----
!---- SWITCHED TO SECOND ORDER (better filling) ----

   norder = 1
   hsign  = 1

   do it = 1, norder

      last = knt == mxknt(n) .and. it == norder

!--------------2-nd order lat/lon contributions-------------------------

      do k = 1, nlev
         do j = js-1, je
         do i = is-1, ie
            rew(i,j) = (rdif(i+1,j,k)-rdif(i,j,k))*few(i,j,k)
            rns(i,j) = (rdif(i,j+1,k)-rdif(i,j,k))*fns(i,j,k)
         enddo
         enddo

         do j = js, je
         do i = is, ie
            rdif(i,j,k) = (rew(i,j)-rew(i-1,j)+ &
                           rns(i,j)-rns(i,j-1))*rap(i,j,k)
         enddo
         enddo
      enddo

!     ---- halo update when necessary (but not on last pass) ----

                                                          do_halos = .false.
      if ( is  == Hgrid%Tmp%is .or.  js == Hgrid%Tmp%js ) do_halos = .true.
      if ( knt == mxknt(n)     .and. it == norder       ) do_halos = .false.

      if ( do_halos ) then
         call update_halo (Hgrid, TEMP, rdif)
         is = is0;  ie = ie0;  js = js0;  je = je0
      else
         call update_halo (Hgrid, TEMP, rdif, POLEONLY)
         is = is+1;  ie = ie-1
         js = js+1;  je = je-1
      endif

   enddo

!  --- update current value of tracer ---

   rcur(:,:,:) = rcur(:,:,:) + hsign*rdif(:,:,:)

!-----------------------------------------------------------------------
 enddo
!-----------------------------------------------------------------------
!---- return the tendency -----

   call update_halo (Hgrid, TEMP, rcur)

   var_dt(:,:,:,n) = (rcur(:,:,:) - var(:,:,:,n)) / dt

 enddo

!-----------------------------------------------------------------------

end subroutine horiz_borrow

!#######################################################################

   subroutine vert_borrow (dt, dpde, var, var_dt, eps, iters)

!-----------------------------------------------------------------------
!
!  This removes negative specific humidity/mixing ratios by borrowing
!  from the grid boxes immediately above and below. If not enough
!  is available to fill the negative then a negative will remain.
!
!-----------------------------------------------------------------------
   real(r8), intent(in)    :: dt, dpde(:,:,:), var(:,:,:,:)
   real(r8), intent(inout) :: var_dt(:,:,:,:)
   real(r8), intent(in), optional :: eps
integer, intent(in), optional :: iters(:)
!-----------------------------------------------------------------------
   real(r8), dimension(size(var,1),size(var,2),size(var,3)) ::  &
                                var_new, deficit, surplus, rdpdt, small
   real(r8), dimension(size(var,1),size(var,2)) :: divid, ratio_dn,  &
                                               ratio_up, var_dp
   real(r8)    :: epsil
   integer :: k, n, kdim, iter, num_iters(size(var,4)), num, num_var
!-----------------------------------------------------------------------

   num_iters = 4; if (present(iters)) num_iters = iters
   if (maxval(num_iters) == 0) return

   num_var = size(var,4)
   if (num_var == 0) return

   kdim = size(var,3)
   rdpdt = 1./(dpde*dt)
   epsil = 0.0; if (present(eps)) epsil = eps
   small = epsil * dpde

!---- variable loop and iteration loop -----

 do n = 1, num_var

 do iter = 1, num_iters(n)

!---- set var_new to minimum value -----
!---- existing negatives will not be corrected ----

   var_new = var(:,:,:,n)
   var_new = var_new + dt * var_dt(:,:,:,n)

!! not safe for mpp reproducibility
!! num = count(var_new < epsil)
!! if ( num == 0 ) exit

   do k = 1, kdim
      var_dp(:,:) = var_new(:,:,k)*dpde(:,:,k)
      deficit(:,:,k) = min(var_dp(:,:),small(:,:,k))-small(:,:,k)
      surplus(:,:,k) = max(var_dp(:,:),small(:,:,k))-small(:,:,k)
   enddo

!---- top level ----

   where (deficit(:,:,1) < 0.0)
      divid(:,:) = max(-surplus(:,:,2)/deficit(:,:,1),1.0_r8)
      ratio_dn(:,:) = surplus(:,:,2)/divid(:,:)
      var_dt(:,:,1,n) = var_dt(:,:,1,n) + ratio_dn(:,:)*rdpdt(:,:,1)
      var_dt(:,:,2,n) = var_dt(:,:,2,n) - ratio_dn(:,:)*rdpdt(:,:,2)
      surplus(:,:,2) = max(surplus(:,:,2)-ratio_dn(:,:),small(:,:,2))
   endwhere

!---- interior levels ----

   do k = 2, kdim-1
     where (deficit(:,:,k) < 0.0)
       divid(:,:) = max(-(surplus(:,:,k-1)+surplus(:,:,k+1))/  &
                        deficit(:,:,k),1.0_r8)
       ratio_up(:,:) = surplus(:,:,k-1)/divid(:,:)
       ratio_dn(:,:) = surplus(:,:,k+1)/divid(:,:)
       var_dt(:,:,k,n)  = var_dt(:,:,k,n)+(ratio_up(:,:)+ratio_dn(:,:))*rdpdt(:,:,k)
       var_dt(:,:,k-1,n) = var_dt(:,:,k-1,n)-ratio_up(:,:)*rdpdt(:,:,k-1)
       var_dt(:,:,k+1,n) = var_dt(:,:,k+1,n)-ratio_dn(:,:)*rdpdt(:,:,k+1)
       surplus(:,:,k+1) = max(surplus(:,:,k+1)-ratio_dn(:,:),small(:,:,k+1))
     endwhere
   enddo

!---- bottom level ----

   where (deficit(:,:,kdim) < 0.0)
      divid(:,:) = max(-surplus(:,:,kdim-1)/deficit(:,:,kdim),1.0_r8)
      ratio_up(:,:) = surplus(:,:,kdim-1)/divid(:,:)
      var_dt(:,:,kdim,n)   = var_dt(:,:,kdim,n)+ratio_up(:,:)*rdpdt(:,:,kdim)
      var_dt(:,:,kdim-1,n) = var_dt(:,:,kdim-1,n)-ratio_up(:,:)*rdpdt(:,:,kdim-1)
   endwhere

 enddo

 enddo

!-----------------------------------------------------------------------

   end subroutine vert_borrow

!#######################################################################

end module bgrid_advection_mod

