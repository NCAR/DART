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

module bgrid_horiz_adjust_mod

!-----------------------------------------------------------------------

use types_mod, only : r8
use bgrid_horiz_mod      , only: horiz_grid_type
use bgrid_vert_mod       , only: vert_grid_type, compute_pres_depth, &
                                 compute_geop_height, compute_pres_half
use bgrid_masks_mod      , only: grid_mask_type
use bgrid_halo_mod       , only: update_halo, TEMP
use bgrid_change_grid_mod, only: mass_to_vel

use         constants_mod, only: OMEGA, RADIUS, RDGAS, RVGAS
use               fms_mod, only: error_mesg, FATAL

implicit none
private

public :: horiz_adjust_vel, horiz_adjust_mass, press_grad, &
          div_damping, compute_grad_pres

!-----------------------------------------------------------------------

   real(r8), parameter :: fcor = 2.*OMEGA,     &
                      curv = 1.0/RADIUS,   &
                      wp   = 0.0625

   real(r8), parameter :: d608 = (RVGAS-RDGAS)/RDGAS
   real(r8), parameter :: rdgas2 = RDGAS*0.50
!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine horiz_adjust_vel ( Hgrid, Masks, dt, &
                              pgfew, pgfns, um, vm, udt, vdt,  &
                              alpha_implicit )

!-----------------------------------------------------------------------
!
!    IN: Hgrid  = horizontal constants
!        Masks  = grid mask constants
!        dt     = time step
!        pgfew,
!         pgfns = pressure gradient force components
!        um,vm  = zonal and meridional wind components
!
! INOUT: udt,vdt = tendency of zonal and meridional wind components
!
! IN (opt): alpha_implicit = coefficient for coriolis/pgf time diff
!                               0.0 = explicit (not recommended)
!                               0.5 = trapezoidal implicit (default)
!                               1.0 = fully implicit
!
!-----------------------------------------------------------------------
type(horiz_grid_type), intent(in)      :: Hgrid
type (grid_mask_type), intent(in)      :: Masks
   real(r8), intent(in)                    :: dt
   real(r8), intent(in),    dimension(Hgrid%ilb:, Hgrid%jlb:, :) :: &
                                           pgfew, pgfns, um, vm
   real(r8), intent(inout), dimension(Hgrid%ilb:, Hgrid%jlb:, :) :: &
                                                        udt, vdt
real(r8), optional, intent(in)              :: alpha_implicit
!-----------------------------------------------------------------------

  real(r8), dimension (Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub) :: &
      pgfu, pgfv, f0, fa, fb, cu, cv, up, vp

  real(r8), dimension (Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub,  &
         size(pgfew,3)) :: un, vn

    integer :: i, j, k
    real(r8)    :: alpha
!-----------------------------------------------------------------------

     alpha = 0.5; if (present(alpha_implicit)) alpha = alpha_implicit
     alpha = min(max(0.0_r8,alpha),1.0_r8)

!    --- force explicit for leapfrog (if not implicit) ----
!    if (nv == 2 .and. alpha <= 0.51) then
!        nl = 2
!        alpha = 0.0
!    endif

!-----------------------------------------------------------------------
!-----------------update u and v (coriolis & pgf)-----------------------

      do k =            1, size(um,3)
      do j = Hgrid%Vel%js, Hgrid%Vel%je
      do i = Hgrid%Vel%is, Hgrid%Vel%ie

!     --------- pressure gradient force components---------------

         pgfu(i,j) = dt*pgfew(i,j,k)
         pgfv(i,j) = dt*pgfns(i,j,k)

!     ------------coriolis & curvature terms----------------------

         f0(i,j) = (um(i,j,k)*curv*Hgrid%tanphv(i,j)+  &
                                   fcor*Hgrid%sinphv(i,j))*dt
         fa(i,j) = f0(i,j)*(1.0-alpha)
         fb(i,j) = f0(i,j)*alpha

         cu(i,j) =  vm(i,j,k)*fa(i,j)
         cv(i,j) = -um(i,j,k)*fa(i,j)

!     ------------compute new u and v (coriolis & pgf)------------

         up(i,j) = pgfu(i,j) + cu(i,j) + um(i,j,k)
         vp(i,j) = pgfv(i,j) + cv(i,j) + vm(i,j,k)

         un(i,j,k) = ((fb(i,j)*vp(i,j)+up(i,j))/   &
                      (fb(i,j)*fb(i,j)+1.0))*Masks%Vel%mask(i,j,k)
         vn(i,j,k) = (vp(i,j)-fb(i,j)*un(i,j,k))*Masks%Vel%mask(i,j,k)

!     ---- return unfiltered tendencies with halos not updated ----

         udt(i,j,k) = udt(i,j,k) + (un(i,j,k)-um(i,j,k))/dt
         vdt(i,j,k) = vdt(i,j,k) + (vn(i,j,k)-vm(i,j,k))/dt
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------

end subroutine horiz_adjust_vel

!#######################################################################

 subroutine horiz_adjust_mass ( nplev, Hgrid, Masks,  &
                                u, v, dpde, cew, cns,        &
                                flew, flns, div, omgalf      )

!-----------------------------------------------------------------------
!
!    IN: nplev  = vertical index of uppermost pure pressure level
!                 (use nplev=0 for sigma models)
!        Hgrid  = horizontal constants
!        Masks  = grid mask constants
!        u, v   = prognostic variables for zonal and meridional wind
!        cew,
!          cns  = zonal and meridional components of grad(p)/p
!
! INOUT: flew,
!          flns = zonal, meridional mass fluxes (summation)
!
!   OUT: div    = mass divergence (kg/m/s3)
!        omgalf = horizontal omega-alpha term (in part) (deg_k/s)
!
!-----------------------------------------------------------------------
integer, intent(in)                   :: nplev
type(horiz_grid_type), intent(inout)  :: Hgrid
type (grid_mask_type), intent(in)     :: Masks
  real(r8), intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: u, v, &
                                                        dpde, cew, cns
  real(r8), intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: flew, flns
  real(r8), intent(out),   dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: div, omgalf
!-----------------------------------------------------------------------

    real(r8), dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) ::   &
                                  few, fns, udy, vdx, tew, tns

    real(r8), dimension(Hgrid%ilb:Hgrid%iub,                &
                    Hgrid%jlb:Hgrid%jub, size(u,3)) ::  adpdxy

    integer :: i, j, k, is, ie, js, je


    is = Hgrid%Tmp%is;  ie = Hgrid%Tmp%ie
    js = Hgrid%Tmp%js;  je = Hgrid%Tmp%je

!-----------------------------------------------------------------------
!  ---- average mass weights for mass fluxes ----

   call mass_to_vel ( Hgrid, dpde, adpdxy )

!!!call press_smooth (Hgrid, nplev, 0.25, adpdxy, adpdxy )

!-----------------------------------------------------------------------

      few = 0.0;  fns = 0.0

   do k = 1, size(u,3)

!---- compute mass fluxes, divergence & horizontal omega-alpha term ----
!-------add in compute flux corrections for grid separation ------------
!---- sum input/output mass fluxes -------

      udy(:,:)=Hgrid%Vel%dy     *u(:,:,k)*adpdxy(:,:,k)
      vdx(:,:)=Hgrid%Vel%dx(:,:)*v(:,:,k)*adpdxy(:,:,k)

    ! fix for fluxes at sub-pole row
    ! if (Hgrid%Vel%js == Hgrid%Vel%jsg) udy(:,js-1)=udy(:,js)
    ! if (Hgrid%Vel%je == Hgrid%Vel%jeg) udy(:,je  )=udy(:,je-1)

!     ---- without grid separation modification ---

      do j = js,   je
      do i = is-1, ie
          few(i,j)   = (udy(i,j)+udy(i,j-1))*0.5
          tew(i,j)   = few(i,j)*cew(i,j,k)
      enddo
      enddo
      do j = js-1, je
      do i = is,   ie
          fns(i,j)   = (vdx(i,j)+vdx(i-1,j))*0.5
          tns(i,j)   = fns(i,j)*cns(i,j,k)
      enddo
      enddo

!-----------------------------------------------------------------------
! ------ sum fluxes (output needed for advection) ------
! ------ (halos will be updated in advection) ------

      flew(:,:,k) = flew(:,:,k) + few
      flns(:,:,k) = flns(:,:,k) + fns

!---------------------------divergence----------------------------------

      do j = js, je
      do i = is, ie
         div(i,j,k) = ((few(i  ,j)+fns(i,j  ))-  &
                       (few(i-1,j)+fns(i,j-1)))  &
                        *Hgrid%Tmp%rarea(i,j)*Masks%Tmp%mask(i,j,k)
      enddo
      enddo

!-------------- horizontal part of omega-alpha -------------------------
!        ------ do not do for pure pressure levels -----

      if (k > nplev) then
         do j = js, je
         do i = is, ie
            omgalf(i,j,k)=(tew(i,j)+tew(i-1,j)+tns(i,j)+tns(i,j-1)) &
                          *0.50*Hgrid%Tmp%rarea(i,j)*Masks%Tmp%mask(i,j,k)
         enddo
         enddo
      else
            omgalf(:,:,k) = 0.0
      endif

!-----------------------------------------------------------------------
!---- end level loop -----

   enddo

!-----------------------------------------------------------------------

 end subroutine horiz_adjust_mass

!#######################################################################

subroutine press_grad ( Hgrid, Vgrid, Masks, fssl, tq,         &
                        dpde, wta, wtb, cew, cns, pgfew, pgfns )

!-----------------------------------------------------------------------
!
!           Computes geopotential height, alpha (RT/p),
!                 and pressure at model levels.
!
!    IN: Hgrid  = horizontal constants
!        Vgrid  = horizontal constants
!        Masks  = grid mask constants
!        fssl   = geopotential height (m2/s2) at eta=1.
!        tq     = virtual temperature
!        dpde   = pressure thickness of model layers
!        wta,
!          wtb  = weights for computing geopotential height
!                 (same as weight for computing full pressures)
!        cew,
!          cns  = zonal and meridional components of grad(p)/p
!
!   OUT: pgfew,
!         pgfns = pressure gradient force components
!
!-----------------------------------------------------------------------
type(horiz_grid_type), intent (inout) :: Hgrid
type (vert_grid_type), intent (in)  :: Vgrid
type (grid_mask_type), intent (in)  :: Masks
   real(r8), intent (in) , dimension(Hgrid%ilb:Hgrid%iub, &
                                 Hgrid%jlb:Hgrid%jub) :: fssl
   real(r8), intent (in),  dimension(Hgrid%ilb:Hgrid%iub, &
                                 Hgrid%jlb:Hgrid%jub,Vgrid%nlev) :: &
                                       tq, dpde, wta, wtb, cew, cns
   real(r8), intent (out), dimension(Hgrid%ilb:Hgrid%iub, &
                                 Hgrid%jlb:Hgrid%jub,Vgrid%nlev) :: &
                                       pgfew, pgfns

!-----------------------------------------------------------------------
  real(r8), dimension (Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: &
                                       pew, pns, tew, tns

  real(r8), dimension (Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub,  &
                   size(tq,3)) :: fim, ppcew, ppcns

!!real(r8), dimension (Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub,  &
!!                 size(tq,3)) :: adpdxy, pwt

    integer :: i, j, k, is, ie, js, je
!-----------------------------------------------------------------------

!  ---- compute pressure depths using appropriate pressure ----
!          dpde  = model layer pressure thicknesses
!
!     note: forward-backward scheme uses fields at tau+1
!           leap-frog scheme uses time averaged fields at tau


!-----------------------------------------------------------------------
!------------- integration of geopotential height-----------------------
!        fim = geopotential height at model levels

   if (Vgrid%nlev /= size(tq,3)) call error_mesg ('geop_height',  &
                              'wrong number of vertical levels', FATAL)

if (Vgrid%nlev > 1) then

   call compute_geop_height (Vgrid, fssl, tq, wta, wtb, fim, &
                             mask=Masks%Tmp%mask)

else

      fim(:,:,:) = dpde(:,:,:)

endif

!----------- lat/lon contributions to pressure gradient ----------------
!   -------- two loops: pressure only, pressure/sigma -------
!     ppcew,ppcns = zonal and meridional auxillary pressure gradient
!                   force components

   is = Hgrid%Vel%is;  ie = Hgrid%Vel%ie
   js = Hgrid%Vel%js;  je = Hgrid%Vel%je

   do k =  1, Vgrid%nplev
     do j = js, je+1
     do i = is, ie
       ppcew(i,j,k) = fim(i+1,j,k)-fim(i,j,k)
     enddo
     enddo
     do j = js, je
     do i = is, ie+1
       ppcns(i,j,k) = fim(i,j+1,k)-fim(i,j,k)
     enddo
     enddo
   enddo

   do k = Vgrid%nplev+1, Vgrid%nlev
     do j = js, je+1
     do i = is, ie
         pew(i,j)   = fim(i+1,j,k)-fim(i,j,k)
         tew(i,j)   = rdgas2*(tq(i+1,j,k)+tq(i,j,k))
       ppcew(i,j,k) = pew(i,j)+tew(i,j)*cew(i,j,k)
     enddo
     enddo
     do j = js, je
     do i = is, ie+1
         pns(i,j)   = fim(i,j+1,k)-fim(i,j,k)
         tns(i,j)   = rdgas2*(tq(i,j+1,k)+tq(i,j,k))
       ppcns(i,j,k) = pns(i,j)+tns(i,j)*cns(i,j,k)
     enddo
     enddo
   enddo

!--------------compute pressure gradient force components---------------

!!!call mass_to_vel  (Hgrid, dpde, adpdxy )
!!!call press_smooth (Hgrid, Vgrid%nplev, 0.25, adpdxy, pwt )
!!!pwt = pwt/adpdxy

   do k =  1, Vgrid%nlev
   do j = js, je
   do i = is, ie
      pgfew(i,j,k)=-0.50*(ppcew(i,j,k)+ppcew(i,j+1,k))*Hgrid%Vel%rdx(i,j)
      pgfns(i,j,k)=-0.50*(ppcns(i,j,k)+ppcns(i+1,j,k))*Hgrid%Vel%rdy
!     ---- smoothing turned off ???? ----
!!!!  pgfew(i,j,k)=pgfew(i,j,k)*pwt(i,j,k)
!!!!  pgfns(i,j,k)=pgfns(i,j,k)*pwt(i,j,k)
   enddo
   enddo
   enddo

!-----------------------------------------------------------------------

 end subroutine press_grad

!#######################################################################

 subroutine compute_grad_pres (Hgrid, nplev, phalf, dpde, &
                               wta, wtb, cew, cns)

   type(horiz_grid_type), intent(in)  :: Hgrid
   integer,               intent(in)  :: nplev
   real(r8), intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) ::  &
                                       phalf, dpde, wta, wtb
   real(r8), intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: cew, cns

   real(r8), dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: &
                                              pewa, pewb, pnsa, pnsb
   real(r8), dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub, &
                                        size(dpde,3)) :: adpdx, adpdy
   integer :: i, j, k

 ! special case for one level model
   if (size(dpde,3) == 1) then
       call error_mesg ('compute_grad_pres', &
                 'single level untested - contact DART developers', FATAL)
       cew(:,:,k) = 0.0  ! DART ... k is unset at this point ... ?value?
       cns(:,:,k) = 0.0  ! DART ... k is unset at this point ... ?value?
       return
   endif

 ! term vanishes on pressure levels 
   do k = 1, nplev
       cew(:,:,k) = 0.0
       cns(:,:,k) = 0.0
   enddo

   call compute_mass_weights (Hgrid, dpde, adpdx, adpdy)

! --- compute pressure differences ----

   do j = Hgrid%Tmp%js,   Hgrid%Tmp%je+1
   do i = Hgrid%Tmp%is-1, Hgrid%Tmp%ie
      pewa(i,j) = phalf(i+1,j,nplev+1)-phalf(i,j,nplev+1)
   enddo
   enddo

   do j = Hgrid%Tmp%js-1, Hgrid%Tmp%je
   do i = Hgrid%Tmp%is,   Hgrid%Tmp%ie+1
      pnsa(i,j) = phalf(i,j+1,nplev+1)-phalf(i,j,nplev+1)
   enddo
   enddo

! --- compute grad pressure term ----

   do k = nplev+1, size(phalf,3)-1
!      ---- east-west contribution ----
       do j = Hgrid%Tmp%js,   Hgrid%Tmp%je+1
       do i = Hgrid%Tmp%is-1, Hgrid%Tmp%ie
           pewb(i,j) = phalf(i+1,j,k+1)-phalf(i,j,k+1)
           cew(i,j,k) = (wta(i+1,j,k)+wta(i,j,k)) * pewa(i,j) + &
                        (wtb(i+1,j,k)+wtb(i,j,k)) * pewb(i,j)
           cew(i,j,k) = cew(i,j,k) / adpdx(i,j,k)
           pewa(i,j) = pewb(i,j)
       enddo
       enddo
!      ---- north-south contribution ----
       do j = Hgrid%Tmp%js-1, Hgrid%Tmp%je
       do i = Hgrid%Tmp%is,   Hgrid%Tmp%ie+1
           pnsb(i,j) = phalf(i,j+1,k+1)-phalf(i,j,k+1)
           cns(i,j,k) = (wta(i,j+1,k)+wta(i,j,k)) * pnsa(i,j) + &
                        (wtb(i,j+1,k)+wtb(i,j,k)) * pnsb(i,j)
           cns(i,j,k) = cns(i,j,k) / adpdy(i,j,k)
           pnsa(i,j) = pnsb(i,j)
       enddo
       enddo
   enddo

 end subroutine compute_grad_pres

!#######################################################################

 subroutine div_damping (Hgrid, Vgrid, Masks, coeff, adpdxy, div,  &
                         u_dt, v_dt)

!-----------------------------------------------------------------------
!
!       Computes divergence damping tendency for momentum fields
!
!    IN: Hgrid  = horizontal constants
!        Vgrid  = horizontal constants
!        Masks  = grid mask constants
!        dpde   = model layer pressure thicknesses
!        div    = mass divergence (units = mass/sec)
!
! INOUT: u_dt,
!          v_dt = zonal and meridional momentum tendencies
!
!-----------------------------------------------------------------------

type(horiz_grid_type), intent (in)  :: Hgrid
type (vert_grid_type), intent (in)  :: Vgrid
type (grid_mask_type), intent (in)  :: Masks
 real(r8), intent (in)                  :: coeff
 real(r8), intent (in) ,   dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: adpdxy, div
 real(r8), intent (inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: u_dt, v_dt

!-----------------------------------------------------------------------

  real(r8), dimension (lbound(div,1):ubound(div,1),    &
                   lbound(div,2):ubound(div,2)) :: &
                            dcoeff, ddx, ddy, dew, dns

  integer :: i, j, k, is, ie, js, je, ke
!-----------------------------------------------------------------------

  is = Hgrid % Vel % is;  ie = Hgrid % Vel % ie
  js = Hgrid % Vel % js;  je = Hgrid % Vel % je
  ke = size(div,3)

  do k = 1, ke

     dcoeff(:,:) = coeff * Hgrid % Vel % rarea(:,:) / adpdxy(:,:,k)
     ddx(:,:) = div(:,:,k) * Hgrid % Tmp % dx(:,:)
     ddy(:,:) = div(:,:,k) * Hgrid % Tmp % dy

     do j = js, je+1
     do i = is, ie
        dew(i,j) = ddx(i+1,j) + ddx(i,j)
     enddo
     enddo

     do j = js, je
     do i = is, ie+1
        dns(i,j) = ddy(i,j+1) + ddy(i,j)
     enddo
     enddo
   
     do j = js, je
     do i = is, ie
        u_dt(i,j,k) = u_dt(i,j,k) + (dns(i+1,j)-dns(i,j)) * dcoeff(i,j)
        v_dt(i,j,k) = v_dt(i,j,k) + (dew(i,j+1)-dns(i,j)) * dcoeff(i,j)
     enddo
     enddo

  enddo

!-----------------------------------------------------------------------

 end subroutine div_damping

!#######################################################################

 subroutine press_smooth (Hgrid, nplev, w, ps, psm )

   type(horiz_grid_type), intent(inout) :: Hgrid
   integer, intent(in)                :: nplev
   real(r8),    intent(in)                :: w
   real(r8),    intent(in)                :: ps (Hgrid%ilb:,Hgrid%jlb:,:)
   real(r8),    intent(out)               :: psm(Hgrid%ilb:,Hgrid%jlb:,:)

!-----------------------------------------------------------------------

    real(r8), dimension(lbound(ps,1):ubound(ps,1),      &
                    lbound(ps,2):ubound(ps,2)) ::   &
                          adx, ady, dpdx, dpdy, dpcew, dpcns, pew, pns

    real(r8) :: wpdt
    integer :: i, j, k, is, ie, js, je

!-----------------------------------------------------------------------

    is = Hgrid%Tmp%is;  ie = Hgrid%Tmp%ie
    js = Hgrid%Tmp%js;  je = Hgrid%Tmp%je

    wpdt = -w*wp

!   --- area weights for fluxes ---
      do j = js-1, je+1
      do i = is-1, ie
         adx(i,j) = Hgrid%Tmp%area(i+1,j)+Hgrid%Tmp%area(i,j)
      enddo
      enddo
      do j = js-1, je
      do i = is-1, ie+1
         ady(i,j) = Hgrid%Tmp%area(i,j+1)+Hgrid%Tmp%area(i,j)
      enddo
      enddo

      do k = nplev+1, size(ps,3)

         do j = js-1, je+1
         do i = is-1, ie
            dpdx(i,j) = (ps(i+1,j,k)-ps(i,j,k))*adx(i,j)
         enddo
         enddo

         do j = js-1, je
         do i = is-1, ie+1
            dpdy(i,j) = (ps(i,j+1,k)-ps(i,j,k))*ady(i,j)
         enddo
         enddo

         do j = js-1,je
         do i = is-1,ie
            dpcew(i,j) = (dpdx(i,j+1)-dpdx(i,j))
            dpcns(i,j) = (dpdy(i+1,j)-dpdy(i,j))
         enddo
         enddo

         do j = js,   je
         do i = is-1, ie
            pew(i,j) = dpcew(i,j)-dpcew(i,j-1)
         enddo
         enddo

         do j = js-1, je
         do i = is,   ie
            pns(i,j) = dpcns(i,j)-dpcns(i-1,j)
         enddo
         enddo

         do j = js, je
         do i = is, ie
            psm(i,j,k) = ps(i,j,k) + wpdt*Hgrid%Tmp%rarea(i,j)* &
                              (pew(i,j)-pew(i-1,j)+pns(i,j)-pns(i,j-1))
         enddo
         enddo

      enddo

      psm(:,:,1:nplev) = ps(:,:,1:nplev)

      call update_halo (Hgrid, TEMP, psm(:,:,nplev+1:size(ps,3)))

!-----------------------------------------------------------------------

 end subroutine press_smooth

!#######################################################################

 subroutine compute_mass_weights (Hgrid, dpde, adpdx, adpdy, adpdxy)

   type(horiz_grid_type), intent(in)     :: Hgrid
   real(r8), intent(in)            :: dpde  (Hgrid%ilb:,Hgrid%jlb:,:)
   real(r8), intent(out)           :: adpdx (Hgrid%ilb:,Hgrid%jlb:,:),  &
                                  adpdy (Hgrid%ilb:,Hgrid%jlb:,:)
   real(r8), intent(out), optional :: adpdxy(Hgrid%ilb:,Hgrid%jlb:,:)

   integer :: i, j, k

!-------- average mass weights -----------------------------------------

   do k = 1, size(dpde,3)

      do j = Hgrid%jlb, Hgrid%jub
      do i = Hgrid%ilb, Hgrid%iub-1
        adpdx(i,j,k) = dpde(i,j,k)+dpde(i+1,j,k)
      enddo
      enddo
      do j = Hgrid%jlb, Hgrid%jub-1
      do i = Hgrid%ilb, Hgrid%iub
        adpdy(i,j,k)=(Hgrid%Tmp%area(i,j  )*dpde(i,j  ,k)+ &
                      Hgrid%Tmp%area(i,j+1)*dpde(i,j+1,k)) &
                     *Hgrid%Vel%rarea(i,j)
      enddo
      enddo
   enddo

!------- average mass weights at velocity points -----------------------

 if (present(adpdxy)) then
   do k = 1, size(dpde,3)
      do j = Hgrid%jlb, Hgrid%jub-1
      do i = Hgrid%ilb, Hgrid%iub-1
         adpdxy(i,j,k) = adpdy(i,j,k)+adpdy(i+1,j,k)
      enddo
      enddo
   enddo
!     --- avoid undef elements ---
      adpdxy(Hgrid%iub,:,:) = 0.0
      adpdxy(:,Hgrid%jub,:) = 0.0
 endif

!-----------------------------------------------------------------------

 end subroutine compute_mass_weights

!#######################################################################

 function time_average2 (coeff, dt, um, u, udt) result (ua)

   real(r8), intent(in)  :: coeff, dt, um(:,:), u(:,:), udt(:,:)
   real(r8) :: ua(size(u,1),size(u,2))
   real(r8) :: coeff2

!  time averages tau-1, tau, tau+1 time levels
!  used for brown-campana (pgf) filter

   coeff2 = 1.-2.*coeff
   ua = coeff * ( 2. * um + udt*dt ) + coeff2 * u

 end function time_average2

!#######################################################################

 function time_average3 (coeff, dt, um, u, udt) result (ua)

   real(r8), intent(in)  :: coeff, dt, um(:,:,:), u(:,:,:), udt(:,:,:)
   real(r8) :: ua(size(u,1),size(u,2),size(u,3))
   real(r8) :: coeff2

!  time averages tau-1, tau, tau+1 time levels
!  used for brown-campana (pgf) filter

   coeff2 = 1.-2.*coeff
   ua = coeff * ( 2. * um + udt*dt ) + coeff2 * u

 end function time_average3

!#######################################################################

end module bgrid_horiz_adjust_mod

