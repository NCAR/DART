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

module bgrid_sponge_mod

!-----------------------------------------------------------------------
!
!   (1) eddy damping of prognostic fields (t,q,u,v) at the top
!       level of the model (i.e., fields are damped with their
!       zonal mean values). THIS OPTION WILL NOT WORK CORRECTLY
!       WITH DOMAIN DECOMPOSITION ALONG THE X-AXIS.
!
!-----------------------------------------------------------------------

use types_mod, only : r8
 use bgrid_horiz_mod,       only: horiz_grid_type
 use bgrid_masks_mod,       only: grid_mask_type
 use bgrid_prog_var_mod,    only: prog_var_type
 use bgrid_change_grid_mod, only: mass_to_vel
 use bgrid_halo_mod,        only: update_halo, vel_flux_boundary, &
                                  TEMP, UWND, VWND, &
                                  NORTH, EAST, NOPOLE

 use         fms_mod, only: error_mesg, FATAL, write_version_number
                            

! use mpp_domains_mod, only: domain2d, mpp_redistribute, &
!                            mpp_define_domains,         &
!                            mpp_get_global_domain,      &
!                            mpp_get_data_domain,        &
!                            CYCLIC_GLOBAL_DOMAIN

 implicit none
 private

 public   sponge_driver, sponge_init

!-----------------------------------------------------------------------

 real(r8), parameter ::  daypsec=1./86400.
 logical :: do_init=.true.

 logical :: do_topsponge
 integer :: nlev_sponge
 real(r8)    :: dfactr   ! coeff for damping eddies at the top level


 character(len=128) :: version='$Revision$'
 character(len=128) :: tagname='$Id$'
 logical :: do_log = .true.

 real(r8) :: small = .000001

 real(r8) :: coeff_vel, coeff_tmp, coeff_trs
 real(r8) :: xaxis_wt
 integer :: numlev

!--- module storage for computing exact/reproducible zonal means ---
  type zsum_type
     !!!type(domain2d) :: Domain
     integer        :: isize, jsize
  end type
  type(zsum_type) :: Zdomain_tmp, Zdomain_vel
  logical :: first = .true.
!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine sponge_driver ( Hgrid, nplev, dt, dpde, Var, Var_dt )
 
!-----------------------------------------------------------------------
 type (horiz_grid_type), intent(inout) :: Hgrid
 integer,                intent(in)    :: nplev
 real(r8),                   intent(in)    :: dt
 real(r8),                   intent(in)    :: dpde(Hgrid%ilb:,Hgrid%jlb:,:)
 type  (prog_var_type),  intent(in)    :: Var
 type  (prog_var_type),  intent(inout) :: Var_dt
!-----------------------------------------------------------------------

 real(r8), dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub,numlev) :: dpxy
 integer :: n, np, nt

!-----------------------------------------------------------------------

   if (do_init) call error_mesg ('bgrid_sponge_mod',  &
                                 'initialization not called', FATAL)

   if ( numlev == 0 ) return

!  do not allow decomposition along x-axis
!  if ( Hgrid%decompx ) call error_mesg ('bgrid_sponge_mod', &
!    'this version of sponge does not allow x-axis decomposition.', FATAL)

!-----------------------------------------------------------------------

   n  = numlev
   np = nplev

!---- temperature and tracers -----

   if ( coeff_tmp > small .or. coeff_trs > small ) then
       if (first) then
          Zdomain_tmp%isize = Hgrid%Tmp%ie - Hgrid%Tmp%is + 1
          Zdomain_tmp%jsize = Hgrid%Tmp%je - Hgrid%Tmp%js + 1
       endif
       !!!if (first) call zsum_init ( Hgrid%Tmp%Domain, Zdomain_tmp )
       nt = Var_dt%ntrace
       call local_filter_mass ( Hgrid, np, coeff_tmp, coeff_trs, dt,       &
                dpde(:,:,1:n), Var   %t(:,:,1:n), Var   %r(:,:,1:n,1:nt),  &
                               Var_dt%t(:,:,1:n), Var_dt%r(:,:,1:n,1:nt) )
   endif

!---- momentum components -----

   if ( coeff_vel > small ) then
        if (first) then
          Zdomain_vel%isize = Hgrid%Vel%ie - Hgrid%Vel%is + 1
          Zdomain_vel%jsize = Hgrid%Vel%je - Hgrid%Vel%js + 1
        endif
        !!!if (first) call zsum_init ( Hgrid%Vel%Domain, Zdomain_vel )
      ! compute pressure weights at velocity points
        dpxy(:,:,:) = dpde(:,:,1:n)
        if (np < n) then
          call mass_to_vel (Hgrid, dpxy(:,:,np+1:n), dpxy(:,:,np+1:n))
          call update_halo (Hgrid, UWND, dpxy(:,:,np+1:n), EAST+NORTH+NOPOLE)
        endif

        call local_filter_vel ( Hgrid, np, coeff_vel, dt, dpxy,       &
                                Var   %u(:,:,1:n), Var   %v(:,:,1:n), &
                                Var_dt%u(:,:,1:n), Var_dt%v(:,:,1:n)  )
   endif

   if (first) first = .false.

!-----------------------------------------------------------------------

 end subroutine sponge_driver

!#######################################################################

 subroutine sponge_init ( damp_vel, damp_tmp, damp_trs, nlev, xwt )

   real(r8),    intent(in), optional :: damp_vel, damp_tmp, damp_trs, xwt
   integer, intent(in), optional :: nlev

!
!   damp_vel, damp_tmp, damp_trs =
!         normalized (0.,1.) damping coefficients for
!         momentum, temperature, and tracers
!         [real(r8), default = 0.]
!
!   nlev   = number of levels at the top of the model where
!            damping is performed
!               [integer, default = 0]
!
!   xwt    = weight (0.,1.) applied to the x-axis
!            the y-axis will have a weight = 1.-xwt
!               [real(r8), default = .50]
!

!-----------------------------------------------------------------------

   if (do_log) then
      call write_version_number (version,tagname)
      do_log = .false.
   endif

! set values for optional arguments

   coeff_vel = 0.; if (present(damp_vel)) coeff_vel = min(max(damp_vel,0.0_r8),1.0_r8)
   coeff_tmp = 0.; if (present(damp_tmp)) coeff_tmp = min(max(damp_tmp,0.0_r8),1.0_r8)
   coeff_trs = 0.; if (present(damp_trs)) coeff_trs = min(max(damp_trs,0.0_r8),1.0_r8)

   numlev   = 0;   if (present(nlev)) numlev   = max(nlev,0)
   xaxis_wt = .50; if (present(xwt))  xaxis_wt = min(max(xwt,0.0_r8),1.0_r8)

!  do not allow more than one sponge layer
   if (numlev > 1) call error_mesg ('bgrid_sponge_mod',  &
                                    'numlev > 1 ', FATAL)

   do_init=.false.

!-----------------------------------------------------------------------

 end subroutine sponge_init

!#######################################################################

 subroutine local_filter_mass ( Hgrid, nplev, coeff_tmp, coeff_trs, dt, &
                                dp, t, tr, tdt, trdt )

   type (horiz_grid_type), intent(inout) :: Hgrid
   integer, intent(in)                   :: nplev
   real(r8),    intent(in)                   :: coeff_tmp, coeff_trs, dt
   real(r8),    intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: dp, t
   real(r8),    intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: tdt
   real(r8),    intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:,:) :: tr
   real(r8),    intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:,:) :: trdt

   real(r8), dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub,size(t,3)) :: tmp, akew, akns, akdp, ztmp
   real(r8), dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: akew2, akns2, tew, tns
   integer :: i, j, k, is, ie, js, je, jsd, jed, n
   real(r8)    :: xwt, ywt

   is  = Hgrid%Tmp%is;   ie  = Hgrid%Tmp%ie
   js  = Hgrid%Tmp%js;   je  = Hgrid%Tmp%je
   jsd = Hgrid%Tmp%jsd;  jed = Hgrid%Tmp%jed

   xwt = 0.125 * xaxis_wt
   ywt = 0.125 * (1.-xaxis_wt)

 ! 2D geometric constants
   do j = js-1, je
   do i = is-1, ie
      akew2(i,j) = xwt * (Hgrid%Tmp%area(i,j)+Hgrid%Tmp%area(i+1,j))
      akns2(i,j) = ywt * (Hgrid%Tmp%area(i,j)+Hgrid%Tmp%area(i,j+1))
   enddo
   enddo

 ! 3D mass weighted constants
   do k = 1, size(t,3)
      akew(is-1:ie,js-1:je,k) = akew2(is-1:ie,js-1:je)
      akns(is-1:ie,js-1:je,k) = akns2(is-1:ie,js-1:je)
      akdp(:,:,k) = Hgrid%Tmp%rarea(:,:)/dt
   enddo
   do k = nplev+1, size(t,3)
      do j = js-1, je
      do i = is-1, ie
         akew(i,j,k) = akew(i,j,k) * (dp(i,j,k)+dp(i+1,j,k))
         akns(i,j,k) = akns(i,j,k) * (dp(i,j,k)+dp(i,j+1,k))
         akdp(i,j,k) = akdp(i,j,k) / (2.*dp(i,j,k))
      enddo
      enddo
   enddo


 ! temperature

   if ( coeff_tmp > small ) then
      tmp = t + dt*tdt
    !!!!! remove zonal mean !!!!!
      call zsum (Zdomain_tmp, tmp(:,jsd:jed,:), ztmp(:,jsd:jed,:))
      tmp(is:ie,js:je,:) = tmp(is:ie,js:je,:) -ztmp(is:ie,js:je,:)/(1.0 * Hgrid%nlon)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call update_halo ( Hgrid, TEMP, tmp )
      do k = 1, size(t,3)
         do j = js-1, je
         do i = is-1, ie
            tew(i,j) = akew(i,j,k) * (tmp(i+1,j,k)-tmp(i,j,k))
            tns(i,j) = akns(i,j,k) * (tmp(i,j+1,k)-tmp(i,j,k))
         enddo
         enddo
         do j = js, je
         do i = is, ie
            tdt(i,j,k) = tdt(i,j,k) + &
                  coeff_tmp*(tew(i,j)-tew(i-1,j)+tns(i,j)-tns(i,j-1))*akdp(i,j,k)
         enddo
         enddo
      enddo
   endif

 ! tracers

   if ( coeff_trs > small ) then
      do n = 1, size(tr,4)
         tmp = tr(:,:,:,n) + dt*trdt(:,:,:,n)
       !!!!! remove zonal mean !!!!!
         call zsum (Zdomain_tmp, tmp(:,jsd:jed,:), ztmp(:,jsd:jed,:))
         tmp(is:ie,js:je,:) = tmp(is:ie,js:je,:) -ztmp(is:ie,js:je,:)/(1.0 * Hgrid%nlon)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call update_halo ( Hgrid, TEMP, tmp )
         do k = 1, size(tr,3)
            do j = js-1, je
            do i = is-1, ie
               tew(i,j) = akew(i,j,k) * (tmp(i+1,j,k)-tmp(i,j,k))
               tns(i,j) = akns(i,j,k) * (tmp(i,j+1,k)-tmp(i,j,k))
            enddo
            enddo
            do j = js, je
            do i = is, ie
               trdt(i,j,k,n) = trdt(i,j,k,n) +  &
                    coeff_trs*(tew(i,j)-tew(i-1,j)+tns(i,j)-tns(i,j-1))*akdp(i,j,k)
            enddo
            enddo
         enddo
      enddo
   endif

 end subroutine local_filter_mass

!#######################################################################

 subroutine local_filter_vel ( Hgrid, nplev, coeff, dt, dp, u, v, udt, vdt )

   type (horiz_grid_type), intent(inout) :: Hgrid
   integer, intent(in)                   :: nplev
   real(r8),    intent(in)                   :: coeff, dt
   real(r8),    intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: dp, u, v
   real(r8),    intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: udt, vdt

   real(r8), dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub,size(u,3)) :: akew, akns, akdp, uu, vv, ztmp
   real(r8), dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) ::  akew2, akns2, uew, vew, uns, vns
   integer :: i, j, k, is, ie, js, je, jsd, jed
   real(r8)    :: xwt, ywt

   is  = Hgrid%Tmp%is;   ie  = Hgrid%Tmp%ie
   js  = Hgrid%Vel%js;   je  = Hgrid%Vel%je
   jsd = Hgrid%Vel%jsd;  jed = Hgrid%Vel%jed

   xwt = 0.125 * coeff * xaxis_wt
   ywt = 0.125 * coeff * (1.-xaxis_wt)

 ! 2D geometric constants
   do j = js, je+1
   do i = is, ie+1
      akew2(i,j) = xwt * (Hgrid%Vel%area(i-1,j)+Hgrid%Vel%area(i,j))
      akns2(i,j) = ywt * (Hgrid%Vel%area(i,j-1)+Hgrid%Vel%area(i,j))
   enddo
   enddo

 ! 3D mass weighted constants
   do k = 1, size(u,3)
      akew(is:ie+1,js:je+1,k) = akew2(is:ie+1,js:je+1)
      akns(is:ie+1,js:je+1,k) = akns2(is:ie+1,js:je+1)
      akdp(:,:,k) = Hgrid%Vel%rarea(:,:)/dt
   enddo
   do k = nplev+1, size(u,3)
      do j = js, je+1
      do i = is, ie+1
         akew(i,j,k) = akew(i,j,k) * (dp(i,j,k)+dp(i-1,j,k))
         akns(i,j,k) = akns(i,j,k) * (dp(i,j,k)+dp(i,j-1,k))
         akdp(i,j,k) = akdp(i,j,k) / (2.*dp(i,j,k))
      enddo
      enddo
   enddo

   uu = u + dt*udt
   vv = v + dt*vdt
 !!!!! remove zonal mean from u comp !!!!!
   call zsum (Zdomain_vel, uu(:,jsd:jed,:), ztmp(:,jsd:jed,:))
   uu(is:ie,js:je,:) = uu(is:ie,js:je,:) -ztmp(is:ie,js:je,:)/(1.0 * Hgrid%nlon)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call update_halo ( Hgrid, UWND, uu )
   call update_halo ( Hgrid, VWND, vv )

   do k = 1, size(u,3)
      do j = js, je+1
      do i = is, ie+1
         uew(i,j) = akew(i,j,k) * (uu(i,j,k)-uu(i-1,j,k))
         vew(i,j) = akew(i,j,k) * (vv(i,j,k)-vv(i-1,j,k))
         uns(i,j) = akns(i,j,k) * (uu(i,j,k)-uu(i,j-1,k))
         vns(i,j) = akns(i,j,k) * (vv(i,j,k)-vv(i,j-1,k))
      enddo
      enddo
     !---- remove meridional averaging at poles ----
      call vel_flux_boundary (Hgrid, uns)
      call vel_flux_boundary (Hgrid, vns)

      do j = js, je
      do i = is, ie
         udt(i,j,k) = udt(i,j,k) + (uew(i+1,j)-uew(i,j)+uns(i,j+1)-uns(i,j))*akdp(i,j,k)
         vdt(i,j,k) = vdt(i,j,k) + (vew(i+1,j)-vew(i,j)+vns(i,j+1)-vns(i,j))*akdp(i,j,k)
      enddo
      enddo
   enddo

 end subroutine local_filter_vel

!#######################################################################
! initializes domain2d type for summation in zonal direction

! subroutine zsum_init ( Domain, Zonal )
! type(domain2d),  intent(in)   :: Domain
! type(zsum_type), intent(out)  :: Zonal
! integer :: isg, ieg, jsg, jeg, layout(2)

! create new domian2d type with 1d zonal decomposition
!  call mpp_get_global_domain  ( Domain, isg, ieg, jsg, jeg )
!  layout = (/1,mpp_npes()/)

! error check
!  if (jeg-jsg+1 < mpp_npes()) call error_mesg ('bgrid_sponge_mod', 'number &
!     &of global latitude rows is less than the number of processors', FATAL)
!
!  call mpp_define_domains ( (/isg,ieg,jsg,jeg/), layout, Zonal%Domain, &
!                            xflags=CYCLIC_GLOBAL_DOMAIN )
!  call mpp_get_data_domain ( Zonal%Domain, xsize=Zonal%isize, ysize=Zonal%jsize )

! end subroutine zsum_init

!#######################################################################

 subroutine zsum (Zonal, array, zarray )
 type(zsum_type), intent(in)   :: Zonal
 real(r8),            intent(in)  ::   array(:,:,:)
 real(r8),            intent(out) ::  zarray(:,:,:)
 real(r8), dimension(Zonal%isize,Zonal%jsize,size(array,3)) :: zdat
 real(r8)    :: c    
 integer :: j, k

!!!   call mpp_redistribute ( Domain, array, Zonal%Domain, zdat )
zdat(:, :, :) = array(2:Zonal%isize + 1, 2:Zonal%jsize + 1, :)

   do k = 1, size(array,3)
   do j = 1, Zonal%jsize
     c = sum(zdat(:,j,k))
     zdat(:,j,k) = c
   enddo
   enddo
!!!   call mpp_redistribute ( Zonal%Domain, zdat, Domain, zarray )
zarray(2:Zonal%isize +1, 2:Zonal%jsize +1, :) = zdat(:, :, :)

 end subroutine zsum

!#######################################################################

end module bgrid_sponge_mod

