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

module vert_advection_mod

!-------------------------------------------------------------------------------

use types_mod, only : r8
use fms_mod, only: error_mesg, FATAL

implicit none
private

public :: vert_advection
integer, parameter, public :: SECOND_CENTERED = 101, FOURTH_CENTERED = 102
integer, parameter, public :: VAN_LEER_LINEAR = 103
integer, parameter, public :: FLUX_FORM = 201, ADVECTIVE_FORM = 202

interface vert_advection
   module procedure vert_advection_1d, vert_advection_2d, vert_advection_3d
end interface

contains

!-------------------------------------------------------------------------------

 subroutine vert_advection_3d ( dt, w, dz, r, rdt, mask, scheme, form )

 real(r8), intent(in)                    :: dt
 real(r8), intent(in),  dimension(:,:,:) :: w, dz, r
 real(r8), intent(out), dimension(:,:,:) :: rdt
 real(r8),    intent(in), optional :: mask(:,:,:)
 integer, intent(in), optional :: scheme, form

! INPUT
!   dt  = time step in seconds
!   w   = advecting velocity at the vertical boundaries of the grid boxes
!         does not assume velocities at top and bottom are zero
!         units = [units of dz / second]
!   dz  = depth of model layers in arbitrary units (usually pressure)
!   r   = advected quantity in arbitrary units
!
! OUTPUT
!   rdt = advective tendency for quantity "r" weighted by depth of layer
!         units = [units of r * units of dz / second]
!
! OPTIONAL INPUT
!   mask   = mask for below ground layers,
!            where mask > 0 for layers above ground
!   scheme = differencing scheme, use one of these values:
!               SECOND_CENTERED = second-order centered
!               FOURTH_CENTERED = fourth-order centered
!               VAN_LEER_LINEAR = piecewise linear, finite volume (van Leer)
!   form   = form of equations, use one of these values:
!               FLUX_FORM      = solves for -d(wr)/dt
!               ADVECTIVE_FORM = solves for -w*d(r)/dt
!
! NOTE
!   size(w,3) == size(dz,3)+1 == size(r,3)+1 == size(rdt,3)+1 == size(mask,3)+1

 real(r8), dimension(size(r,1),size(r,2),size(r,3)) :: slp
 real(r8), dimension(size(w,1),size(w,2),size(w,3)) :: flux, rst
 real(r8)    :: c1, c2
 real(r8)    :: small = 1.e-6
 integer :: k, ks, ke
 integer :: diff_scheme, eqn_form

 ! set default values for optional arguments
   diff_scheme = VAN_LEER_LINEAR
   eqn_form    = FLUX_FORM
   if (present(scheme)) diff_scheme = scheme
   if (present(form))   eqn_form    = form

 ! note: size(r,3)+1 = size(w,3)
   if (size(w,3) /= size(r,3)+1) &
      call error_handler ('vertical dimension of input arrays inconsistent')

   ks = 1; ke = size(r,3)

 ! determine value of R at flux boundaries

   rst(:,:,ks)   = r(:,:,ks)   ! most likely w = 0 at these points
   rst(:,:,ke+1) = r(:,:,ke)

   select case (diff_scheme)
      case (SECOND_CENTERED)
         do k = ks+1, ke
            rst(:,:,k) = 0.5*(r(:,:,k)+r(:,:,k-1))
         enddo

      case (FOURTH_CENTERED)
         c1 = 7./12.;  c2 = 1./12.
         if (present(mask)) then
          ! second order if adjacent to ground
            do k = ks+2, ke-1
               where (mask(:,:,k+1) > small)
                  rst(:,:,k) = c1*(r(:,:,k)+r(:,:,k-1)) - c2*(r(:,:,k+1)+r(:,:,k-2))
               elsewhere
                  rst(:,:,k) = 0.5*(r(:,:,k)+r(:,:,k-1))
               endwhere
            enddo
         else
          ! no mask: always fourth order
            do k = ks+2, ke-1
               rst(:,:,k) = c1*(r(:,:,k)+r(:,:,k-1)) - c2*(r(:,:,k+1)+r(:,:,k-2))
            enddo
         endif
         ! second order at top and bottom
         rst(:,:,ks+1) = 0.5*(r(:,:,ks+1)+r(:,:,ks  ))
         rst(:,:,ke)   = 0.5*(r(:,:,ke  )+r(:,:,ke-1))

      case (VAN_LEER_LINEAR)
       ! slope along the z-axis
         call slope_z ( r, dz, slp )
         do k = ks+1, ke
            where (w(:,:,k) >= 0.)
               rst(:,:,k) = (r(:,:,k-1) + 0.5*slp(:,:,k-1)*(1.-dt*w(:,:,k)/dz(:,:,k-1)))
            elsewhere
               rst(:,:,k) = (r(:,:,k  ) - 0.5*slp(:,:,k  )*(1.+dt*w(:,:,k)/dz(:,:,k  )))
            endwhere
         enddo

      case default
        ! ERROR
          call error_handler ('invalid value for optional argument scheme')
   end select

 ! compute fluxes
   flux(:,:,ks:ke+1) = w(:,:,ks:ke+1) * rst(:,:,ks:ke+1)

 ! vertical advective tendency

   select case (eqn_form)
      case (FLUX_FORM)
         do k = ks, ke
            rdt (:,:,k) = - (flux(:,:,k+1) - flux (:,:,k)) 
!del        rdt (:,:,k) = - (flux(:,:,k+1) - flux (:,:,k)) / dz(:,:,k)
         enddo
      case (ADVECTIVE_FORM)
         do k = ks, ke
            rdt (:,:,k) = - (flux(:,:,k+1) - flux (:,:,k) - &
                             r(:,:,k)*(w(:,:,k+1)-w(:,:,k)))
!del                         r(:,:,k)*(w(:,:,k+1)-w(:,:,k))) / dz(:,:,k)
         enddo
      case default
        ! ERROR
          call error_handler ('invalid value for optional argument form')
   end select


 end subroutine vert_advection_3d

!-------------------------------------------------------------------------------

 subroutine slope_z ( r, dz, slope )
 real(r8), intent(in),  dimension(:,:,:) :: r, dz
 real(r8), intent(out), dimension(:,:,:) :: slope

 real(r8)    :: grad(size(r,1),size(r,2),2:size(r,3))
 real(r8)    :: rmin, rmax
 integer :: i, j, k, n

  n = size(r,3)

! compute slope (weighted for unequal levels)
  do k = 2, n
    grad(:,:,k) = (r(:,:,k)-r(:,:,k-1))/(dz(:,:,k)+dz(:,:,k-1))
  enddo
  do k = 2, n-1
    slope(:,:,k) = (grad(:,:,k+1)+grad(:,:,k))*dz(:,:,k)
  enddo
    slope(:,:,1) = 2.*grad(:,:,2)*dz(:,:,1)
    slope(:,:,n) = 2.*grad(:,:,n)*dz(:,:,n)

! apply limiters to slope
  do k = 1, n
  do j = 1, size(r,2)
  do i = 1, size(r,1)
    if (k >= 2 .and. k <= n-1) then
      rmin = min(r(i,j,k-1), r(i,j,k), r(i,j,k+1))
      rmax = max(r(i,j,k-1), r(i,j,k), r(i,j,k+1))
    else if (k == 1) then
      rmin = min(r(i,j,k), r(i,j,k+1))
      rmax = max(r(i,j,k), r(i,j,k+1))
    else if (k == n) then
      rmin = min(r(i,j,k-1), r(i,j,k))
      rmax = max(r(i,j,k-1), r(i,j,k))
    endif
    slope(i,j,k) = sign(1.0_r8,slope(i,j,k)) *  &
             min( abs(slope(i,j,k)), 2.*(r(i,j,k)-rmin), 2.*(rmax-r(i,j,k)) )
  enddo
  enddo
  enddo

 end subroutine slope_z

!-------------------------------------------------------------------------------
!--------------------------- overloaded versions -------------------------------

 subroutine vert_advection_1d ( dt, w, dz, r, rdt, mask, scheme, form )
 
 real(r8), intent(in)                :: dt
 real(r8), intent(in),  dimension(:) :: w, dz, r
 real(r8), intent(out), dimension(:) :: rdt
 real(r8),    intent(in), optional :: mask(:)
 integer, intent(in), optional :: scheme, form

 real(r8), dimension(1,1,size(r,1)) :: dz3, r3, rdt3, mask3
 real(r8), dimension(1,1,size(w,1)) :: w3

  ! input
    w3 (1,1,:) = w
    dz3(1,1,:) = dz
    r3 (1,1,:) = r

    if (present(mask)) then
       mask3(1,1,:) = mask
       call vert_advection_3d ( dt, w3, dz3, r3, rdt3, mask=mask3, scheme=scheme, form=form )
    else
       call vert_advection_3d ( dt, w3, dz3, r3, rdt3, scheme=scheme, form=form )
    endif

  ! output
    rdt = rdt3(1,1,:)

 end subroutine vert_advection_1d

!-------------------------------------------------------------------------------

 subroutine vert_advection_2d ( dt, w, dz, r, rdt, mask, scheme, form )

 real(r8), intent(in)                  :: dt
 real(r8), intent(in),  dimension(:,:) :: w, dz, r
 real(r8), intent(out), dimension(:,:) :: rdt
 real(r8),    intent(in), optional :: mask(:,:)
 integer, intent(in), optional :: scheme, form

 real(r8), dimension(size(r,1),1,size(r,2)) :: dz3, r3, rdt3, mask3
 real(r8), dimension(size(w,1),1,size(w,2)) :: w3

  ! input
    w3 (:,1,:) = w
    dz3(:,1,:) = dz
    r3 (:,1,:) = r

    if (present(mask)) then
       mask3(:,1,:) = mask
       call vert_advection_3d ( dt, w3, dz3, r3, rdt3, mask=mask3, scheme=scheme, form=form )
    else
       call vert_advection_3d ( dt, w3, dz3, r3, rdt3, scheme=scheme, form=form )
    endif

  ! output
    rdt = rdt3(:,1,:)

 end subroutine vert_advection_2d

!-------------------------------------------------------------------------------

 subroutine error_handler ( message )
 character(len=*), intent(in) :: message

   call error_mesg ('vert_advection', trim(message), FATAL)

!  print *, 'FATAL ERROR in vert_advection'
!  print *, trim(message)
!  stop 111

 end subroutine error_handler

!-------------------------------------------------------------------------------

end module vert_advection_mod

