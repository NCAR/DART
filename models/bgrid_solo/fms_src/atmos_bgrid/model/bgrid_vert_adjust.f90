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

module bgrid_vert_adjust_mod

!-----------------------------------------------------------------------

use bgrid_vert_mod, only:  vert_grid_type
use  constants_mod, only:  CP

implicit none
private

!-----------------------------------------------------------------------

real, parameter :: ef4t = 1./CP

!-----------------------------------------------------------------------
!------------- interfaces -------------

public   vert_adjust

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine vert_adjust (res, div, wta, wtb, mask, Vgrid,  &
                        omgalf, etadot, psdt)

   real, intent(in), dimension(:,:)   :: res
   real, intent(in), dimension(:,:,:) :: div, wta, wtb, mask
type(vert_grid_type), intent(in) :: Vgrid
   real, intent(inout) :: omgalf(:,:,:)
   real, intent(inout) :: etadot(:,:,:)
   real, intent(out)   :: psdt(:,:)

!-----------------------------------------------------------------------
   real, dimension(size(div,1),size(div,2),size(div,3)) ::  sdiv
!-----------------------------------------------------------------------

!------------------ vertical Adjustments -------------------------------

   call ps_tendency (div, psdt, sdiv)
   call vert_omgalf (Vgrid, res, wta, wtb, sdiv, mask, omgalf)

!----- compute "eta-dot", the vertical velocity -----

   call vert_velocity (res, sdiv, mask, Vgrid%eta, etadot)

!-----------------------------------------------------------------------

end subroutine vert_adjust

!#######################################################################

subroutine ps_tendency (div, psdt, sdiv)

!-----------------------------------------------------------------------
!
!          Routine for calculating the vertical integral
!          of divergence and surface pressure tendency
!
!-----------------------------------------------------------------------

   real, intent(in)  :: div(:,:,:)
   real, intent(out) :: psdt(:,:), sdiv(:,:,:)

   integer  k, kdim
!-----------------------------------------------------------------------

      kdim = size(div,3)

!---------integrate divergence to get surface pressure tendency---------

      sdiv(:,:,1) = div(:,:,1)
   do k = 2, kdim
      sdiv(:,:,k) = sdiv(:,:,k-1) + div(:,:,k)
   enddo
      psdt(:,:) = -sdiv(:,:,kdim)

!-----------------------------------------------------------------------

end subroutine ps_tendency

!#######################################################################

subroutine vert_velocity (res, sdiv, mask, eta, etadot)

!-----------------------------------------------------------------------
!
!    Routine for calculating the vertical velocity ("eta-dot")
!
!-----------------------------------------------------------------------

   real, intent(in)  :: res(:,:), sdiv(:,:,:),  &
                        mask(:,:,:), eta(:)
   real, intent(inout) :: etadot(:,:,:)

   real, dimension(size(res,1),size(res,2)) :: pret
   integer  k, kdim

!---- computation of etadot (add onto previous value) ----

   kdim = size(sdiv,3)

   pret(:,:) = -sdiv(:,:,kdim)*res(:,:)

   do k = 2, kdim
      etadot(:,:,k) = etadot(:,:,k)  &
                      - (pret(:,:)*eta(k)+sdiv(:,:,k-1))*mask(:,:,k)
   enddo

!-----------------------------------------------------------------------

end subroutine vert_velocity

!#######################################################################

subroutine vert_omgalf (Vgrid, res, wta, wtb, sdiv, mask, omgalf)

!-----------------------------------------------------------------------
!
!   Routine for calculating the vertical part of the omega-alpha term
!
!-----------------------------------------------------------------------

   type(vert_grid_type), intent(in)       :: Vgrid
   real, intent(in),    dimension(:,:)    :: res
   real, intent(in),    dimension(:,:,:)  :: wta, wtb, sdiv, mask
   real, intent(inout), dimension(:,:,:)  :: omgalf

   real, dimension(size(res,1),size(res,2)) :: pret
   integer  k, kdim

!-----------------------------------------------------------------------
!-----kinetic energy generation terms in the thermodynamic equation-----

   kdim = size(sdiv,3)

      omgalf(:,:,1) = omgalf(:,:,1) - wtb(:,:,1)*sdiv(:,:,1)*mask(:,:,1)

!!!do k=2,kdim-1
   do k=2,kdim
      omgalf(:,:,k) = omgalf(:,:,k) - mask(:,:,k)* &
                    (wta(:,:,k)*sdiv(:,:,k-1)+wtb(:,:,k)*sdiv(:,:,k))
   enddo

!  --- may need to add/modify for eta coordinate ---
!  if (kdim > 1) then
!     pret(:,:) = -sdiv(:,:,kdim)*res(:,:)
!     omgalf(:,:,kdim) = omgalf(:,:,kdim) + mask(:,:,kdim)* &
!                   (wtb(:,:,k)*pret(:,:)-wta(:,:,k)*sdiv(:,:,kdim-1))
!  endif

!-----------------------------------------------------------------------

end subroutine vert_omgalf

!#######################################################################

end module bgrid_vert_adjust_mod

