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

module bgrid_change_grid_mod

!-----------------------------------------------------------------------
!   utility programs specific to the b-grid configuration.
!-----------------------------------------------------------------------

use types_mod, only : r8
use bgrid_horiz_mod, only:  horiz_grid_type

implicit none
private

!-----------------------------------------------------------------------

  public  mass_to_vel, vel_to_mass

interface mass_to_vel
  module procedure mass_to_vel_2d, mass_to_vel_area_2d,  &
                   mass_to_vel_3d, mass_to_vel_area_3d
end interface

interface vel_to_mass
  module procedure vel_to_mass_2d, vel_to_mass_area_2d,  &
                   vel_to_mass_3d, vel_to_mass_area_3d
end interface

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine mass_to_vel_3d (fm, fv)

!-----------------------------------------------------------------------
!
!        converts data from mass points to velocity points
!        by doing a 4-pt average.
!
!    input:   fm    2-d or 3-d array located at mass points
!
!    output:  fv    2-d or 3-d array located at velocity points
!
!    note:  fm and fv CAN NOT be the same array
!

      real(r8), intent(in),  dimension(:,:,:) :: fm
      real(r8), intent(out), dimension(:,:,:) :: fv

!-----------------------------------------------------------------------
      real(r8), dimension(size(fm,1),size(fm,2)) :: wk
      integer  ::  i, j, k, ie, je
!-----------------------------------------------------------------------

      ie = size(fm,1);  je = size(fm,2)

      do k = 1, size(fm,3)
           wk(ie,:) = fm(ie,:,k)
           wk(:,je) = fm(:,je,k)
      do j = 1, je-1
      do i = 1, ie-1
           wk(i,j) = (fm(i,j  ,k) + fm(i+1,j  ,k) +  &
                      fm(i,j+1,k) + fm(i+1,j+1,k)) * 0.25
      enddo
      enddo
           fv(:,:,k) = wk(:,:)
      enddo

!-----------------------------------------------------------------------

end subroutine mass_to_vel_3d

!#######################################################################

subroutine mass_to_vel_2d (fm,fv)

!-----------------------------------------------------------------------
!
!        converts data from mass points to velocity points
!        by doing a 4-pt average.
!
!    input:   fm    2-d or 3-d array located at mass points
!                     dimension ix x jx (x nl)
!
!    output:  fv    2-d or 3-d array located at velocity points
!                     dimension ix x jx (x nl)
!

          real(r8), intent(in),  dimension(:,:) :: fm
          real(r8), intent(out), dimension(:,:) :: fv

!-----------------------------------------------------------------------
      real(r8), dimension(size(fm,1),size(fm,2),1) :: fm3,fv3
!-----------------------------------------------------------------------

      fm3(:,:,1) = fm(:,:)
      call mass_to_vel_3d (fm3,fv3)
      fv(:,:) = fv3(:,:,1)

!-----------------------------------------------------------------------

end subroutine mass_to_vel_2d

!#######################################################################

subroutine vel_to_mass_3d (u,v,um,vm,mask)

!-----------------------------------------------------------------------
!
!        converts data from velocity points to mass points
!        by doing a 4-pt average.
!
!    input:  u,v    2-d or 3-d arrays of velocity components array 
!                     dimensioned ix x jx (x nl)
!            mask   2-d or 3-d topography mask array located at 
!                     velocity points dimensioned ix x jx (x nl)
!
!    output: um,vm  2-d or 3-d array averaged to mass points
!                     dimensioned ix x jx (x nl)
!
!    note:  u,um and v,vm CAN be the same array
!

        real(r8), intent(in),  dimension(:,:,:) :: u,v
        real(r8), intent(out), dimension(:,:,:) :: um,vm
        real(r8), intent(in),  dimension(:,:,:) :: mask

!-----------------------------------------------------------------------
        real(r8), dimension(size(u,1),size(u,2)) :: u1, v1, wt
        integer :: i, j, k, ie, je
!-----------------------------------------------------------------------

     ie = size(u,1);  je = size(u,2)

     do k = 1, size(u,3)

         u1(1,:) = u(1,:,k)
         u1(:,1) = u(:,1,k)
         v1(1,:) = v(1,:,k)
         v1(:,1) = v(:,1,k)

     do j = 2, je
     do i = 2, ie
         wt(i,j) = mask(i-1,j-1,k) + mask(i,j-1,k) +  &
                   mask(i-1,j  ,k) + mask(i,j  ,k)
     enddo
     enddo

     where (wt(2:ie,2:je) > 0.50)
         wt(2:ie,2:je) = 1./wt(2:ie,2:je)
     elsewhere
         wt(2:ie,2:je) = 0.
     endwhere

     do j = 2, je
     do i = 2, ie
        u1(i,j) = (u(i-1,j-1,k)+u(i,j-1,k)+  &
                   u(i-1,j  ,k)+u(i,j  ,k))*wt(i,j)
        v1(i,j) = (v(i-1,j-1,k)+v(i,j-1,k)+  &
                   v(i-1,j  ,k)+v(i,j  ,k))*wt(i,j)
     enddo
     enddo

         um(:,:,k) = u1
         vm(:,:,k) = v1

     enddo

!-----------------------------------------------------------------------

end subroutine vel_to_mass_3d

!#######################################################################

subroutine vel_to_mass_2d (u,v,um,vm,mask)

!-----------------------------------------------------------------------
!
!        converts data from velocity points to mass points
!        by doing a 4-pt average.
!
!    input:  u,v    2-d or 3-d arrays of velocity components at
!                     velocity points
!            mask   2-d or 3-d topography mask array located at 
!                     velocity points (optional)
!
!    output: um,vm  2-d or 3-d array averaged to mass points
!
!    note:  u,um and v,vm CAN be the same array
!

        real(r8), intent(in),  dimension(:,:) :: u,v
        real(r8), intent(out), dimension(:,:) :: um,vm
        real(r8), intent(in),  dimension(:,:), optional :: mask

!-----------------------------------------------------------------------
      real(r8), dimension(size(u,1),size(u,2),1) :: u3,v3,um3,vm3,mask3
!-----------------------------------------------------------------------

      u3(:,:,1) = u(:,:);  v3(:,:,1) = v(:,:)

      if (present(mask)) then
         mask3(:,:,1) = mask(:,:)
      else
         mask3(:,:,1) = 1.0
      endif

         call vel_to_mass_3d (u3,v3, um3,vm3, mask3)

      um(:,:) = um3(:,:,1);  vm(:,:) = vm3(:,:,1)

!-----------------------------------------------------------------------

end subroutine vel_to_mass_2d

!#######################################################################

subroutine mass_to_vel_area_3d (Hgrid, fm, fv)

!-----------------------------------------------------------------------
!
!        converts data from mass points to velocity points
!        by doing a 4-pt area-weighted average.
!
!    input:   fm    2-d or 3-d array located at mass points
!
!    output:  fv    2-d or 3-d array located at velocity points
!
!    note:  fm and fv CAN NOT be the same array
!

    type(horiz_grid_type), intent(in) :: Hgrid
    real(r8), intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: fm
    real(r8), intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: fv

!-----------------------------------------------------------------------
    real(r8), dimension(Hgrid%ilb:Hgrid%iub, &
                    Hgrid%jlb:Hgrid%jub) :: wk
    integer :: i, j, k, ilb, iub, jlb, jub
!-----------------------------------------------------------------------

    ilb = Hgrid%ilb;  iub = Hgrid%iub
    jlb = Hgrid%jlb;  jub = Hgrid%jub

    do k = 1, size(fm,3)

       wk( : , jub) = fm( : , jub, k)
       wk(iub,  : ) = fm(iub,  : , k)

!      --- interpolate mass from pole to pole ---
       do j = jlb, jub-1
       do i = ilb, iub-1
            wk(i,j) = (Hgrid%Tmp%area(i  ,j  )*fm(i  ,j  ,k)  &
                     + Hgrid%Tmp%area(i+1,j  )*fm(i+1,j  ,k)  &
                     + Hgrid%Tmp%area(i  ,j+1)*fm(i  ,j+1,k)  &
                     + Hgrid%Tmp%area(i+1,j+1)*fm(i+1,j+1,k)) &
                     * Hgrid%Vel%rarea(i,j)*0.25
       enddo
       enddo




         fv(:,:,k) = wk(:,:)

   end do


!-----------------------------------------------------------------------

end subroutine mass_to_vel_area_3d

!#######################################################################

subroutine mass_to_vel_area_2d (Hgrid, fm, fv)

!-----------------------------------------------------------------------
!
!        converts data from mass points to velocity points
!        by doing a 4-pt area-weighted average.
!
!    input:   fm    2-d or 3-d array located at mass points
!                     dimension ix x jx (x nl)
!
!    output:  fv    2-d or 3-d array located at velocity points
!                     dimension ix x jx (x nl)
!

    type(horiz_grid_type), intent(in) :: Hgrid
    real(r8), intent(in),  dimension(:,:) :: fm
    real(r8), intent(out), dimension(:,:) :: fv

!-----------------------------------------------------------------------
    real(r8), dimension(size(fm,1),size(fm,2),1) :: fm3, fv3
!-----------------------------------------------------------------------

    fm3(:,:,1) = fm(:,:)
    call mass_to_vel_area_3d (Hgrid, fm3, fv3)
    fv(:,:) = fv3(:,:,1)

!-----------------------------------------------------------------------

end subroutine mass_to_vel_area_2d

!#######################################################################

subroutine vel_to_mass_area_3d (Hgrid, u, v, um, vm, mask)

!-----------------------------------------------------------------------
!
!        converts data from velocity points to mass points
!        by doing a 4-pt area-weighted average.
!
!    input:  u,v    2-d or 3-d arrays of velocity components array 
!                     dimensioned ix x jx (x nl)
!            mask   2-d or 3-d topography mask array located at 
!                     velocity points dimensioned ix x jx (x nl)
!
!    output: um,vm  2-d or 3-d array averaged to mass points
!                     dimensioned ix x jx (x nl)
!
!    note:  u,um and v,vm CAN be the same array
!

    type(horiz_grid_type), intent(in)   :: Hgrid
    real(r8), intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: u, v
    real(r8), intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: um, vm
    real(r8), intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: mask

!-----------------------------------------------------------------------
    real(r8), dimension(Hgrid%ilb:Hgrid%iub, &
                    Hgrid%jlb:Hgrid%jub) :: wt, u1, v1
    integer :: i, j, k, ilb, iub, jlb, jub
!-----------------------------------------------------------------------

    ilb = lbound(u,1);  iub = ubound(u,1)
    jlb = lbound(u,2);  jub = ubound(u,2)

    do k = 1, size(u,3)

      u1(ilb,:) = u(ilb,:,k)
      u1(:,jlb) = u(:,jlb,k)
      v1(ilb,:) = v(ilb,:,k)
      v1(:,jlb) = v(:,jlb,k)

      do j = jlb+1, jub
      do i = ilb+1, iub
         wt(i,j) = mask(i-1,j-1,k) + mask(i  ,j-1,k) &
                 + mask(i-1,j  ,k) + mask(i  ,j  ,k)
      enddo
      enddo

      where (wt(ilb+1:iub,jlb+1:jub) > 0.50)
         wt(ilb+1:iub,jlb+1:jub) = Hgrid%Tmp%rarea(ilb+1:iub,jlb+1:jub) / &
                                             wt(ilb+1:iub,jlb+1:jub)
      elsewhere
         wt(ilb+1:iub,jlb+1:jub) = 0.
      endwhere

      do j = jlb+1, jub
      do i = ilb+1, iub
        u1(i,j) = (Hgrid%Vel%area(i-1,j-1)*u(i-1,j-1,k) &
                 + Hgrid%Vel%area(i  ,j-1)*u(i  ,j-1,k) &
                 + Hgrid%Vel%area(i-1,j  )*u(i-1,j  ,k) &
                 + Hgrid%Vel%area(i  ,j  )*u(i  ,j  ,k)) * wt(i,j)
        v1(i,j) = (Hgrid%Vel%area(i-1,j-1)*v(i-1,j-1,k) &
                 + Hgrid%Vel%area(i  ,j-1)*v(i  ,j-1,k) &
                 + Hgrid%Vel%area(i-1,j  )*v(i-1,j  ,k) &
                 + Hgrid%Vel%area(i  ,j  )*v(i  ,j  ,k)) * wt(i,j)
      enddo
      enddo

        um (:,:,k) = u1 (:,:)
        vm (:,:,k) = v1 (:,:)

    enddo

!-----------------------------------------------------------------------

end subroutine vel_to_mass_area_3d

!#######################################################################

subroutine vel_to_mass_area_2d (Hgrid, u, v, um, vm, mask)

!-----------------------------------------------------------------------
!
!        converts data from velocity points to mass points
!        by doing a 4-pt average.
!
!    input:  u,v    2-d or 3-d arrays of velocity components at
!                     velocity points
!            mask   2-d or 3-d topography mask array located at 
!                     velocity points (optional)
!
!    output: um,vm  2-d or 3-d array averaged to mass points
!
!    note:  u,um and v,vm CAN NOT be the same array
!

    type(horiz_grid_type), intent(in)   :: Hgrid
    real(r8), intent(in),  dimension(:,:) :: u, v
    real(r8), intent(out), dimension(:,:) :: um, vm
    real(r8), intent(in),  dimension(:,:), optional :: mask

!-----------------------------------------------------------------------
    real(r8), dimension(size(u,1),size(u,2),1) :: u3, v3, um3, vm3, mask3
!-----------------------------------------------------------------------

      u3(:,:,1) = u(:,:);  v3(:,:,1) = v(:,:)

      if (present(mask)) then
         mask3(:,:,1) = mask(:,:)
      else
         mask3(:,:,1) = 1.0
      endif

      call vel_to_mass_area_3d (Hgrid, u3, v3, um3, vm3, mask3)

      um(:,:) = um3(:,:,1);  vm(:,:) = vm3(:,:,1)

!-----------------------------------------------------------------------

end subroutine vel_to_mass_area_2d

!#######################################################################

end module bgrid_change_grid_mod

