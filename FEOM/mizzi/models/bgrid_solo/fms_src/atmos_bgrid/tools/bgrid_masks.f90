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

module bgrid_masks_mod

!-----------------------------------------------------------------------
!
!    allocates storage and initializes masks and vertical indexing
!    associated with the step mountain vertical coordinate.
!
!-----------------------------------------------------------------------

use types_mod, only : r8
use bgrid_horiz_mod, only:  horiz_grid_type
use  bgrid_halo_mod, only:  update_halo, NOPOLE, UWND
use  bgrid_vert_mod, only:  vert_grid_type
use         fms_mod, only:  write_version_number
!!!use         mpp_mod, only:  mpp_max

implicit none
private

!-----------------------------------------------------------------------
!    ---- public data types ----

   public      mask_type
   public grid_mask_type

   type mask_type
      real(r8),    pointer, dimension(:,:,:) :: mask
      integer, pointer, dimension(:,:)   :: kbot
      integer                            :: kbotmin
   end type mask_type

!  mask  = step-mountain topography mask (0.0 or 1.0) for
!            mass (height) grid points
!  kbot  = lowest model level above ground
!
!   note:  for the sigma coordinate, mask = 1.0 everywhere, and
!          kbot = number of vertical model levels

   type grid_mask_type
      type(mask_type) :: Tmp, Vel
      logical :: sigma
   end type grid_mask_type

!  Tmp = grid masking values for the temperature/mass grid
!  Vel = grid masking values for the velocity/momentum grid
!  sigma = logical flag that specific whether vertical coordinate is
!            the step-mountain (eta) or sigma coordinate
!
!-----------------------------------------------------------------------
!    ---- public interfaces (and some that may become public) ----

   public  grid_masks_init

   private compute_mass_mask, compute_vel_mask, compute_lowest_level

   logical :: sigma  ! local variable

!  version number info

   character(len=128) :: version='$Revision$'
   character(len=128) :: tagname='$Id$'
   logical :: do_log = .true.

!-----------------------------------------------------------------------

contains

!#######################################################################

   function grid_masks_init (Hgrid, Vgrid, res) result (Mask)

!-----------------------------------------------------------------------
!
!     arguments (intent in)
!
!     res  =  reciprical of eta at the surface
!
!-----------------------------------------------------------------------

   type (horiz_grid_type), intent(inout) :: Hgrid
   type  (vert_grid_type), intent(in)    :: Vgrid
   real(r8),                   intent(in)    :: res(Hgrid%ilb:,Hgrid%jlb:)

   type (grid_mask_type) :: Mask

!-----------------------------------------------------------------------
!   sigma = logical flag for terrian-following versus step-mountain
!           vertical coordinate.

   logical :: sigma
   integer :: kx
   real(r8)    :: maxres

   if (do_log) then
      call write_version_number (version,tagname)
      do_log = .false.
   endif

   kx = Vgrid % nlev

!--------------allocate global storage----------------------------------

 allocate ( Mask%Tmp%mask (Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub, kx), &
            Mask%Vel%mask (Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub, kx), &
            Mask%Tmp%kbot (Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub),     &
            Mask%Vel%kbot (Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub)      )

!--------------is this an eta of sigma coordinate mountain? ------------

      maxres = maxval(res)
      !!!call mpp_max (maxres)

      if (maxres > 1.0001) then
         Mask % sigma=.false.
         sigma=.false.
            if (Vgrid%hybrid) then
               write (*,100) 'eta/hybrid'
            else
               write (*,100) 'eta'
            endif
             ! print gaudy banner for untested eta option
               write (*,'(4(a/))') '*******************************************************', &
                                   'WARNING: The eta coordinate option has not been tested.', &
                                   '         Proceed with caution.',                          &
                                   '*******************************************************'
      else
         Mask % sigma=.true.
         sigma=.true.
            if (Vgrid%hybrid) then
               write (*,100) 'sigma/hybrid'
            else
               write (*,100) 'sigma'
            endif
      endif

 100  format (/,'B-grid dynamical core has been initialized with the ',a,' vertical coordinate.')

!--------------topography masks ----------------------------------------

      Mask % Tmp % mask = compute_mass_mask (res, Vgrid % aeta)
      Mask % Vel % mask = compute_vel_mask  (res, Vgrid % aeta)
      call update_halo (Hgrid, UWND, Mask%Vel%mask, flags=NOPOLE)
!!!!! call update_halo (Hgrid, UWND, Mask%Vel%mask)   ! sets mask=0 at poles

!------------- compute the lowest model level --------------------------

      Mask % Tmp % kbot = compute_lowest_level (Mask % Tmp % mask)
      Mask % Vel % kbot = compute_lowest_level (Mask % Vel % mask)

!     ----- global values -----

      Mask % Tmp % kbotmin = minval(Mask % Tmp % kbot)
      Mask % Vel % kbotmin = minval(Mask % Vel % kbot)

!-----------------------------------------------------------------------

end function grid_masks_init

!#######################################################################

   function compute_mass_mask (res, aeta) result (mask)

   real(r8), intent(in) :: res(:,:), aeta(:)
   real(r8), dimension(size(res,1),size(res,2),size(aeta)) :: mask
   integer  i, j, k

      mask = 1.0

      if (.not.sigma) then
         do j=1,size(res,2); do i=1,size(res,1)
         do k=1,size(aeta)
            if (aeta(k) > (1.0/res(i,j))) mask(i,j,k) = 0.0
         enddo; enddo; enddo
      endif

   end function compute_mass_mask

!#######################################################################

   function compute_vel_mask (res, aeta) result (mask)

   real(r8), intent(in) :: res(:,:), aeta(:)
   real(r8), dimension(size(res,1),size(res,2),size(aeta)) :: mask
   integer  i, j, k

      mask = 1.0

      if (.not.sigma) then
         do j=2,size(res,2); do i=2,size(res,1)
         do k=1,size(aeta)
            if (aeta(k) > (1.0/res(i,j))) then
                mask(i-1,j-1,k) = 0.0
                mask(i  ,j-1,k) = 0.0
                mask(i-1,j  ,k) = 0.0
                mask(i  ,j  ,k) = 0.0
            endif
         enddo; enddo; enddo
      endif

   end function compute_vel_mask

!#######################################################################

   function compute_lowest_level (mask) result (kbot)

   real(r8), intent(in) :: mask(:,:,:)
   integer, dimension(size(mask,1),size(mask,2)) :: kbot
   integer   i, j, k, kdim

      kdim = size(mask,3)
      kbot = kdim

      if (.not.sigma) then
         do j=1,size(mask,2); do i=1,size(mask,1)
         do k=1,kdim
            if (mask(i,j,k) < 0.50) then
               kbot(i,j)=k-1; exit
            endif
         enddo; enddo; enddo
       ! must not be zero
         kbot = max(kbot,1)
      endif

   end function compute_lowest_level

!#######################################################################

end module bgrid_masks_mod

