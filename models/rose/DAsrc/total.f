      subroutine total
!------------------------------------------------------------------------
!
!    calculates ozone and molecular oxygen columns [cm-2]
!
!    ozone column specified at the top of the model using climatology
!    vertical integration at all other levels
!------------------------------------------------------------------------

      use params

      use dynam, only : h, dz
      use phys, only : cl_o3
      use chem, only : hnm, o3t, o2t, qn1

      implicit none

      real :: tt
      integer :: i, j, k
      real, dimension(nz) :: o2prof, o3prof

      real, save :: dzcm, hcm
      logical, save :: entered

      data entered /.false./

!... initialize layer depth and scale height

      if( .not. entered ) then
	 entered = .true.

         dzcm = 1.e2 * dz      ! depth of layer
         hcm = 1.e2 * h        ! scale height in cm

      end if

!... loop over all grid points to calculate columns

      do j=1,ny

         tt = cl_o3(nz,j) * hcm

         do i = 1,nx

	    o3prof = qn1(:,i,j,19) * hnm(:,i,j)
	    o2prof = qn1(:,i,j,26) * hnm(:,i,j)

!...        assume N = H * n at top of model
            o3t(nz,i,j) = tt * hnm(nz,i,j)
            o2t(nz,i,j) = o2prof(nz) * hcm
	    
            do k = nz-1,1, -1

	       o2t(k,i,j) = o2t(k+1,i,j) + 
     $	              dzcm * ( o2prof(k+1) - o2prof(k)) /
     $                   (alog( o2prof(k+1)/o2prof(k) ) ) 

!               o3t(k,i,j) = o3t(k+1,i,j) +
!     $	              dzcm * ( o3prof(k+1) - o3prof(k)) /
!     $                   (alog( o3prof(k+1)/o3prof(k) ) ) 

               o3t(k,i,j) = o3t(k+1,i,j) +
     $	              0.5 * dzcm * ( o3prof(k+1) + o3prof(k))

            end do
         end do
      end do

      end

!-----------------------------------------------------------
