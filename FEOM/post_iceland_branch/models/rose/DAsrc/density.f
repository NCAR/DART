      subroutine density (tfull)
!------------------------------------------------------------------------
!  calculate atmospheric density, needed for chemistry and radiation
!------------------------------------------------------------------------
!     calculation of number/cm**3 from ideal gas law
!     hnm = p/kT = (pmb*1.e-4) / (T * 1.38e-23)
!--------------------------------------------------------------

      use params
      use dynam
      use chem, only : hnm

      implicit none

      real, dimension(nz,nx,ny) :: tfull
      integer i, j, k

!-------------------------------------------------------------------
!     number/cm^3 needed for chemistry routines
!-------------------------------------------------------------------
      do j=1,ny
         do i = 1,nx
            do k = 1,nz
               hnm(k,i,j) = pmb(k) /(tfull(k,i,j)*1.38e-19)
            end do
         end do
      end do

      end
