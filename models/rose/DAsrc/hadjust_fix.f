!-------------------------------------------------------------------
      subroutine hadjust_fix
!-------------------------------------------------------------------
!   adjust cooling rate so that global mean heating is zero
!-------------------------------------------------------------------

      use params
      use dynam
      use phys

      implicit none

      integer :: i, j, k, k2
      real :: ysum, divx, sum, tglobal, eps(nz), qtemp(ny)

      ysum = 0.
      do j=1,ny
         ysum = ysum + cosfi(j)
      end do
      divx = 1./float(nx)

      do k = 1,nz
         k2 = k + nztrop
         eps(k) = 0.
         do j=1,ny
            sum = 0.
            do i=1,nx
               sum = sum + solht(k,i,j) + qir(k2,i,j) + ch_heat(k,i,j)
            end do
            eps(k) = eps(k) + sum * cosfi(j)
         end do
         eps(k) = eps(k) * divx / ysum
     $          + dtglob(k)
     $          + dtadv(k)
         do j=1,ny
            do i=1,nx

!... index changed from k2 to k
!              q_resid(k2,i,j) = - eps(k)
!...
               q_resid(k,i,j) = - eps(k)

            end do
         end do
      end do

!... adjust perturbation temperature so that global mean is zero

      do k = 1,nz
         tglobal = 0.
         do j=1,ny
            sum = 0.
            do i=1,nx
               sum = sum + tn1(k,i,j)
            end do
            tglobal = tglobal + sum * cosfi(j)
         end do
         tglobal = tglobal * divx / ysum
         do j=1,ny
            do i=1,nx
               tn1(k,i,j) = tn1(k,i,j) - tglobal
            end do
         end do
      end do

      end
