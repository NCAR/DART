      subroutine hydrost
c  integrate hydrostatic eq to get fi
      use params
      use dynam

      implicit none

      integer :: i, j, k

      do j=1,ny
         do i=1,nx
            fin1(1,i,j) = fibc(i,j) + (tn1(1,i,j)+tbc(i,j))*dzrh2
            do k=1,nz-1
               fin1(k+1,i,j) = fin1(k,i,j) 
     $                       + (tn1(k,i,j)+tn1(k+1,i,j))*dzrh2
            end do
         end do
      end do
      end
