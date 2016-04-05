      subroutine shapir_w

c-----------------------------------------------------------------------
c     Shapiro filter in the zonal and meridional direction for w 
c-----------------------------------------------------------------------

      use params
      use dynam

      implicit none

      integer i, j, k, n, iopp

      real :: div
      real, dimension(nx,ny) :: axw, bxw
      real, dimension(ny,nx) :: ayw, byw

      do k = 1,nz
	 div = 1./tey(k)
c-----------------------------------------------------------------------
c  zonal
         do j=1,ny
            do i=1,nx
               bxw(i,j)=wn1(k,i,j)
            end do
         
            do n=1,nifx(k)
               do i=2,nx-1
                  axw(i,j)=bxw(i-1,j)-2.*bxw(i,j)+bxw(i+1,j)
               end do
               axw(1,j)=bxw(nx,j)-2.*bxw(1,j)+bxw(2,j)
               axw(nx,j)=bxw(nx-1,j)-2.*bxw(nx,j)+bxw(1,j)

               do i=2,nx-1
                  bxw(i,j)=axw(i-1,j)-2.*axw(i,j)+axw(i+1,j)
               end do
               bxw(1,j)=axw(nx,j)-2.*axw(1,j)+axw(2,j)
               bxw(nx,j)=axw(nx-1,j)-2.*axw(nx,j)+axw(1,j)
            end do
         
            do i = 1,nx
               byw(j,i) = wn1(k,i,j) - bxw(i,j)*div
            end do
         end do
c-----------------------------------------------------------------------
c  meridional 
         do n = 1,nify(k)
            do i=1,nx
               do j = 2,ny-1
                  ayw(j,i) = byw(j-1,i) - 2.*byw(j,i) + byw(j+1,i)
               end do
               if (i.gt.nxhlf) iopp = i - nxhlf
               if (i.le.nxhlf) iopp = i + nxhlf
               ayw(1,i) =  byw(1,iopp) - 2.*byw(1,i) + byw(2,i)
               ayw(ny,i) = byw(ny-1,i) - 2.*byw(ny,i) + byw(ny,iopp)
            end do

            do i=1,nx
               do j = 2,ny-1
                  byw(j,i) = ayw(j-1,i) - 2.*ayw(j,i) + ayw(j+1,i)
               end do
               if (i.gt.nxhlf) iopp = i - nxhlf
               if (i.le.nxhlf) iopp = i + nxhlf
               byw(1,i) =  ayw(1,iopp) - 2.*ayw(1,i) + ayw(2,i)
               byw(ny,i) = ayw(ny-1,i) - 2.*ayw(ny,i) + ayw(ny,iopp)
            end do
         end do

         do i=1,nx
            do j = 1,ny
               wn1(k,i,j) = wn1(k,i,j) - div*byw(j,i)
	    end do
         end do
      end do

      end




