      subroutine shapiry

c-----------------------------------------------------------------------
c     Shapiro filter in the meridional direction - with polar wraparound
c-----------------------------------------------------------------------

      use params
      use dynam

      implicit none

      integer i, j, k, n, iopp

      real :: div
      real, dimension(ny,nx) :: au2, av2, at2, bu2, bv2, bt2

      do k = 1,nz
	 div = 1./tey(k)
         do i = 1,nx
            do j = 1,ny
               bu2(j,i) = un2(k,i,j)
               bv2(j,i) = vn2(k,i,j)
               bt2(j,i) = tn2(k,i,j)
            end do
         end do

         do n = 1,nify(k)
            do i=1,nx
               do j = 2,ny-1
                  au2(j,i) = bu2(j-1,i) - 2.*bu2(j,i) + bu2(j+1,i)
                  av2(j,i) = bv2(j-1,i) - 2.*bv2(j,i) + bv2(j+1,i)
                  at2(j,i) = bt2(j-1,i) - 2.*bt2(j,i) + bt2(j+1,i)
               end do
               if (i.gt.nxhlf) iopp = i - nxhlf
               if (i.le.nxhlf) iopp = i + nxhlf
               au2(1,i) = -bu2(1,iopp) - 2.*bu2(1,i) + bu2(2,i)
               av2(1,i) = -bv2(1,iopp) - 2.*bv2(1,i) + bv2(2,i)
               at2(1,i) =  bt2(1,iopp) - 2.*bt2(1,i) + bt2(2,i)
               au2(ny,i) = bu2(ny-1,i) - 2.*bu2(ny,i) - bu2(ny,iopp)
               av2(ny,i) = bv2(ny-1,i) - 2.*bv2(ny,i) - bv2(ny,iopp)
               at2(ny,i) = bt2(ny-1,i) - 2.*bt2(ny,i) + bt2(ny,iopp)
            end do

            do i=1,nx
               do j = 2,ny-1
                  bu2(j,i) = au2(j-1,i) - 2.*au2(j,i) + au2(j+1,i)
                  bv2(j,i) = av2(j-1,i) - 2.*av2(j,i) + av2(j+1,i)
                  bt2(j,i) = at2(j-1,i) - 2.*at2(j,i) + at2(j+1,i)
               end do
               if (i.gt.nxhlf) iopp = i - nxhlf
               if (i.le.nxhlf) iopp = i + nxhlf
               bu2(1,i) = -au2(1,iopp) - 2.*au2(1,i) + au2(2,i)
               bv2(1,i) = -av2(1,iopp) - 2.*av2(1,i) + av2(2,i)
               bt2(1,i) =  at2(1,iopp) - 2.*at2(1,i) + at2(2,i)
               bu2(ny,i) = au2(ny-1,i) - 2.*au2(ny,i) - au2(ny,iopp)
               bv2(ny,i) = av2(ny-1,i) - 2.*av2(ny,i) - av2(ny,iopp)
               bt2(ny,i) = at2(ny-1,i) - 2.*at2(ny,i) + at2(ny,iopp)
            end do
         end do

         do i=1,nx
            do j = 1,ny
               un2(k,i,j) = un2(k,i,j) - div*bu2(j,i)
               vn2(k,i,j) = vn2(k,i,j) - div*bv2(j,i)
               tn2(k,i,j) = tn2(k,i,j) - div*bt2(j,i)
	    end do
         end do
      end do

      end




