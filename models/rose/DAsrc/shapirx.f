      subroutine shapirx

      use params
      use dynam

      implicit none

      integer i, j, k, n
      real :: div
      real, dimension(nx,ny) :: au2, av2, at2, bu2, bv2, bt2

      do k=1,nz
         div=1./tex(k)
         do j=1,ny
            do i=1,nx
               bu2(i,j)=un2(k,i,j)
               bv2(i,j)=vn2(k,i,j)
               bt2(i,j)=tn2(k,i,j)
            end do
         
            do n=1,nifx(k)
               do i=2,nx-1
                  au2(i,j)=bu2(i-1,j)-2.*bu2(i,j)+bu2(i+1,j)
                  av2(i,j)=bv2(i-1,j)-2.*bv2(i,j)+bv2(i+1,j)
                  at2(i,j)=bt2(i-1,j)-2.*bt2(i,j)+bt2(i+1,j)
               end do
               au2(1,j)=bu2(nx,j)-2.*bu2(1,j)+bu2(2,j)
               av2(1,j)=bv2(nx,j)-2.*bv2(1,j)+bv2(2,j)
               at2(1,j)=bt2(nx,j)-2.*bt2(1,j)+bt2(2,j)
               au2(nx,j)=bu2(nx-1,j)-2.*bu2(nx,j)+bu2(1,j)
               av2(nx,j)=bv2(nx-1,j)-2.*bv2(nx,j)+bv2(1,j)
               at2(nx,j)=bt2(nx-1,j)-2.*bt2(nx,j)+bt2(1,j)

               do i=2,nx-1
                  bu2(i,j)=au2(i-1,j)-2.*au2(i,j)+au2(i+1,j)
                  bv2(i,j)=av2(i-1,j)-2.*av2(i,j)+av2(i+1,j)
                  bt2(i,j)=at2(i-1,j)-2.*at2(i,j)+at2(i+1,j)
               end do
               bu2(1,j)=au2(nx,j)-2.*au2(1,j)+au2(2,j)
               bv2(1,j)=av2(nx,j)-2.*av2(1,j)+av2(2,j)
               bt2(1,j)=at2(nx,j)-2.*at2(1,j)+at2(2,j)
               bu2(nx,j)=au2(nx-1,j)-2.*au2(nx,j)+au2(1,j)
               bv2(nx,j)=av2(nx-1,j)-2.*av2(nx,j)+av2(1,j)
               bt2(nx,j)=at2(nx-1,j)-2.*at2(nx,j)+at2(1,j)
            end do
         
            do i=1,nx
               un2(k,i,j)=un2(k,i,j)-bu2(i,j)*div
               vn2(k,i,j)=vn2(k,i,j)-bv2(i,j)*div
               tn2(k,i,j)=tn2(k,i,j)-bt2(i,j)*div
            end do
         end do
      end do
      return
      end
