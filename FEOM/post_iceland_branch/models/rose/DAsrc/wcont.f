       subroutine wcont
c  w from continuity
c    note vertical grid is offset by 1/2

      use params
      use dynam

      implicit none

      integer :: i, j, k, jp1, jm1, ip1, im1
      real, dimension(nz) :: diverg

      do j=1,ny
         jp1=j+1
         if(j.eq.ny)jp1=ny
         jm1=j-1
         if(j.eq.1)jm1=1
         do i=1,nx
            ip1=mod(i,nx)+1
            im1=i-1
            if(i.eq.1) im1=nx
            do k=1,nz
               diverg(k)=((un1(k,ip1,j)-un1(k,im1,j))/dx(j)
     +                   + ((vn1(k,i,jp1)+vn1(k,i,j))*cosf2(j+1)
     +                   - (vn1(k,i,j)+vn1(k,i,jm1))*cosf2(j))/dy
     +                   /cosfi(j))*(dz*.5)*rou(k)
            end do
            wn1(nz,i,j)=diverg(nz)
            do k=nz-1,1,-1
               wn1(k,i,j) = wn1(k+1,i,j) + diverg(k)
            end do
         end do
      end do

c  filter in horizontal
      call shapir_w

c   interpolate w to model grid and convert units for transport
      do j=1,ny
         do i=1,nx
            do k = 1,nz-1
               ww(k,i,j) = (wn1(k+1,i,j)/row(k+1) + wn1(k,i,j)
     $                         /row(k))*.5
            end do
            ww(nz,i,j) = wn1(nz,i,j)/row(nz)*.5
         end do
      end do
      return
      end
