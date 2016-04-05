      subroutine filter

      use params
      use dynam

      implicit none

      integer, parameter :: nx2=nx*2

      complex x(nx2)
      integer :: i, j, k, nv, nvar, inv, jlow, jhigh, nfilt(ny)
      integer :: ifi, ifiend

      data nfilt/2,34*0,2/, jlow/1/, jhigh/ny/
      nvar=3

      do j=1,ny
         if(j.le.jlow.or.j.ge.jhigh)then
            do k=1,nz
               do nv = 1,nvar
                  if(nv.eq.1)then
                     do i=1,nx
                        x(i) = un2(k,i,j)
                     end do
                  endif
                  if(nv.eq.2)then
                     do i=1,nx
                        x(i) = vn2(k,i,j)
                     end do
                  endif
                  if(nv.eq.3)then
                     do i=1,nx
                        x(i) = tn2(k,i,j)
                     end do
                  endif

                  inv = -1
                  call fft(x,inv,wfft)
                  ifi = nfilt(j) + 2
                  ifiend = nx - ifi + 2
                  do i=ifi,ifiend
                     x(i) = 0.
                  end do
                  inv = 1
                  call fft(x,inv,wfft)
                  if(nv.eq.1)then
                     do i=1,nx
                        un2(k,i,j) = real(x(i))
                     end do
                  endif
                  if(nv.eq.2)then
                     do i=1,nx
                        vn2(k,i,j) = real(x(i))
                     end do
                  endif
                  if(nv.eq.3)then
                     do i=1,nx
                        tn2(k,i,j) = real(x(i))
                     end do
                  endif
               end do
            end do
         endif
      end do
      return
      end
