      subroutine fft (x,inv,w)
c   fft
c   (note that variable iter = log(base 2) of nx

      use params
      implicit none

      integer, parameter :: nx2=nx*2

      integer :: iter, i, j, j1, j2, k, nxp2, it, n, n1, n2, nxp, m, mxp
      integer :: inv, itab(6)
      complex :: x(nx2), w(nxhlf), t, wk

      data iter/5/
      data itab/1, 2, 4, 8, 16, 32/

      nxp2 = nx
      do it=1,iter
         n = 1
         nxp = nxp2
         nxp2 = nxp/2

         do m=1,nxp2
            wk = w(n)
            if(inv.gt.0)wk = conjg(wk)
            do mxp = nxp,nx,nxp
               j1 = mxp - nxp + m
               j2 = j1 + nxp2
               t = x(j1) - x(j2)
               x(j1) = x(j1) + x(j2)
               x(j2) = t*wk
            end do
            n = n +itab(it)
         end do
      end do

      n2 = nx/2
      n1 = nx - 1
      j = 1
 
      do i=1,n1
         if (i.lt.j) then
            t = x(j)
            x(j) = x(i)
            x(i) = t
         end if
         k = n2
 42      if(k.lt.j) then
            j = j - k
            k = k/2
            go to 42
         end if
         j = j + k
      end do

      if(inv.lt.0) return
      do i=1,nx
         x(i) = x(i)/float(nx)
      end do

      return
      end
