      subroutine dmolec
c  molecular thermal diffusion from Banks and Kockarts, adapted from 2D
c  assuming constant mixing ratio for O2 and N2
      use params
      use dynam
      use phys

      implicit none

      integer :: lowestz, k
      real :: g, cp, r, po, wm, const(nz)

      data g/9.8/, cp/1005./, r/8314./, po/1.013e05/, wm/28.9/

      lowestz = nz-4  !  for altitudes 100 km or higher

      do k=lowestz,nz
         const(k) = (g*h)**2 * wm / (po * cp * exp(-z(k)/h) * r)
         dtzz(k) = const(k) * (56. * tref(k)**.69) / 
     $                       tref(k) * 1.e-5
      end do
      do k=1,lowestz-1
         dtzz(k) = 0.
      end do

      end
