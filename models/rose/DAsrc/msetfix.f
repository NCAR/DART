c---------------------------------------------------------------------
      subroutine msetfix
c---------------------------------------------------------------------
c   initialize geometric variables
c---------------------------------------------------------------------

      use params
      use chem
      use dynam
      use phys

      implicit none

      integer :: i, j, k, n, nif(nz)
      real :: a, arg, zbot, xdiv, phirad(ny), co2prof(nzz)
 
      data nif/3*4,21*8,4*4,10*2/
      data a /6366197./
      data zbot/17500./
      data co2prof /35*355., 353., 350., 335., 300., 260., 200., 145., 
     $              105., 80., 50./

      dy = 555555.55
      defi = 5.
      pi = 3.141592654
      h = 7000.
      dz = 2500.
      dzrh2 = 287.04*dz*.5/h
      xdiv = 1./float(nx)

      do j=1,ny
         phideg(j) = (-85.+(j-1.5)*defi)
         phirad(j) = phideg(j)/180.*pi
         sinfi(j) = sin(phirad(j))
         cosfi(j) = cos(phirad(j))
         tgfia(j) = tan(phirad(j))/a
         cor(j) = 2.*7.292116e-5*sinfi(j)
         dx(j) = 2.*pi*xdiv*a*cosfi(j)
      end do
      do j=2,ny
         cosf2(j) = cos((-90.+(j-1)*defi)/180.*pi)
      end do
      cosf2(1) = 0.
      cosf2(nyp1) = 0.
 
      do k=1,nz
         z(k) = (k-1)*dz+zbot
         zkm(k) = z(k)*1.e-3
         pmb(k) = 1013.*exp(-z(k)/7.e3)
         rou(k) = exp(-z(k)/h)
         row(k) = exp(-(z(k)-0.5*dz)/h)
         thconv(k) = (1013./pmb(k))**.28
      end do

c  initialize array for fft
      do i=1,nxhlf
         arg = (i-1)*pi/float(nxhlf)
         wfft(i) = cmplx(cos(arg),sin(arg))
      end do

c   initialize for vertical diffusion
      gravit = 9.8
      rair = 287.04
      do k=1,nz
c         zkmin(k) = .5
         zkmin(k) = .01
         xml2(k) = 30.0**2
      end do

c   initialize for slt transport
      dxrad = 2.*pi/float(nx)
      dyrad = defi*pi/180.
      dzz = dz
      do j=1,ny
         gcx(j) = 1./dxrad/a/cosfi(j)
      end do
      gcy = 1./dyrad/a
      gcz = 1./dz
      ymin = phirad(1)

c----------------------------------------------------------------------
c  input for ir cooling
c----------------------------------------------------------------------
      co2mr = co2prof * 1.e-6
      tcltop = 270.

c  filter parameters
      do k=1,nz
         tex(k)=float(4**nif(k))
         tey(k)=float(4**nif(k))
         nifx(k)=nif(k)/2
         nify(k)=nif(k)/2
      end do

      end










