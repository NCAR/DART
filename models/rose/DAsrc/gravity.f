      subroutine gravity
C-----------------------------------------------------------------------
c  Lindzen's gravity wave drag for 3-d zonal wind
C-----------------------------------------------------------------------
      use params
      use dynam
      use phys

      implicit none
C-----------------------------------------------------------------------
      integer, parameter :: nfreq=7

      integer i, j, k, mm, nn, ic, ib, nsmfx, nsmdz, nsmfa

      real, dimension(nz) :: umc, uz, umclimit, alpha, zalpha
      real, dimension(nz) :: bv, rnk, hnk, bv2, ff
      real, dimension(nfreq) ::  cgrav, amp
      real utilda, gcp, rh, umc0, phi0, afact, xc, sw, zb, xb, zint
      real amplat(ny), ampnew, efact

      data cgrav /-15.,-10.,-5., 0., 5.,10.,15./
      data amp /1.e-11, 1.e-11, 1.e-11, 1.e-11, 1.e-11, 1.e-11, 1.e-11/

      data utilda / 5. /
      data gcp /9.8e-03/, rh/.04/

c   initialize gravity wave terms

      do k=1,nz
         alpha(k) = 3.e-6 * (1.2 + tanh((zkm(k)-50.)/12.))
      end do
      umc0 = 20.
      phi0 = 2.5

      afact = 8.

      do j = 1,ny
         amplat(j) = (sin((j-18.5)*10./180.*pi)) ** 2
         if(j.le.9)amplat(j) = 1.
         if(j.ge.28)amplat(j) = 1.
         if(j.le.18)amplat(j) = amplat(j) * .7  ! "orographic" gw weaker in SH
         do i=1,nx
            do k=1,nz
               fdzz(k,i,j) = 0.
               fgr(k,i,j) = 0.
               fcgr(k,i,j) = 0.
               falph(k,i,j) = 0.
            end do
            do k=2,nz-1
               uz(k) = (un1(k+1,i,j)-un1(k-1,i,j))/(2.*dz)
               bv2(k) = rh * ((tn1(k+1,i,j)-tn1(k-1,i,j))/(2.*dz)
     $                       + ssref(k))
               bv(k) = sqrt(amax1(bv2(k), 1.e-6))
               hnk(k) = 7.e3 * bv(k) / 8.e-5
               rnk(k) = bv(k) / 8.e-5
            end do
            uz(1) = (un1(2,i,j)-un1(1,i,j))/dz
            bv2(1) = bv2(2)
            hnk(1) = hnk(2)
            rnk(1) = rnk(2)
            uz(nz) = (un1(nz,i,j)-un1(nz-1,i,j))/dz
            bv2(nz) = bv2(nz-1)
            hnk(nz) = hnk(nz-1)
            rnk(nz) = rnk(nz-1)

            do nn = 1,nfreq
               ampnew = (.6 + 1. * amplat(j)) * amp(nn)
               do k = 1,nz
                  umc(k) =  un1(k,i,j) - cgrav(nn)
                  if(abs(umc(k)).gt.1.) umclimit(k) = umc(k)
                  if(abs(umc(k)).lt.1.) umclimit(k) = sign(1.,umc(k))
               end do

c     calculate critical level height
               ic = nz
               xc =  umc(1)
               do k = nz,1,-1
                  sw = xc * umc(k)
                  if (sw.lt.0.) ic = k
               end do

c     adjustment to breaking level from damping
               zalpha(1) = hnk(1)*alpha(1) / umclimit(1)**2*dz
               do k = 2,ic
                  zalpha(k) = zalpha(k-1) + hnk(k)*alpha(k)
     $                                    / umclimit(k)**2*dz
               end do
c     calculate breaking level - adjusted for damping
               ib = nz
               do k = ic,1,-1
                  zb = 21.e3 * alog (abs(umclimit(k))/utilda )
     $                  + zalpha(k)
                  if (k.lt.nz) xb = z(k+1) - zb
                  if (k.eq.nz) xb = z(k) + dz - zb
                  if (xb.ge.0.) ib = k
               end do

c     calculate momentum drag and eddy diffusion coefficient.
c     note that modifications to fdzz for Pr>1 are applied in subr. vertdiff
               if( ib .le. ic ) then
                  do k = ib,ic
                     ff(k) = - ampnew * afact * umc(k) * 
     $                            (umc(k)-21.e3*uz(k))
                     fgr(k,i,j) = fgr(k,i,j) + ff(k)
                     fcgr(k,i,j) = fcgr(k,i,j) 
     $                              + ff(k) * cgrav(nn)
                     fdzz(k,i,j) = fdzz(k,i,j) + abs(ff(k)*umc(k)**2
     $                                        /bv2(k))
                  end do
c     exponential decay below zb
c     also include damping term from Holton and Zhu to account for momentum
c          forcing at critical levels
                  zint = 0.
                  do k = 1,ib
                     ff(k) = ff(ib) * exp((z(k)-z(ib))/h)
                     fgr(k,i,j) = fgr(k,i,j) + ff(k)
                     fcgr(k,i,j) = fcgr(k,i,j)+ff(k)*cgrav(nn)
                     fdzz(k,i,j) = fdzz(k,i,j) + abs(ff(k)*umc(k)**2
     $                                 /bv2(k))

                     zint = zint + rnk(k) * alpha(k)/
     $                                         (umclimit(k)**2)*dz
                     efact = exp(z(k)/h - zint)
                     falph(k,i,j) = falph(k,i,j) 
     $                         - phi0 * alpha(k) * efact /
     $                     (2. * umc0 * umclimit(k) * abs(umclimit(k)))
                  end do
               end if
            end do

         end do
      end do

!c-----------------------------------------------------------------------
!c  smooth in the zonal direction
!      nsmfx = 3
!      nsmdz = 3
!      nsmfa = 3
!      do nn=1,3
!         call smoothx( fgr, 1, nz, nsmfx )
!         call smoothx( fcgr, 1, nz, nsmfx )
!         call smoothx( falph, 1, nz, nsmfa )
!c  smooth in the meridional direction
!         call smoothl( fgr, 1, nz, nsmfx )
!         call smoothl( fcgr, 1, nz, nsmfx )
!         call smoothl( falph, 1, nz, nsmfa )
!c  smooth in the vertical direction
!         call smoothv ( fgr, 1, nz, nsmfx )
!         call smoothv ( fcgr, 1, nz, nsmfx )
!         call smoothv ( falph, 1, nz, nsmfa )
!      end do
!      call smoothx( fdzz, 1, nz, nsmdz )
!      call smoothl( fdzz, 1, nz, nsmdz )
!      call smoothv ( fdzz, 1, nz, nsmdz )

      return
      end

      subroutine smoothx( y, k1, k2, nsmol)
      use params
      implicit none

      real, dimension(nx) :: x
      real, dimension(nz,nx,ny) :: y
      integer i, j, k, k1, k2, n, nsmol

      do j = 1,ny
         do n = 1,nsmol
            do k = k1,k2
               x(1) = (y(k,nx,j)+2.*y(k,1,j)+y(k,2,j))/4.
               do i = 2,nx-1
                  x(i) = (y(k,i-1,j)+2.*y(k,i,j)+y(k,i+1,j))/4.
               end do
               x(nx) = (y(k,nx-1,j)+2.*y(k,nx,j)+y(k,1,j))/4.
               do i = 1,nx
                 y(k,i,j) = x(i)
               end do
            end do
         end do
      end do
      end

      subroutine smoothl( y, k1, k2, nsmol)
      use params
      implicit none

      real, dimension(ny) :: x
      real, dimension(nz,nx,ny) :: y
      integer i, j, k, k1, k2, n, nsmol

      do i = 1,nx
         do n = 1,nsmol
            do k = k1,k2
               x(1) = (2.*y(k,i,1)+y(k,i,2))/3.
               do j = 2,ny-1
                  x(j) = (y(k,i,j-1)+2.*y(k,i,j)+y(k,i,j+1))/4.
               end do
               x(ny) = (y(k,i,ny-1)+2.*y(k,i,ny))/3.
               do j = 1,ny
                  y(k,i,j) = x(j)
               end do
            end do
         end do
      end do
      end

      subroutine smoothv( y, k1, k2, nsmov)
      use params
      implicit none

      real, dimension(nz) :: x
      real, dimension(nz,nx,ny) :: y
      integer i, j, k, k1, k2, n, nsmov

      do i = 1,nx
         do n = 1,nsmov
            do j = 1,ny
               do k = k1+1,k2-1
                  x(k) = (y(k-1,i,j)+2.*y(k,i,j)+y(k+1,i,j))/4.
               end do
               x(k2) = (y(k2,i,j)*2.+y(k2-1,i,j))/3.
               x(k1) = (y(k1,i,j)*2.+y(k1+1,i,j))/3.
               do k = k1,k2
                  y(k,i,j) = x(k)
               end do
            end do
         end do
      end do
      end


