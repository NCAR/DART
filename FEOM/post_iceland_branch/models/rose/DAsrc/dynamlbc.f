
      subroutine dynamlbc( gmt_frac )

!  interpolate dynamical tropospheric and lower boundary data in time
!  13-day averages from Oct 1, 1992 to Dec 27, 1994
!  yearno is 0 for 1992, etc.

      use params, only : nx, ny, nz, ntime
      use dynam

      use DArose_mod, only : h_tune, z_tune

      implicit none

      real, intent(in) :: gmt_frac

      integer :: icall, i, j, nday4, nday24, lengdy, nt0
      real :: a, zlbc, zubc, thlbc, thubc, amplify, zt, zzt
      real :: xdiv, gmt_rad, drad, dzfctr1
      real :: rad_lt, cct1, sst1, cct2, sst2, zfactr, hfactr, rad_adj
      real, dimension(nx) :: cos1, sin1, cos2, sin2, xlamb

      save icall,nt0,cos1,sin1,cos2,sin2,xdiv,thlbc,thubc
      data icall/1/,a/6366197./

      print *, 'dynamlbc - gmt_frac: ', gmt_frac

      if(icall.eq.1) then
         icall=icall+1
         xdiv=1./float(nx)
         lengdy = ntime * 24
         nt0 = 3. * lengdy
         zlbc = z(1) - dz
         zubc = z(nz) + dz
         thlbc = exp(.28*zlbc/h)
         thubc = exp(.28*zubc/h)
         zfactr = exp(-(16.218-15.)/7.) * z_tune  ! tuning amplitude by an arbitrary factor "z_tune"
!         zfactr = exp(-(16.218-15.)/7.)       ! amplitude and phase adjustments
         hfactr = (16.218-15.)/25. * 2. * pi  ! ..  for GSWM forcing altitude
!         h_tune = pi    ! tuning to get observed phase at 95 km
      endif

      print *, 'dynamlbcV2: h_tune = ', h_tune

      ubc = u_lbc
      vbc = v_lbc
      tbc = t_lbc
      fibc = fi_lbc

      tubc = tnzp1 * thubc 
      uubc = unzp1

c  potential temperature
      do j=1,ny
         do i=1,nx
            thbc(i,j) = (tbc(i,j) + treflb) * thlbc
         end do
      end do

c  with tides
      gmt_rad = 2. * pi * gmt_frac
      drad = 2.*pi / float(nx)

      do i=1,nx
         rad_lt = gmt_rad + float(i-1) * drad
         rad_adj = rad_lt + hfactr + h_tune
         cct1 = cos(rad_adj)
         sst1 = sin(rad_adj)
         cct2 = cos(2.*rad_lt)
         sst2 = sin(2.*rad_lt)
         do j=1,36
            ubc(i,j)  = ubc(i,j) + (dau(j)*sst1 + dbu(j)*cct1) * zfactr
     $                           + (sau(j)*sst2 + sbu(j)*cct2)
            vbc(i,j)  = vbc(i,j) + (dav(j)*sst1 + dbv(j)*cct1) * zfactr
     $                           + (sav(j)*sst2 + sbv(j)*cct2)
            tbc(i,j)  = tbc(i,j) + (dat(j)*sst1 + dbt(j)*cct1) * zfactr
     $                           + (sat(j)*sst2 + sbt(j)*cct2)
            fibc(i,j) = fibc(i,j)+ (daz(j)*sst1 + dbz(j)*cct1) * zfactr
     $                           + (saz(j)*sst2 + sbz(j)*cct2)
c  potential temperature
            thbc(i,j) = (tbc(i,j) + treflb) * thlbc
         end do
      end do

      return
      end

