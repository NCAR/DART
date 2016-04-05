      subroutine hines_driver
c==================================================================
c    driver for hines gravity wave drag

      use params
      use dynam
      use phys

      implicit none

      integer, parameter :: nlons=nx
      integer, parameter :: nlevs=nz+2
      integer, parameter :: nazmth = 8

      integer  :: icall, i, j, k, nn, nsmfx, nsmdz, nsmfa
      integer  :: naz, il1, il2, levbot, levtop, nmessg, ierror
      integer  :: icutoff, nsmax, iheatcal, iprnt1, iprnt2
      logical  :: dragil(nlons,nlevs), drag(nlons)
      logical  :: do_alpha(nlons,nazmth)
      real  :: kstar, m_min, f1, f2, f3, f5, f6, slope, alt_cutoff, 
     $         smco, vv, rh, denref, sch, zz, bv2
      real  :: drag_u(nlons,nlevs),   drag_v(nlons,nlevs) 
      real  :: gwheat(nlons,nlevs),   gwdiff(nlons,nlevs)
      real  :: flux_u(nlons,nlevs),   flux_v(nlons,nlevs)
      real  :: vel_u(nlons,nlevs),    vel_v(nlons,nlevs)
      real  :: bvfreq(nlons,nlevs),   density(nlons,nlevs)
      real  :: visc_mol(nlons,nlevs), alt(nlons,nlevs)
      real  :: rms_wind(nlons),       bvfb(nlons),   densb(nlons)
      real  :: ubot(nlons),           vbot(nlons)
      real  :: sigma_t(nlons,nlevs)
      real  :: sigma_alpha(nlons,nlevs,nazmth)
      real  :: sigsqh_alpha(nlons,nlevs,nazmth)
      real  :: m_alpha(nlons,nlevs,nazmth), v_alpha(nlons,nazmth)
      real  :: ak_alpha(nlons,nazmth),      k_alpha(nlons,nazmth)
      real  :: mmin_alpha(nlons,nazmth) ,   i_alpha(nlons,nazmth)
      real  :: specfac(nlons,nazmth),  kstar2(ny)
      real :: dtrop(nztrop),tprof(nzz),sprof(nzz),ptrop(nztrop)
      real :: denzz(nzz)
      save iheatcal, icutoff, iprnt1, iprnt2, nsmax, smco, 
     $     alt_cutoff, kstar, m_min, slope, f1, f2, f3, f5, f6, naz,
     $     il1, il2, levbot, levtop, visc_mol, kstar2

      data icall/0/, rh/.04/

c********************************************************************
c  setup for Hines gravity wave drag
      if(icall.eq.0) then
         nmessg = 6
         call hines_setup (ierror, naz, slope, f1, f2, f3, f5, f6, 
     $                  kstar, m_min, icutoff, alt_cutoff, smco, nsmax, 
     $                  iheatcal, nmessg, nazmth)
         if (ierror.ne.0) then
            print 10,ierror
 10         format('error in setup of Hines gwd')
            stop
         end if
         icall = icall + 1

c  coefficient for molecular viscosity - from SOCRATES
         do k=1,nz
            denref = pmb(k)/(tref(k)*1.38e-19)
            vv =  1.52e14 * sqrt(2./28.9*tref(k))/denref
            do i=1,nlons
               visc_mol(i,k+2) = vv
            end do
         end do
         do k=6,nztrop
            ptrop(k) = 1013. * exp (-(k-1)*2.5/7.)
            denref = ptrop(k)/(troptref(k)*1.38e-19)
            vv =  1.52e14 * sqrt(2./28.9*troptref(k)) / denref
            do i=1,nlons
               visc_mol(i,k-5) = vv
            end do
         end do
c  mean density (kg/m^3)
         do k=1,nz
            sch = 287. * tref(k) / 9.8
            denzz(k+2) = 1.293 * exp(-z(k)/7.e3)
         end do
         do k=6,nztrop
            zz = (k-1)*dz
            sch = 287. * troptref(k) / 9.8
            denzz(k-5) = 1.293 * exp(-zz/7.e3)
         end do
c  other inputs
         il1 = 1
         il2 = nlons
         levbot = 1
         levtop = nz+2
         do j=1,ny
c            kstar2(j) = 1.5e-6 / (.7*cosfi(j)+.4)
            kstar2(j) = 0.8e-6 * (.7*cosfi(j)+.4)
         end do
      end if

c********************************************************************

      do j=1,ny
         kstar = kstar2(j)
         do i=1,nlons
c            rms_wind(i) = 1.5
            rms_wind(i) = 1.8
            do k=1,nz
               vel_u(i,k+2) = un1(k,i,j)
               vel_v(i,k+2) = vn1(k,i,j)
               density(i,k+2) = denzz(k+2)
               alt(i,k+2) = z(k)
               tprof(k+2) = tn1(k,i,j)
               sprof(k+2) = ssref(k)
            end do
            do k=6,nztrop
               vel_u  (i,k-5) = unmc(k,j)
               vel_v  (i,k-5) = 0.
               density(i,k-5) = denzz(k)
               alt    (i,k-5) = (k-1)*dz
               tprof  (k-5) = tpnmc(k,j)
               sprof  (k-5) = tropsref(k)
            end do
            do k=2,nz+1
               bv2 = rh * ((tprof(k+1)-tprof(k-1))/(2.*dz)
     $                       + sprof(k))
               bvfreq(i,k) = sqrt(amax1(bv2, 1.e-6))
            end do
            bvfreq(i,1) = bvfreq(i,2)
            bvfreq(i,nz+2) = bvfreq(i,nz+1)
            do nn=1,nazmth
               k_alpha(i,nn) = kstar
            end do
         end do

         call hines_dsp1 (drag_u, drag_v, gwheat, gwdiff,         ! output
     $                    flux_u, flux_v,                         ! .
     $                    vel_u, vel_v, bvfreq, density,          ! input
     $                    visc_mol, alt, rms_wind, k_alpha,       ! .
     $                    v_alpha, m_alpha, sigma_alpha,          ! work arrays
     $                    sigsqh_alpha, ak_alpha, specfac,        ! .
     $                    do_alpha, drag, dragil, mmin_alpha,     ! .
     $                    i_alpha, sigma_t, densb, bvfb, ubot,    ! .
     $                    vbot,                                   ! .
     $                    iheatcal, icutoff,                      ! hines_setup
     $                    iprnt1, iprnt2,                         ! input
     $                    nsmax, smco, alt_cutoff, kstar, m_min,  ! hines_setup
     $                    slope, f1, f2, f3, f5, f6, naz,         ! .
     $                    il1, il2, levbot, levtop, nlons, nlevs, ! input
     $                    nazmth)                                 ! .

         do k=1,nz
            do i=1,nlons
               fxh(k,i,j) = drag_u(i,k+2)
               fyh(k,i,j) = drag_v(i,k+2)
               fdzzh(k,i,j) = amax1(gwdiff(i,k+2),0.1)
            end do
         end do
      end do

      nsmfx = 2
      nsmdz = 2
      nsmfa = 2
c  smooth in the zonal direction
c      call smoothx (fxh, 1, nz, nsmfx)
c      call smoothx (fyh, 1, nz, nsmfx)
c      call smoothx (fdzzh, 1, nz, nsmdz)
c  smooth in the meridional direction
c      call smoothl (fxh, 1, nz, nsmfx)
c      call smoothl (fyh, 1, nz, nsmfx)
c      call smoothl (fdzzh, 1, nz, nsmdz)
c  smooth in the vertical direction
c      call smoothv (fxh, 1, nz, nsmfx)
c      call smoothv (fyh, 1, nz, nsmfx) 
c      call smoothv (fdzzh, 1, nz, 4)

      return
      end
