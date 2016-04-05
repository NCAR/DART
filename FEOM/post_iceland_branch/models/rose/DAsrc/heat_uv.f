!----------------------------------------------------------------------
      subroutine heat_uv (temp)
!----------------------------------------------------------------------
!
!  solar heating by o3 and o2 (k sec^-1) based on code 'heatmd.for' by Xun Zhu
!  used to calculate solar heating rate for middle atmosphere modeling (8/95)
!
!----------------------------------------------------------------------
!
!  input
!
!     temp    = temperature array
!
!  output
!     
!     solht() in phys.mod.f is modified
!
!  changes made by A. K. Smith and D. Marsh
!
!     - single precision
!     - add model loops to driver program
!     - move zenith angle calculation to zenith.f
!     - remove tropospheric water heating
!     - convert heating rate to K/sec
!     - functions incorporated into main subroutines
!     - heating efficiencies (function of pressure only) calculated just once 
!
!----------------------------------------------------------------------

      use params
      use phys, only : solht
      use chem, only : qn1,     ! constituent mixing ration array
     $                 hnm,     ! total density
     $                 o3t,     ! ozone column (3D)
     $                 o2t,     ! o2 column (3D)
     $                 sc2d     ! sec SZA (2D)

      implicit none

      real, dimension(nz,nx,ny) :: temp 

      integer :: i, j

      real, dimension(nz) :: tz, o3z, o2z, bo3, bo2, huvir
      real :: coszen 

      real :: ome0              ! scattering albedo

      real, dimension(ny) :: omega = (/
     $      .882, .783, .696, .619, .552, .495, .446, .404, .370,
     $      .342, .321, .304, .292, .283, .278, .274, .273, .273,
     $      .273, .273, .274, .278, .283, .292, .304, .321, .342,
     $      .370, .404, .446, .495, .552, .619, .696, .783, .882 /)


      solht(:,:,:) = 0. ! zero heating rates

      do j=1,ny

      ome0 = omega(j)           ! set albedo

!... select column data

         do i=1,nx

	    coszen = 1./sc2d(i,j)

            if (coszen .gt. 0.02) then          !  daytime  

              tz(:) = temp(:,i,j)

              !... convert o3 and o2 to molecules/m^3

              o3z(:) = qn1(:,i,j,19) * hnm(:,i,j) * 1.e6 

              o2z(:) = qn1(:,i,j,26) * hnm(:,i,j) * 1.e6 

              bo3(:) = o3t(:,i,j) * 1.e4     ! convert to 1/m^2

	      bo2(:) = o2t(:,i,j) * 1.e4     ! convert to 1/m^2

              !... calculate column heating rate

              call heat3(huvir, coszen, o2z, o3z, bo2, bo3, tz, ome0)

              solht(:,i,j) = huvir(:)

	   endif

         end do
      end do

      return
      end

!----------------------------------------------------------------------
      subroutine heat3(qh, bmu1, o2, o3, bo2, bo3, tem, ome0)
!----------------------------------------------------------------------
!
!       to calculate the heating rate (qh: k/sec) at time (hour) for o3 
!       and o2, following the parameterization scheme by Strobel (1978).
!
!       bmu1    = cos of zenith angle
!       o2, o3  = number density (m^-3)
!       bo2,bo3 = column density (m^-2).
!       pre     = pressure (pa)
!       tem     = temperature
!
!----------------------------------------------------------------------

      use params, only : nz
      use dynam,  only : pmb
      use phys,   only : rair

      implicit none

      real :: bmu1, ome0
      real, dimension(nz) :: qh, o2, o3, bo2, bo3, tem 

      real, external :: htmu
      real, external :: effcy 

      integer :: i
      real :: taust, p1
      real :: rhoi

      real, parameter :: cp_dry = 1004. ! specific heat of dry air

!... flag for vertical array increases upward
      integer, parameter :: itop = 0

      real, save, dimension(nz) :: eff_o2, eff_o3

      logical, save :: firstcall
      data firstcall /.true./

!... calculate heating efficiency factors (done once/run)

      if (firstcall)then
	 firstcall = .false.

	 print *, 'init efficiencies'

         do i=1,nz

           p1=pmb(i)

	   eff_o2(i) = effcy(p1,2)
	   eff_o3(i) = effcy(p1,3)

	 end do

      endif

!... zero heating rate

      qh(:) = 0.0

      if (itop.eq.0) then
        taust = bo3(1)
      else if (itop.eq.1) then
        taust = bo3(nz)
      else
        print *, 'heat_uv: itop not defined'
	stop
      end if

      do i=1,nz

        qh(i) = htmu(bmu1, o2(i), o3(i), bo2(i), bo3(i),
     $               ome0, taust, eff_o2(i), eff_o3(i) )
	
        rhoi = (pmb(i) * 1.0e2) / (rair * tem(i))
	
        qh(i) = qh(i) / (cp_dry * rhoi)

      end do 

      return
      end
      
!----------------------------------------------------------------------
      function htmu(bmu, o2, o3, bo2, bo3, ome0, bo3t, eff_o2, eff_o3)
!----------------------------------------------------------------------
!
!     heating rate (mks units) for o3 and o2, following the
!     parameterization scheme by strobel (1978) over six bands: 
!
!          chappius            (7500-4500 A)
!          huggins             (3550-2825 A)
!          hartley             (2725-2450 A)
!          herzberg continuum  (2400-2060 A)
!          schumann-runge      (2025-1750 A)
!          schumann-runge cont (1740-1260 A)
!
!     inputs:
!         bmu = cos(theta)
!         o2, o3 = number density (m^-3)
!         bo2, bo3 = column density (m^-2). 
!
!----------------------------------------------------------------------

      implicit none

      real :: htmu
      real :: bmu, o2, o3, bo2, bo3, ome0, bo3t, eff_o2, eff_o3 

      real :: xrm, sno3, sno2

      real :: hart, hugg2, hugg3, hubb
      real :: chap, htsum1, att, htsum2

      real :: src1, src21, src22, src23,  htsrc, htscat
      real :: fact1, fact2, htsrb, tau, taus, fac1 

!... magnification factor (vertical to slant column)
!    note:  previous formulation was designed to work for cos(SZA) = 0
!           xrm = 35.0/sqrt(1224.0*bmu**2+1.0)
!           However, this is not consistent with photolysis calculation
!           and so the standard sec(SZA) method is now used since sec(SZA)
!           is always less than or equal to 50.

      xrm = 1.0/bmu

      sno3 = bo3*xrm
      sno2 = bo2*xrm
      
!----------------------------------------------------------------------
!... ozone absorption
!----------------------------------------------------------------------

!... hartley band

      hart = 5.13 * 8.70e-22 * exp(-amin1(2.0e2, 8.70e-22*sno3))
     $        * eff_o3                              ! >> note efficieny <<

!... huggins band: long and short

      hugg2 = -2.0e-2 * exp(-amin1(2.0e2, 2.470e-23*sno3))
     $         / (sno3*0.01273)                          !   long

      hugg3 = -5.0e-2 * exp(-amin1(2.0e2, 3.573e-22*sno3))
     $         / (sno3*0.01273)                          !   short

      hubb=5.49882/sno3+hugg2+hugg3

!... chappius band

      chap = 370.0 * 2.85e-25 * exp(-amin1(2.0e2, 2.85e-25*sno3))

!... total ozone 

      htsum1=o3*(chap+hart+hubb)

!... herzberg band

      att=exp( -amin1(2.0e2, 5.0e-28 * sno2 + 4.0e-22 * sno3) )
      htsum2 = 1.5 * (5.0e-28 * o2 + 4.0e-22 * o3) * att

!----------------------------------------------------------------------
!... oxygen absorption
!----------------------------------------------------------------------

!... schumann-runge continuum 

      src1 = 0.65e-3 * 1.1e-21 * exp(-amin1(2.0e2, 1.1e-21*sno2))

      src21 = 1.23e-3 * exp(-amin1(2.0e2, 3.0e-23*sno2))/sno2

      src22 = -8.3e-4 * exp(-amin1(2.0e2, 2.0e-22*sno2))/sno2

      src23 = -4.0e-4 * exp(-amin1(2.0e2, 1.5e-21*sno2))/sno2

      htsrc = o2 * (src1+src21+src22+src23) * eff_o2  ! >> note efficieny <<
       

!... schumann-runge band

      fact1=sqrt(1.0+1.734e-22*sno2)     ! 4*sigma/(pi*y) = 4*2.07e-24/(pi*0.0152)
      
      fact2=0.023876*(fact1-1.0)         ! pi*y/2 = (pi*0.0152)/2 = 0.023876

      htsrb= 2.6496e-26 * exp( -amin1(2.0e2,fact2) ) / fact1
                                         ! f*sigma=0.0128*2.07e-24=2.6496e-26
      htsrb=htsrb*o2

!... heating rate due to diffuse scatted solar radiation

      tau  = 2.85e-25 * bo3
      taus = 2.85e-25 * bo3t

      fac1 = taus/bmu + 1.9*(taus-tau)

      htscat = 2.109e-22 * o3 * ome0 * bmu * exp(-amin1(2.0e2, fac1))
                                          ! 2.11e-22= 2*370*2.85e-25

!----------------------------------------------------------------------
!... total heating rate
!----------------------------------------------------------------------

      htmu = htsum1 + htsum2 + htsrc + htsrb + htscat

      return
      end
      
      
!----------------------------------------------------------------------
      function effcy(p, mode)
!----------------------------------------------------------------------
!
!  calculates solar heating efficiency factors based on Mlynczak
!  and Solomon (1993)
!
!----------------------------------------------------------------------

      implicit none

      real :: effcy    
      real :: p         ! pressure (mbar) 
      integer :: mode   ! (mode=3) efficiency of o3 hartley band 
                        ! (mode=2) efficiency of o2 sr continuum 
      
      real :: pmb
      real :: c0, c1, c2, c3
      real :: xp, x, x2, x3 
      real :: effcy1

      pmb = p

      if (pmb.gt.1.0) pmb=1.0
      if (pmb.lt.1.0e-4) pmb=1.0e-4

      xp=alog10(pmb)

!----------------------------------------------------------------------
!... ozone efficiency coefficients
!----------------------------------------------------------------------

      if(mode.eq.3) then
        if(pmb.le.1.0e-2) then
	x=xp+3.0
	  c0=0.66965
	  c1=-0.009682
	  c2=0.033093
	  c3=0.017938
	else
	  x=xp+1.0
	  c0=0.92621
	  c1=0.13396
	  c2=-0.076863
	  c3=0.006897
	endif
      endif

!----------------------------------------------------------------------
!... o2 efficiency coefficients 
!----------------------------------------------------------------------

      if(mode.eq.2) then
        if(pmb.le.1.0e-2) then
	  x=xp+3.0
	  c0=0.75349
	  c1=0.0036
	  c2=0.059468
	  c3=-0.022795
	else
	  x=xp+1.0
	  c0=0.92621
	  c1=0.13396
	  c2=-0.076863
	  c3=0.006897
	endif
      endif

!----------------------------------------------------------------------
!... calculate efficiency using polynomial coeffs 
!----------------------------------------------------------------------

      x2=x*x
      x3=x2*x
      effcy1 = c0 + c1*x + c2*x2 + c3*x3

      if (mode.eq.3) effcy = 1.0 - 0.78831 * (1.0 - effcy1)

      if (mode.eq.2) effcy = 1.0 - 0.41 * (1.0 - effcy1)

      if ((effcy .lt. 0) .or. (effcy .gt. 1.0)) then
         print *, 'error in effcy', effcy
	 stop
      endif

      return
      end

!----------------------------------------------------------------------

