      subroutine rates (tfull)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     This routine calculates the chemical rate coefficients          c
c     for gas phase reactions                                         c
c     Notation from Brasseur and Solomon, 1984                        c
c     rates from JPL, 1994                                            c
c                                                                     c
c     no ozone hole reaction rates                                    c
c----------------------------------------------------------------------
c
c...hox chemistry rates
c           a1:    h + o2 + m  -> ho2 + m
c           a1et:  o1d + h20   -> oh + oh
c           a2:    h + o3      -> oh + o2
c           a3et:  o1d + h2    -> oh + h
c           a5:    o + oh      -> o2 + h
c           a6:    oh + o3     -> ho2 + o2
c           a6b:   ho2 + o3    -> oh + 2o2
c           a7:    o + ho2     -> oh + o2
c           a17:   oh + ho2    -> h2o + o2
c           a19:   oh + h2     -> h2o + h
c           a23a:  h + ho2     -> oh + oh
c           a23b:  h + ho2     -> h2 + o2
c           a23c:  h + ho2     -> h2o + o
c           a26:   no + ho2    -> no2 + oh
c           a27:   ho2 + ho2   -> h2o2 + o2
c           a30:   oh + h2o2   -> h2o + ho2
c           a36:   oh + co     -> co2 + h
c           a81:   h2o2 + o    -> oh + ho2
c           a82:   oh + oh     -> o + h2o
c           a83:   oh + oh + m -> h2o2 + m 

c...nox
c           b3 :   o + no2     -> no + o2
c           b4 :   o3 + no     -> no2 + o2
c           b6 :   n + no      -> n2 + o
c           b7 :   n + o2      -> no + o
c           b7a:   n(2d) + o2  -> no + o
c           b9 :   o3 + no2    -> no3 + o2
c           b12:   no2 + no3 + m -> n2o5 + m
c           b32:   n2o5 + m    -> no2 + no3 + m
c           b22:   oh + no2 + m -> hno3 + m
c           b23:   ho2 + no2 + m -> ho2no2 + m
c           b24:   ho2no2 + m  -> ho2 + no2 + m
c           b27:   hno3 + oh   -> h2o + no3
c           b28:   oh + ho2no2 -> products
c           b38:   o1d + n2o   -> n2 + o2
c           b39:   o1d + n2o   -> no + no
c           b71:   o + no3     -> o2 + no2
c           b72:   oh + no3    -> ho2 + no2
c           b73a:  ho2 + no3   -> oh + no2 + o2
c           b73b:  ho2 + no3   -> hno3 + o
c           b81:   o + no2 + m -> no3 + m
c           b82:   no + o + m  -> no2 + m
c           b84:   no + no3    -> no2 + no2

c...hydrocarbon
c           c1:    ch4 + o1d   -> ch3 + oh
c           c1et:  o1d + ch4   -> oh + ch3
c           c2:    ch4 + oh    -> ch3 + h2o
c           c8:    ch2o + oh   -> cho + h2o
c           c9:    ch2o + o    -> cho + oh

c...clx
c           d2:    cl + o3      -> clo + o2
c           d3:    clo + o      -> cl + o2
c           d4:    clo + no     -> no2 + cl
c           d5:    cl + ch4     -> hcl + ch3
c           d6:    cl + h2      -> hcl + h
c           d7:    cl + ho2     -> hcl + o2
c           d8:    clo + oh     -> cl + ho2
c           d10:   ch2o + cl     -> hcl + hco
c           d11:   oh + hcl     -> h20 + cl
c           d31:   clo + no2 + m -> clono2 + m
c           d32:   o + clono2   -> products
c           d33:   clo + ho2    -> hocl + o2
c           d34:   oh + hocl    -> h20 + clo
c           d35:   o + hocl     -> oh + clo
c        *  d36:   cl + no2 + m -> clno2 + m
c        *  d37:   cl + hocl    -> oh + cl2
c           d46:   clo + oh     -> hcl + o2
c        *  d47:   clo + clo    -> cl + oclo
c        *  d48:   clo + clo    -> cl2 + o2
c        *  d60:   clo + clo + m -> cl2o2 + m
c        *  d61:   cl2o2 + m    -> clo + clo + m
c        *  d62:   oclo + oh    -> hocl + o2
c        *  d63:   cl + oclo    -> clo + clo
c        *  d64:   oclo + o     -> clo + o2
c        *  d65:   oclo + no    -> no2 + clo
c        *  d71:   cl2 + o1d    -> cl + clo
c        *  d72:   cl2o2 + cl   -> cl2 + clo2
c           d73:   no3 + cl     -> clo + no2
c        *  d74:   clo + no3    -> no2 + clo2
c           d75:   hcl + o1d    -> cl + oh
c        *  d81:   cl2 + oh     -> hocl + cl
c        *  d82:   cl + clono2  -> no3 + cl2
c           d83:   ho2 + cl     -> clo + oh
c           d84:   h2o2 + cl    -> hcl + ho2
c           d85:   hcl + o      -> cl + oh
c           d87:   clono2 + oh  -> hocl + no3

c...brx
c        *  e2 :   br + o3    -> bro + o2
c        *  e3 :   bro + o    -> br + o2
c        *  e4 :   bro + no   -> no2 + br
c        *  e5a:   bro + clo  -> oclo + br
c        *  e5b:   bro + clo  -> br + cl + o2
c        *  e5c:   bro + clo  -> brcl + o2
c        *  e6 :   bro + bro  -> 2br + o2
c        *  e7 :   br + ho2   -> hbr + o2
c        *  e8 :   br + oclo  -> bro + clo
c        *  e9 :   br + ch2o  -> hbr + hco
c        *  e11:   oh + hbr   -> h2o + br
c        *  e13:   bro + no2 + m -> brono2 + m
c        *  e15:   bro + ho2  -> hobr + o2
c        *  e71:   hbr + o1d  -> oh + br
c        *  e72:   oh + bro   -> ho2 + br
c        *  e81:   hbr + o    -> br + oh
c
c...ox 
c           hk1:   o + o + m  -> o2 + m
c           hk2:   o + o2 + m -> o3 + m
c           hk3:   o + o3     -> o2 + o2
c           hk4:   o1d + n2   -> o + n2
c           hk5:   o1d + o2   -> o + o2
c           hk7:   o1d + o3   -> 2 o2
c           hk21:  o1d + n2   -> n2o
c
c
c...heterogeneous reactions - calculated in subroutine polchem
c           het1:   clono2 + h2o -> hocl + hno3
c           het2:   clono2 + hcl -> cl2  + hno3
c           het3:   n2o5   + h2o ->       2hno3
c           het4:   n2o5   + hcl -> clno2 +hno3
c           het5:   hocl   + hcl -> cl2 + h2o  
c
c
c...aerosol reactions - calculated in subroutine aerosols
c           g1:   clono2 + h2o ->  hocl + hno3
c           g2:   n2o5 + h2o   ->  2 hno3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      use params
      use chem
      use dynam

      implicit none

      integer :: i, j, k
      real tfull(nz,nx,ny), xlnt(nz,nx,ny), xpo(nz,nx,ny), aux(nz,nx,ny)
      real pm(nz,nx,ny), tinv(nz,nx,ny), ak0(nz,nx,ny), ak1(nz,nx,ny)
      real :: xln6, a0hox,  a1hox, a0hox1, a1hox1, a0nox1, a1nox1, 
     $        a0nox2, a1nox2, a0nox3, a1nox3, a0nox4, a1nox4,
     $        a0nox5, a1nox5, a0clx1, a1clx1, a0clx2, a1clx2,
     $        a0clx3, a1clx3,  a0brx,  a1brx

      logical firstcall

      save xln6, a0hox,  a1hox, a0hox1, a1hox1, a0nox1, a1nox1, 
     $          a0nox2, a1nox2, a0nox3, a1nox3, a0nox4, a1nox4,
     $          a0nox5, a1nox5, a0clx1, a1clx1, a0clx2, a1clx2,
     $          a0clx3, a1clx3,  a0brx,  a1brx,
     $          pm

      data firstcall /.true./

      if (firstcall)then
         firstcall = .false.
	 print *, '*****************rates*****************'
c----------------------------------------------------------------------
c
c        Rate constants which are not temperature dependent
c
c----------------------------------------------------------------------
c
c...ox
         hk7 = 1.2e-10
c
c----------------------------------------------------------------------
c
c...hox
         a1et = 2.2e-10

         a3et = 1.0e-10

         a23a = 8.1e-11*0.9   

         a23b = 8.1e-11*0.08   

         a23c = 8.1e-11*0.02  

         a23  = a23a + a23b + a23c

         a23p = a23a - (a23b + a23c)

         a0hox = 5.7e-32*300.**1.6
         a1hox = 7.5e-11
         a0hox1 = 6.9e-31*300.**0.8
         a1hox1 = 1.5e-11
c
c----------------------------------------------------------------------
c
c...nox
         b7a = 5.0e-12

         b38 = 4.9e-11

         b39 = 6.7e-11

         b71 = 1.0e-11

         b72 = 2.2e-11

         b73a = 3.5e-12*0.7

         b73b = 3.5e-12*0.3

         a0nox1 = 2.2e-30*300.**3.9
         a1nox1 = 1.5e-12*300.**0.7
         a0nox2 = 2.6e-30*300.**3.2
         a1nox2 = 2.4e-11*300.**1.3
         a0nox3 = 1.8e-31*300.**3.2
         a1nox3 = 4.7e-12*300.**1.4
         a0nox4 = 9.0e-32*300.**2.0
         a1nox4 = 2.2e-11 
         a0nox5 = 9.0e-32*300.**1.5
         a1nox5 = 3.0e-11
c
c----------------------------------------------------------------------
c
c...hydrocarbons
         c1 = 1.5e-10

         c8 = 1.e-11
c
c----------------------------------------------------------------------
c
c...clx
         d35 = 1.7e-13                        ! new for JPL-97 

         d73 = 2.4e-11

         d75 = 1.5e-10

         a0clx1 = 1.8e-31*300.**3.4
         a1clx1 = 1.5e-11*300.**1.9
c
c----------------------------------------------------------------------
c
c...other constants
         xln6 = alog(.6)

         do j=1,ny
            do k = 1,nz
               do i = 1,nx
                  pm(k,i,j) = pmb(k)
               end do
            end do
         end do
c----------------------------------------------------------------------

      end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c----------------------------------------------------------------------
c
c  temperature and/or density dependent rate coefficients
c
c----------------------------------------------------------------------

      xlnt = alog(tfull)
      tinv = 1./tfull
      where (tfull <= 150.) 
         tinv = 1./150.
      end where

c----------------------------------------------------------------------
c...hox chemistry rates
      ak0 = a0hox * exp(xlnt*(-1.6))
      xpo = 1./(1. + alog10((ak0 * hnm)/a1hox)**2)
      a1 = (ak0 * hnm) /(1. + ak0 * hnm/a1hox) * exp(xln6*xpo)

      a2 = 1.4e-10 * exp(-470.*tinv)

      a5 = 2.2e-11 * exp(120.*tinv)

      a6 = 1.5e-12 * exp(-880.*tinv)          ! new for JPL00

      a6b = 2.0e-14 * exp(-680.*tinv)         ! new for JPL00

      a7 = 3.0e-11 * exp(200.*tinv)

      a17 = 4.8e-11 * exp(250.*tinv)

      a19 = 5.5e-12 * exp(-2000.*tinv)

      a26 = 3.5e-12 * exp(250.*tinv)

      a27 = (2.3e-13 * exp(600.*tinv)
     $    + 1.7e-33 * exp(1000.*tinv) * hnm)
     $              * (1. + 1.4e-21 * (qn1(:,:,:,3)*hnm) 
     $              * exp(2200.*tinv))

      a30 = 2.9e-12 * exp(-160.*tinv)

      a36 = 1.5e-13 * (1. + .6*pm/1013.)

      a81 = 1.4e-12 * exp(-2000.*tinv)

      a82 = 4.2e-12 * exp(-240.*tinv)

      ak0 = a0hox1 * exp(xlnt*(-0.8))
      xpo = 1./(1. + alog10((ak0 * hnm)/a1hox1)**2)
      a83 = (ak0 * hnm) /(1. + ak0 * hnm/a1hox1) * exp(xln6*xpo)

c
c----------------------------------------------------------------------
c...the nox chemistry reaction rates

      b3 = 5.6e-12 * exp(180.*tinv)           ! new for JPL00
            
      b4 = 3.e-12 * exp(-1500.*tinv)          ! new for JPL00

      b6 = 2.1e-11 * exp(100.*tinv)

      b7 = 1.5e-11 * exp(-3600.*tinv)

      b9 = 1.2e-13 * exp(-2450.*tinv)

      ak0 = a0nox1 * exp(xlnt*(-3.9))
      ak1 = a1nox1 * exp(xlnt*(-.7))
      xpo = 1./(1. + alog10((ak0 * hnm)/ak1)**2)
      b12 = (ak0 * hnm) /(1. + ak0 * hnm/ak1) * exp(xln6*xpo)

      ak0 = a0nox2 * exp(xlnt*(-3.2))
      ak1 = a1nox2 * exp(xlnt*(-1.3))
      xpo = 1./(1. + alog10((ak0 * hnm)/ak1)**2)
      b22 = (ak0 * hnm) /(1. + ak0*hnm/ak1) * exp(xln6*xpo)

      ak0 = a0nox3 * exp(xlnt*(-3.2))
      ak1 = a1nox3 * exp(xlnt*(-1.4))
      xpo = 1./(1. + alog10((ak0*hnm)/ak1)**2)
      b23 = (ak0 * hnm)/(1. + ak0 * hnm/ak1) * exp(xln6*xpo)

      b24 = b23 / (2.1e-27 * exp(10900.*tinv))

      aux = 1.9e-33 * exp(725.*tinv) * hnm
      b27 = 7.2e-15 * exp(785.*tinv)
     $    + aux / (1. + aux/(4.1e-16 * exp(1440.*tinv)))

      b28 = 1.3e-12 * exp(380.*tinv)

      b32 = b12 / (2.7e-27 * exp((11000.*tinv)))

      ak0 = a0nox4 * exp(xlnt*(-2.0))
      ak1 = a1nox4
      xpo = 1./(1. + alog10((ak0 * hnm)/ak1)**2)
      b81 = (ak0 * hnm) /(1. + ak0 * hnm/ak1) * exp(xln6*xpo)

      ak0 = a0nox5 * exp(xlnt*(-1.5))
      ak1 = a1nox5
      xpo = 1./(1. + alog10((ak0 * hnm)/ak1)**2)
      b82 = (ak0 * hnm) /(1. + ak0 * hnm/ak1) * exp(xln6*xpo)

      b84 = 1.5e-11 * exp(170.*tinv)

c----------------------------------------------------------------------
c...hydrocarbon reaction rates

      c2 = 2.65e-12 * exp(-1800.*tinv)

      c9 = 3.40e-11 * exp(-1600.*tinv)

c----------------------------------------------------------------------
c...the clx chemistry reaction rates

      d2 = 2.3e-11 * exp(-200.*tinv)          ! new for JPL00

      d3 = 3.e-11 * exp(70.*tinv)

      d4 = 6.4e-12 * exp(290.*tinv)

      d5 = 9.6e-12 * exp(-1360.*tinv)         ! new for JPL00

      d6 = 3.7e-11 * exp(-2300.*tinv)

      d7 = 1.8e-11 * exp(170.*tinv)

      d8 = 7.4e-12 * exp(270.*tinv)           ! new for JPL00

      d10 = 8.1e-11 * exp(-30.*tinv)

      d11 = 2.6e-12 * exp(-350.*tinv)
 
      ak0 = a0clx1 * exp(xlnt*(-3.4))
      ak1 = a1clx1 * exp(xlnt*(-1.9))
      xpo = 1./(1. + alog10((ak0 * hnm)/ak1)**2)
      d31 = (ak0 * hnm) /(1. + ak0 * hnm/ak1) * exp(xln6*xpo)

      d32 = 2.9e-12 * exp(-800.*tinv)

      d33 = 4.8e-13 * exp(700.*tinv)

      d34 = 3.0e-12 * exp(-500.*tinv)

      d46 = 3.2e-13 * exp(320.*tinv)          ! new for JPL00

      d83 = 4.1e-11 * exp(-450.*tinv)

      d84 = 1.1e-11 * exp(-980.*tinv)

      d85 = 1.0e-11 * exp(-3300.*tinv)

      d87 = 1.2e-12 * exp(-330.*tinv)

c----------------------------------------------------------------------
c...the ox chemistry reaction rates

      hk1 = 4.23e-28 * hnm/(tfull * tfull)

      hk2 = 6.0e-34 * (300.*tinv)**2.4 * hnm  ! new for JPL00

      hk3 = 8.e-12 * exp(-2060.*tinv)
            
!      hk4 = 1.8e-11 * exp(110.*tinv)
      hk4 = 2.1e-11 * exp(115.*tinv)    ! new from Ravishankara et al., 2002
            
      hk5 = 3.2e-11 * exp(70.*tinv)

      hk21 = 3.5e-37 * (300.*tinv)**0.6 * hnm

c  initialize aerosol rates
c      g1  = 0.
c      g2  = 0.

      return
      end










