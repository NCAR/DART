
      subroutine moldiff

!---------------------------------------------------------------------
!	... Compute molecular diffusion coefficients on a per
!	    species basis
!---------------------------------------------------------------------
      use params
      use chem
      use dynam
      use phys

      implicit none

!---------------------------------------------------------------------
!	... Parameters
!---------------------------------------------------------------------
      real, parameter    :: const = 1.52e14             ! xmz : m**2/s

!---------------------------------------------------------------------
!	... Local variables
!---------------------------------------------------------------------
      integer ::  i, j, k, l, m
      real :: hscal, bb, bb_h, alphat
      real, dimension(nz) :: dtdz, densty, Hair
      real :: w_h, w_h2, w_o3p, w_o2, w_o3, w_n2, 
     $        w_no, w_no2, w_n, w_co

!---------------------------------------------------------------------
!	... Molecular weight of long-lived and short-lived species.
!     n2o 44.     ch4 16.     h2o 18.     co 28.     hcl 36.45
!     cl2 70.9    hbr 80.9    co2 44.     h2 2.      hf  20.
!     o1d 16.     oh 17.      cl 35.45    o3p 16.    o2 32.
!     o3 48.      ho2 34.     no2 46.     no 30.     br 79.9
!     clo 51.45   bro 95.9    no3 62.     h 1.
!---------------------------------------------------------------------

      w_h = 1.
      w_h2 = 2.
      w_o3p = 16.
      w_o2 = 32.
      w_o3 = 48.
      w_n = 14.
      w_n2 = 28.
      w_no = 30.
      w_no2 = 46.
      w_co = 28.

      do k=1,nz
         wmole(k) = 0.
         do j=1,ny
            do i=1,nx
               wmole(k) = wmole(k) 
     $                  + w_n2  * (1.0 - qn1(k,i,j,26) - qn1(k,i,j,18))  
     $                  + w_o2  * qn1(k,i,j,26) 
     $                  + w_o3p * qn1(k,i,j,18)
     $                  + w_h2  * qn1(k,i,j,16)
     $                  + w_o3  * qn1(k,i,j,19)
     $                  + w_h   * qn1(k,i,j,7)
     $                  + w_no  * qn1(k,i,j,23)
     $                  + w_no2 * qn1(k,i,j,24)
            end do
         end do
         wmole(k) = wmole(k) / (float(ny) * float(nx))
      end do

!---------------------------------------------------------------------
!	... Initialize xmz and wdif
!---------------------------------------------------------------------
      wdif = 0.
      dmolch = 0.
      do k=2,nz-1
         dtdz(k) = (tref(k+1) - tref(k-1))/dz/2.
      end do
      dtdz(nz) = (tref(nz) - tref(nz-1))/dz

      do k=1,nz
         densty(k) = pmb(k) /(tref(k)*1.38e-19)
         Hair(k) = 8.314 * tref(k) / (wmole(k)*9.806)        ! km
      end do

      alphat = -.38             ! Thermal diff coeff for H

!---------------------------------------------------------------------
!	... Formula of xmz from Bank and Kockart (1972)
!
!   Diffusion coefficient Di = 1.52e18 cm^2/s * ( 1/mi + 1/m )^0.5 * T^0.5 / N 
!
!   mi = mass of species i (amu)
!   m  = mean molecular mass (amu)
!   T  = temperature (K)
!   N  = number density (cm^-3)
!
!---------------------------------------------------------------------

      do k = 22,nz
!---------------------------------------------------------------------
!	... Formula of xmz from Bank and Kockart (1972)
!---------------------------------------------------------------------
         dmolch(k,7) = const *               ! H
     $              SQRT( (1./wmole(k) + 1./w_h)*tref(k)) / densty(k)
         dmolch(k,16) = const *              ! H2
     $              SQRT( (1./wmole(k) + 1./w_h2)*tref(k)) / densty(k)
         dmolch(k,18) = const *              ! O
     $              SQRT( (1./wmole(k) + 1./w_o3p)*tref(k)) / densty(k)
         dmolch(k,19) = const *              ! O3
     $              SQRT( (1./wmole(k) + 1./w_o3)*tref(k)) / densty(k)
         dmolch(k,22) = const *              ! N
     $              SQRT( (1./wmole(k) + 1./w_n)*tref(k)) / densty(k)
         dmolch(k,23) = const *              ! NO
     $              SQRT( (1./wmole(k) + 1./w_no)*tref(k)) / densty(k)
         dmolch(k,24) = const *              ! NO2
     $              SQRT( (1./wmole(k) + 1./w_no2)*tref(k)) / densty(k)
         dmolch(k,9) = const *               ! CO
     $              SQRT( (1./wmole(k) + 1./w_co)*tref(k)) / densty(k)

!---------------------------------------------------------------------
!	... Molecular diffusion "velocity" equivalent
!---------------------------------------------------------------------

         hscal = 8.314 * tref(k) / (w_h*9.8)                   ! km
         bb_h = 1./hscal + 7. * alphat * dtdz(k) * 1000.0 
     $                            /(Hair(k) * tref(k))         ! 1/km
         wdif(k,7) = .001 * dmolch(k,7) * (1./Hair(k) - bb_h)  ! m/s 

         hscal = 8.314 * tref(k) / (w_h2*9.8)
         bb = 1./hscal
         wdif(k,16) = .001 * dmolch(k,16) * (1./Hair(k) - bb)

         hscal = 8.314 * tref(k) / (w_o3p*9.8)
         bb = 1./hscal
         wdif(k,18) = .001 * dmolch(k,18) * (1./Hair(k) - bb)

         hscal = 8.314 * tref(k) / (w_n*9.8)
         bb = 1./hscal
         wdif(k,22) = .001 * dmolch(k,22) * (1./Hair(k) - bb)

         hscal = 8.314 * tref(k) / (w_no*9.8)
         bb = 1./hscal
         wdif(k,23) = .001 * dmolch(k,23) * (1./Hair(k) - bb)

         hscal = 8.314 * tref(k) / (w_no2*9.8)
         bb = 1./hscal
         wdif(k,24) = .001 * dmolch(k,24) * (1./Hair(k) - bb)

         hscal = 8.314 * tref(k) / (w_co*9.8)
         bb = 1./hscal
         wdif(k,9) = .001 * dmolch(k,9) * (1./Hair(k) - bb)
      end do

      end subroutine MOLDIFF
