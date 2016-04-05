      subroutine zenith(dayofyr, gmt_frac)

!---------------------------------------------------------------------
!     calculates at each longitude and latitude the secant of solar
!     zenith angle (sc2d) for a given day of year and GMT
!---------------------------------------------------------------------

      use params, only: nx, ny
      use dynam,  only: pi, cosfi, sinfi
      use chem,   only: sc2d    ! secant of solar  zenith angle

      implicit none

      integer :: dayofyr
      real :: gmt_frac

!... local variables

      real :: rnx

      real :: theta0            ! day of year in radians
      real :: ahor              ! hour angle
      real :: cosc              ! cosine of sza

      real :: decl              ! solar declination
      real :: cosdec, sindec    ! cos and sin of declination
      real :: coscos, sinsin
      integer :: i, j

!---------------------------------------------------------------------
!     calculation of solar declination 
!---------------------------------------------------------------------
!  A series expansion is used to calculate declination. Taken from
!  "Radiative Processes in Meteorology and Climatology" by Paltridge
!  and Platt (Elsevier pub.). This replace the simpler formula below, 
!  which ignores orbital eccentricity: 
!  
!      real :: day_time
!      day_time = float(dayofyr) + gmt_frac
!      decl = -pi/180.*23.5*sin(2.*pi*(day_time - 0.5 - 264.)/365.)
!---------------------------------------------------------------------

      theta0 = 2. * pi * (float(dayofyr-1) + gmt_frac) / 365.

      decl = 0.006918 - 0.399912*cos(theta0) + 0.070257*sin(theta0)
     $     - 0.006758*cos(2*theta0) + 0.000907*sin(2*theta0)
     $     - 0.002697*cos(3*theta0) + 0.001480*sin(3*theta0)
      
      print *, 'zenith, decl:', decl

      rnx = 1./float(nx)

!---------------------------------------------------------------------
!     calculation of secant of solar zenith angle
!---------------------------------------------------------------------

      cosdec = cos(decl)
      sindec = sin(decl)

      do j=1,ny
         coscos = cosdec*cosfi(j)
         sinsin = sindec*sinfi(j)
         do i = 1,nx
            ahor = 2.*pi*(gmt_frac - 0.5 + (float(i) - 1.)*rnx)

            cosc = amax1( amin1( cos(ahor)*coscos + sinsin,
     $                    0.9999999 ), -0.9999999 )

            sc2d(i,j) = 1./amax1(cosc, 1.8e-2)
         end do
      end do

      end

!---------------------------------------------------------------------

