c-----------------------------------------------------------------------
      subroutine phot(dayofyr)
c-----------------------------------------------------------------------
c
c    calculates photodissociation frequencies (J coefficients) [sec-1]
c    by loglinear interpolation 
c
c-----------------------------------------------------------------------
c
c    input : pmb     = pressure ( in millibars )
c            sc2d    = secant of the zenith angle at surface
c            o3t     = ozone column
c            dayofyr = used to correct for Sun-Earth distance
c
c    output: tj      = photodiss. freq. (for nphot species)
c
c    CAUTION: The interpolation program is turned off for some
c             species, so tj is not completely filled.  
c             Change values of "index" to turn on/off calculation
c             of specific J coefficients.
c
c    A warning message will appear if the ozone column is greater than 
c    the maximum or less than the minimum standard column used for the 
c    table.  To recompute the table with different standard columns, see 
c    ../phot/tuv/jtable.mod.f or ../phot/stam/j.stam.mlt.f
c
c    The table can also be modified to include additional species, 
c    such as halocarbons.
c
c    The photolysis table automatically includes the Chapman function 
c    to account for dependence of solar zenith angle on altitude.
c
c-----------------------------------------------------------------------

      use params

      use chem, only: tj, sc2d, o3t 
      use dynam, only: namf30, pi

      use utilities_mod, only : open_file, close_file

      implicit none

      integer, intent(in) :: dayofyr

      integer, parameter :: nzen=10         ! number of zenith angles
      integer, parameter :: ncol=15         ! number of ozone columns

      integer :: iunit

      real :: v3std(nz)                     ! standard ozone column
      real :: ajl(nphot, nz, nzen, ncol)    ! look-up table

      real :: v3rat(nz)                     ! ratio of model ozone column to standard

      real :: dd2, dd3                      ! variables used in bi-linear
      real :: dels(4), weight(2,2)          ! interpolation
      integer :: is, iv(nz)
      real :: g, dist                       ! used to correct for Sun-Earth distance

      logical :: entered = .false. 
      real :: ln10,  sum_ajl
      integer :: i, j, k, l
      integer :: ll, mm, nn

      real :: vsec(nzen)  = (/ 1., 1.064, 1.155, 1.305, 1.556,
     $            2., 2.924, 5.759, 11.474, 50. /)

      real :: xv3(ncol) = (/ .04, .1, .2, .3, .5, .75, 1., 1.25, 
     $          1.5, 1.75, 2., 3., 4., 6., 12. /)

      integer :: index(nphot) = (/
     $          1, 1, 1, 1,  ! o2,          o2->o1d,     o3,         o3->o1d
     $          1, 1, 1, 1,  ! ch4,         no2,         hno3,       hocl 
     $          1, 1, 1, 1,  ! ho2no2->ho2, ho2no2->oh,  n2o5->no2,  n2o5->no,
     $          1, 0, 0, 1,  ! h2o2,        oclo,        cl2o2,      hcl 
     $          0, 1, 1, 0,  ! cl2,         co2,         n2o,        brono2
     $          0, 0, 1, 1,  ! brcl,        hobr,        hno2,       ch2o->2h
     $          1, 1, 1, 1,  ! ch2o->h2,    clono2->clo, clono2->cl, no  
     $          1, 1, 1, 1,  ! no3->no2,    no3->no,     h2o->h+oh,  h2o->h2+o(1d) 
     $          1  /)        ! h2o->2h+o

      save   entered, ajl, v3std, ln10

c--------------------------------------------------------------------
c     first call - read in look-up table
c--------------------------------------------------------------------

      if( .not. entered ) then

         print *, 'reading '//namf30

         iunit = open_file(namf30,form='formatted',action='read')         
  
c        open (unit = 30,
c    $         file = namf30,
c    $         form = 'formatted')

         read (iunit,'(6e13.5)') v3std     ! read in standard ozone column

         do i = 1,nphot                     ! read in look-up table 
            do l = 1,ncol
               do k = 1,nzen
                  read(iunit,'(8e12.4)') (ajl(i, j, k, l), j=1, 8)
                  read(iunit,'(8e12.4)') (ajl(i, j, k, l), j=9, 16)
                  read(iunit,'(8e12.4)') (ajl(i, j, k, l), j=17, 24)
                  read(iunit,'(8e12.4)') (ajl(i, j, k, l), j=25, 32)
                  read(iunit,'(8e12.4)') (ajl(i, j, k, l), j=33, 38)
               end do
            end do
         end do

         call close_file(iunit)

         ln10 = alog(10.)
         entered = .true.

      end if

      tj = 0.0
      
c----------------------------------------------------------------------
c  loop over latitude (j) and longitude (k)
c----------------------------------------------------------------------

      do j=1,ny
        do k = 1,nx

c... calculate only selected J's in daytime

          if (sc2d(k,j) .le. 50) then 

            if( sc2d(k,j) .lt. vsec(1)-.01 ) then
              print 195, j, k, sc2d(k,j)
 195          format('zenith angle out of bounds',2i3,f10.3)
              stop
            end if
   
c... find zenith angle in table (altitude independent)

            do mm = 1,nzen
              if( vsec(mm) .gt. sc2d(k,j) ) goto 200
            end do
            mm = mm - 1
 200        continue
            is = max0( mm - 1, 1)

c... find ozone column in table
            do i = 1,nz
               v3rat(i) = o3t(i,k,j)/v3std(i)
	       
               if (v3rat(i).lt.xv3(1) .or. v3rat(i).gt.xv3(ncol)) then
c                  print 210,i,j,k,v3rat(i),v3std(i),o3t(i,k,j)
 210              format('ozone column out of bounds',3i3,f10.3,2e15.5)
               end if
	       
               do mm = 1,ncol
                  if( xv3(mm) .gt. v3rat(i) ) goto 300
               end do
               mm = mm - 1
 300           continue
               iv(i) = max0( mm - 1, 1)
            end do

c... calculate weights for bi-linear interpolation
 
            do i = 1,nz
               ll = iv(i)
               dels(2) = (sc2d(k,j) - vsec(is))/(vsec(is+1) 
     $                       - vsec(is))
               dels(3) = (v3rat(i) - xv3(ll))/(xv3(ll+1) - xv3(ll))

c...  do not allow extrapolation of ozone column ratio

               if(dels(3).gt.1.) dels(3) = 1.
               if(dels(3).lt.0.) dels(3) = 0.

               dd2 = 1. - dels(2)
               dd3 = 1. - dels(3)
               weight(1,1) = dd2*dd3
               weight(1,2) = dd2*dels(3)
               weight(2,1) = dels(2)*dd3
               weight(2,2) = dels(2)*dels(3)

c... interpolate J's from table

               do  nn = 1,nphot
	         if (index(nn).eq.1) then
                   sum_ajl = weight(1,1) * ajl(nn, i, is, ll)
     $                     + weight(1,2) * ajl(nn, i, is, ll+1)
     $                     + weight(2,1) * ajl(nn, i, is+1, ll)
     $                     + weight(2,2) * ajl(nn, i, is+1, ll+1)
                   tj(i,k,j,nn) = exp( ln10*sum_ajl )
		 endif
               end do
            end do

          else                               ! zero J's at night

            do i = 1,nz
               do  nn = 1,nphot
                 tj(i,k,j,nn) = 0.
               end do
            end do

	  endif

        end do                               ! end longitude loop
      end do                                 ! end latitude loop

!... correct J's for variation with season of distance to the Sun
!    formula from Paltridge & Platt
!    dist = (Rbar/R)**2

      g = 2. * pi * dayofyr / 365.
      dist = 1.000110 + 0.034221*cos(g)   + 0.001280*sin(g) 
     $                + 0.000719*cos(2*g) + 0.000077*sin(2*g)

!      print *, 'Sun-Earth distance correction:', dist

      tj = tj*dist
      
      return
      end
