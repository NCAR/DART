c-------------------------------------------------------------------
      subroutine tinterp_rad (dayofyr)
c-------------------------------------------------------------------
c  interpolate climatological data in time
c  ozone & other gase for heating/cooling
c  1  o3
c  2  o2
c  3  h2o
c  4  o
c  5  co2
c-------------------------------------------------------------------

      use params
      use dynam
      use phys

      use utilities_mod, only : open_file, close_file

      implicit none

      integer, intent(in) :: dayofyr

      integer :: iunit
      integer :: j, jj, k, kk, ll, m, mm, m1, m2
      integer :: iday, icall, mon, month, nbar
      real :: f1, f2

      integer, dimension (16) :: rd_dates
      real, dimension (nzz,ny,16,5) :: rad_g

      save rd_dates, rad_g

      data icall/0/

c-------------------------------------------------------------------
c  read all input files once and save for interpolation
c-------------------------------------------------------------------
      icall = icall+1
      if(icall.le.1)then

c  read water data

         iunit = open_file(namf28,form='unformatted',action='read')

c        open (unit = 28,
c    $        file = namf28,
c    $        form = 'unformatted',
c    $         status = 'old')

         read(iunit) rd_dates, rad_g

         call close_file(iunit)

         print *,rd_dates

      end if
c______________________________________________________________________
c______________________________________________________________________

c  Now interpolate all data to day of year

c______________________________________________________________________
c______________________________________________________________________

c  interpolation parameters
      iday = dayofyr
      if(iday.eq.rd_dates(1))then
         m1 = 1
         m2 = 1
         f1 = 1.
      end if
      if(iday.gt.rd_dates(1))then
         do mm=1,15
            if(iday.gt.rd_dates(mm).and.iday.le.rd_dates(mm+1))then
               m1 = mm
               m2 = mm+1
               f1 = float(rd_dates(m2)-iday)/25.
            end if
         end do
      end if
      f2 = 1.-f1

      print *, 'tinterp - dayofyr: ', dayofyr
      print *, 'tinterp - f1, f2: ', f1, f2 

c______________________________________________________________________
c______________________________________________________________________

c  fields returned to model

      do j=1,ny
         do k=1,nzz
            cl_o3zz(k,j) = (f1*rad_g(k,j,m1,1) + f2*rad_g(k,j,m2,1))
            cl_h2o(k,j)  = (f1*rad_g(k,j,m1,3) + f2*rad_g(k,j,m2,3))
            cl_o3p(k,j)  = (f1*rad_g(k,j,m1,4) + f2*rad_g(k,j,m2,4))
            cl_co2(k,j)  = (f1*rad_g(k,j,m1,5) + f2*rad_g(k,j,m2,5)) 
         end do
         do k=1,nz
            cl_o3(k,j) = cl_o3zz(k+nztrop,j)
            cl_o2(k,j) = (f1*rad_g(k,j,m1,2) + f2*rad_g(k,j,m2,2))
         end do
      end do


      end
