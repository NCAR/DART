c-------------------------------------------------------------------
      subroutine tinterp_dyn (dayofyr)
c-------------------------------------------------------------------
c  interpolate climatological data in time
c  u, v, t, phi lower boundary
c-------------------------------------------------------------------

      use params
      use dynam
      use phys

      use utilities_mod, only : open_file, close_file

      implicit none

      integer, intent(in) :: dayofyr

      integer, parameter :: nmonths = 12
      integer :: iunit
      integer :: j, jj, k, kk, ll, m, mm, m1, m2
      integer :: iday, nfirst(nmonths), icall, mon, month, nbar
      real :: ddays(nmonths), f1, f2
      real, dimension(ny,nmonths) :: tuin, uuin
      real, dimension(nz,ny,nmonths) :: ubarin

      real :: tinit(nzp1,ny),uinit(nzp1,ny),phbot(ny),tupperbc(ny)

      save tuin, uuin

      data nfirst/15,46,74,105,135,166,196,226,258,288,319,349/
      data ddays/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./
      data icall/0/

c-------------------------------------------------------------------
c  read all input files once and save for interpolation
c-------------------------------------------------------------------
      icall = icall+1
      if(icall.le.1)then

c  read zonal mean climatology

         iunit = open_file(namf12,form='formatted',action='read')     
 
c        open (unit = 12,
c    $         file = namf12,
c    $         form = 'formatted',
c    $         status = 'old')

         do mm = 1,nmonths
            read (iunit,'(i2)'    ) month
            read (iunit,'(16f5.1)') uinit
            read (iunit,'(16f5.0)') tinit
            read (iunit,'(10f8.0)') phbot
            read (iunit,'(16f5.0)') tupperbc
            do j=1,ny
               do k=1,nz
                  ubarin(k,j,mm) = uinit(k+1,j)
               end do
            end do
         end do
         call close_file(iunit)

c  read zonal mean boundary data
c  upper dynamical

         iunit = open_file(namf23,form='formatted',action='read')     

c        open (unit = 23,
c    $         file = namf23,
c    $         form = 'formatted',
c    $         status = 'old')

         do mm=1,nmonths
            read(iunit,'(i2)')mon
            if(mon.ne.mm)then
               print '(1x,"zonal ubc read error",2i6)',mon,mm
            end if
            read(iunit,'(16f5.1)')(uuin(ll,mm),ll=1,ny)
            read(iunit,'(16f5.0)')(tuin(ll,mm),ll=1,ny)
         end do
         call close_file(iunit)
      end if
c______________________________________________________________________
c______________________________________________________________________

c  Now interpolate all data to day of year
c  Also, convert temperature to perturbation temperature

c______________________________________________________________________
c______________________________________________________________________

c  interpolation parameters
      iday = dayofyr
      if(iday.le.nfirst(1))then
         m1 = nmonths
         m2 = 1
         f1 = float(15-iday)/ddays(1)
      end if
      if(iday.gt.nfirst(nmonths))then
         m1 = nmonths
         m2 = 1
         f1 = float(380-iday)/ddays(nmonths)
      end if
      if(iday.gt.nfirst(1).and.iday.le.nfirst(nmonths))then
         do mm=1,11
            if(iday.gt.nfirst(mm).and.iday.le.nfirst(mm+1))then
               m1 = mm
               m2 = mm+1
               f1 = float(nfirst(m2)-iday)/ddays(m1)
            end if
         end do
      end if
      f2 = 1.-f1

      print *, 'tinterp - dayofyr: ', dayofyr
      print *, 'tinterp - f1, f2: ', f1, f2 

c______________________________________________________________________
c______________________________________________________________________

c  fields returned to model

c______________________________________________________________________
c  climatology
      do j=1,ny
         do k=1,nz
            uclimo(k,j) = (f1*ubarin(k,j,m1) + f2*ubarin(k,j,m2))
         end do
      end do

c______________________________________________________________________
c  dynamical upper boundary

      do j=1,ny
         unzp1(j) = (f1*uuin(j,m1) + f2*uuin(j,m2))
         tnzp1(j) = (f1*tuin(j,m1) + f2*tuin(j,m2))         !  full T
      end do

      end
