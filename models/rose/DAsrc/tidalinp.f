      subroutine tidalinp (dayofyr)
c  read tidal lower boundary fields and interpolate in time
      use params
      use dynam
      use phys

      implicit none

      integer :: m, m1, m2, mm, j, dayofyr, nfirst(4), montid, mindex
      integer icall, iday
      real :: f1, f2, g, ddays(4)
      real, dimension(ny,4) :: dauin, dbuin, davin, dbvin,
     $                         datin, dbtin, dazin, dbzin, 
     $                         sauin, sbuin, savin, sbvin,
     $                         satin, sbtin, sazin, sbzin
      save dauin, dbuin, davin, dbvin, datin, dbtin, dazin, dbzin, 
     $     sauin, sbuin, savin, sbvin, satin, sbtin, sazin, sbzin

      data nfirst/15,105,196,288/
      data ddays/90.,91.,92.,92./
      data icall/0/, g/9.8/

      icall = icall+1
      if(icall.le.1)then

c  tropospheric tidal forcing from GSWM
         open (unit = 52,
     $         file = namf52,
     $         form = 'formatted',
     $         status = 'old')
         do m=1,4
            read (52,200)montid
            mindex = 1 + (montid-1)/3
            if(m.ne.mindex) then
               print 100,m,montid,mindex
 100           format('problems with tidal b.c. input',3i5)
               stop
            end if
 200        format(i2)
            read(52,220) (dauin(j,m),j=1,ny),(dbuin(j,m),j=1,ny),
     $                   (davin(j,m),j=1,ny),(dbvin(j,m),j=1,ny),
     $                   (datin(j,m),j=1,ny),(dbtin(j,m),j=1,ny),
     $                   (dazin(j,m),j=1,ny),(dbzin(j,m),j=1,ny)
 220        format(8e10.3)
         end do
         close (52)
         open (unit = 53,
     $         file = namf53,
     $         form = 'formatted',
     $         status = 'old')
         do m=1,4
            read (53,200)montid
            mindex = 1 + (montid-1)/3
            read(53,220)
     $              (sauin(j,mindex),j=1,ny),(sbuin(j,mindex),j=1,ny),
     $              (savin(j,mindex),j=1,ny),(sbvin(j,mindex),j=1,ny),
     $              (satin(j,mindex),j=1,ny),(sbtin(j,mindex),j=1,ny),
     $              (sazin(j,mindex),j=1,ny),(sbzin(j,mindex),j=1,ny)
         end do
      end if

c  interpolate in time
      iday = dayofyr
      if(iday.le.nfirst(1))then
         m1 = 4
         m2 = 1
         f1 = float(15-iday)/ddays(1)
      end if
      if(iday.gt.nfirst(4))then
         m1 = 4
         m2 = 1
         f1 = float(380-iday)/ddays(4)
      end if
      if(iday.gt.nfirst(1).and.iday.le.nfirst(4))then
         do mm=1,3
            if(iday.gt.nfirst(mm).and.iday.le.nfirst(mm+1))then
               m1 = mm
               m2 = mm+1
               f1 = float(nfirst(m2)-iday)/ddays(m1)
            end if
         end do
      end if

      f2 = 1.-f1
      do j=1,ny
         dau(j) = (f1*dauin(j,m1) + f2*dauin(j,m2))
         dbu(j) = (f1*dbuin(j,m1) + f2*dbuin(j,m2))
         dav(j) = (f1*davin(j,m1) + f2*davin(j,m2))
         dbv(j) = (f1*dbvin(j,m1) + f2*dbvin(j,m2))
         dat(j) = (f1*datin(j,m1) + f2*datin(j,m2))
         dbt(j) = (f1*dbtin(j,m1) + f2*dbtin(j,m2))
         daz(j) = (f1*dazin(j,m1) + f2*dazin(j,m2)) * g
         dbz(j) = (f1*dbzin(j,m1) + f2*dbzin(j,m2)) * g

         sau(j) = (f1*sauin(j,m1) + f2*sauin(j,m2))
         sbu(j) = (f1*sbuin(j,m1) + f2*sbuin(j,m2))
         sav(j) = (f1*savin(j,m1) + f2*savin(j,m2))
         sbv(j) = (f1*sbvin(j,m1) + f2*sbvin(j,m2))
         sat(j) = (f1*satin(j,m1) + f2*satin(j,m2))
         sbt(j) = (f1*sbtin(j,m1) + f2*sbtin(j,m2))
         saz(j) = (f1*sazin(j,m1) + f2*sazin(j,m2)) * g
         sbz(j) = (f1*sbzin(j,m1) + f2*sbzin(j,m2)) * g
      end do

      end



