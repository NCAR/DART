
      module MO_CALENDAR

      implicit none

      save

      character(len=16) :: &
        type = 'gregorian       '  ! calendar type.  Currently '365' or 'gregorian'

      CONTAINS

      subroutine INICALENDAR( xtype )
!-----------------------------------------------------------------------
! 	... Initialize calendar module common block.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy argument
!-----------------------------------------------------------------------
      character(len=16), intent(in) :: &
        xtype                     ! calendar type.  Currently '365' or 'gregorian'

      type = xtype

      end subroutine INICALENDAR

      real function DIFFDAT( dat1, sec1, dat2, sec2 )
!-----------------------------------------------------------------------
! 	... Compute the difference: (dat2,sec2) - (dat1,sec1)  in days.
!           Return value:
!           (dat2,sec2) - (dat1,sec1)  in days.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Input arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
        dat1, &   ! date in yyyymmdd format
        sec1, &   ! seconds relative to dat1
        dat2, &   ! date in yyyymmdd format
        sec2      ! seconds relative to dat2

      if( type(1:3)  ==  '365' ) then
         DIFFDAT = DIFFDAT365( dat1, sec1, dat2, sec2 )
      else if ( type(1:9)  ==  'gregorian' ) then
         DIFFDAT = DIFFDATGRG( dat1, sec1, dat2, sec2 )
      end if

      end function DIFFDAT

      real function DIFFDAT365( dat1, sec1, dat2, sec2 )
!-----------------------------------------------------------------------
! 	... Compute the difference: (dat2,sec2) - (dat1,sec1)  in days.
! 	    N.B. Assume 1 year = 365 days.
!           Return value:
!           (dat2,sec2) - (dat1,sec1)  in days.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
        dat1, &   ! date in yyyymmdd format
        sec1, &   ! seconds relative to dat1
        dat2, &   ! date in yyyymmdd format
        sec2      ! seconds relative to dat2


!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer :: lastd, lasts, firstd, firsts, &
                 yr2, mo2, dy2, doy2, yr1, mo1, dy1, doy1, ndays
      real :: sign, days

!-----------------------------------------------------------------------
!     	... Dates equal
!-----------------------------------------------------------------------
      if( dat1 == dat2 .and. sec1 == sec2 ) then
         DIFFDAT365 = 0.
         return
      end if

!-----------------------------------------------------------------------
!     	... Which date is later?
!-----------------------------------------------------------------------
      if( dat2 > dat1 ) then
         sign = 1.
         lastd  = dat2
         lasts  = sec2
         firstd = dat1
         firsts = sec1
      else if( dat2 < dat1 ) then
         sign   = -1.
         lastd  = dat1
         lasts  = sec1
         firstd = dat2
         firsts = sec2
      else
         if( sec2 > sec1 ) then
            sign   = 1.
            lastd  = dat2
            lasts  = sec2
            firstd = dat1
            firsts = sec1
         else
            sign   = -1.
            lastd  = dat1
            lasts  = sec1
            firstd = dat2
            firsts = sec2
         end if
      end if

!-----------------------------------------------------------------------
!     	... Compute number of days between lastd and firstd
!-----------------------------------------------------------------------
      yr2  = lastd / 10000
      mo2  = MOD( lastd, 10000 ) / 100
      dy2  = MOD( lastd, 100 )
      doy2 = DOY( mo2, dy2 )

      yr1  = firstd / 10000
      mo1  = MOD( firstd, 10000 ) / 100
      dy1  = MOD( firstd, 100 )
      doy1 = DOY( mo1, dy1 )

      ndays = 365*(yr2 - yr1) + doy2 - doy1

!-----------------------------------------------------------------------
!     	... Adjust for remaining seconds
!-----------------------------------------------------------------------
      days = REAL( ndays ) + REAL( lasts - firsts )/86400.

!-----------------------------------------------------------------------
!     	... Adjust sign
!-----------------------------------------------------------------------
      DIFFDAT365 = sign * days

      end function DIFFDAT365

      real function DIFFDATGRG( dat1, sec1, dat2, sec2 )
!-----------------------------------------------------------------------
! 	... Compute the difference: (dat2,sec2) - (dat1,sec1)  in days.
!           N.B. Assume Gregorian calendar.
!           Return value:
!           (dat2,sec2) - (dat1,sec1)  in days.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
        dat1, &   ! date in yyyymmdd format
        sec1, &   ! seconds relative to dat1
        dat2, &   ! date in yyyymmdd format
        sec2      ! seconds relative to dat2


!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer :: lastd, lasts, firstd, firsts, ndays
      real    :: sign, days

!-----------------------------------------------------------------------
!     	... Dates equal?
!-----------------------------------------------------------------------
      if( dat1 == dat2 .and. sec1 == sec2 ) then
         DIFFDATGRG = 0.
         return
      end if

!-----------------------------------------------------------------------
!     	... Which date is later?
!-----------------------------------------------------------------------
      if( dat2 > dat1 ) then
         sign   = 1.
         lastd  = dat2
         lasts  = sec2
         firstd = dat1
         firsts = sec1
      else if( dat2 < dat1 ) then
         sign   = -1.
         lastd  = dat1
         lasts  = sec1
         firstd = dat2
         firsts = sec2
      else
         if( sec2 > sec1 ) then
            sign   = 1.
            lastd  = dat2
            lasts  = sec2
            firstd = dat1
            firsts = sec1
         else
            sign   = -1.
            lastd  = dat1
            lasts  = sec1
            firstd = dat2
            firsts = sec2
         end if
      end if

!-----------------------------------------------------------------------
!     	... Compute number of days between lastd and firstd
!-----------------------------------------------------------------------
      ndays = GREG2JDAY( lastd ) - GREG2JDAY( firstd )

!-----------------------------------------------------------------------
!     	... Adjust for remaining seconds
!-----------------------------------------------------------------------
      days = REAL( ndays ) + REAL( lasts - firsts )/86400.

!-----------------------------------------------------------------------
!     	... Adjust sign
!-----------------------------------------------------------------------
      DIFFDATGRG = sign * days

      end function DIFFDATGRG

      subroutine ADDSEC2DAT( increment, dat1, sec1 )
!-----------------------------------------------------------------------
! 	... Add seconds to date and compute new date.
!           (dat1,sec1) + sec = (dat1,sec1).
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
        increment    ! seconds to add to (dat1,sec1) (can be negative)
      integer, intent(inout) :: &
        dat1, &      ! date in yyyymmdd format
        sec1         ! seconds relative to dat1

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer :: sec2

      if( increment /= 0 ) then
         sec2 = sec1 + increment
         sec1 = MOD( sec2, 86400 )
         dat1 = NEWDATE( dat1, sec2/86400 )
         if( sec1 < 0 ) then
            sec1 = sec1 + 86400
            dat1 = NEWDATE( dat1, -1 )
         end if
      end if

      end subroutine ADDSEC2DAT

      integer function DOY( mon, cday )
!-----------------------------------------------------------------------
! 	... Compute day of year ignoring leap years.
!           Returns values in the range [1,365].
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: mon, cday

      integer, save :: jdcon(12) = &
            (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/)

      doy = jdcon(mon) + cday

      end function DOY

      integer function NEWDATE( date, iday )
!-----------------------------------------------------------------------
! 	... Find the date a specified number of days (possibly negative)
!           from the given date.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
         date, &               ! starting date in yyyymmdd format (can be negative)
         iday                  ! number of days to increment date (can be negative)

      if( type(1:3) == '365' ) then
         NEWDATE = NEWDATE365( date, iday )
      else if ( type(1:9)  ==  'gregorian' ) then
         NEWDATE = NEWDATEGRG( date, iday )
      end if

      end function NEWDATE

      integer function NEWDATE365( date, iday )
!-----------------------------------------------------------------------
! 	... Find the date a specified number of days (possibly negative)
!           from the given date.
!           N.B. Assume 1 year = 365 days.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
         date, &               ! starting date in yyyymmdd format (can be negative)
         iday                  ! number of days to increment date (can be negative)

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer :: &
        idate, &               ! intial value of date (may change sign)
        year, &                ! year of date
        month, &               ! month of date
        imonth, &              ! save the initial value of month
        day, &                 ! day of date
        ndays, &               ! running count of days to increment date by
        i, &                   ! index
        cmonth, &              ! current month as we increment date
        nd                     ! days to end of month
      integer, save :: mdays(24) = &
                (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                   31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

      if( iday == 0 ) then
         NEWDATE365 = date
         return
      end if

      idate = date
      year  = idate / 10000
      if( idate < 0 ) then
         idate = -idate
      end if
      month  = MOD( idate,10000 ) / 100
      imonth = month
      day    = MOD( idate,100 )
      
      if( iday  > 0 ) then
         ndays = iday
         year  = year + ndays/365
         ndays = MOD( ndays,365 )
         do i = imonth,imonth+12
            cmonth = i  
            if( i > 12 ) then
	       cmonth = i - 12
	    end if
            nd = ND2ENDM( cmonth, day )
            if( ndays > nd ) then
               month = month + 1
               if( month > 12 ) then
                  month = 1
                  year = year + 1
               end if
               day = 1
               ndays = ndays - nd - 1
               if( ndays == 0 ) then
	          exit
	       end if
            else
               day = day + ndays
               exit
            end if
         end do
      else if( iday < 0 ) then
         ndays = -iday
         year  = year - ndays/365
         ndays = MOD( ndays,365 )
         imonth = month
         do i = imonth+12,imonth, -1
            if( ndays >= day ) then
               month = month - 1
               if( month < 1 ) then
                  month = 12
                  year = year - 1
               end if
               ndays = ndays - day
               day = mdays( month )
               if( ndays == 0 ) then
	          exit
	       end if
            else
               day = day - ndays
               exit
            end if
         end do
      end if

      if( year >= 0 ) then
         NEWDATE365 = year*10000 + month*100 + day
      else
         NEWDATE365 = -year*10000 + month*100 + day
         NEWDATE365 = -NEWDATE365
      end if

      end function NEWDATE365

      integer function NEWDATEGRG( date, iday )
!-----------------------------------------------------------------------
! 	... Find the date a specified number of days (possibly negative)
!           from the given date.
!           N.B. Assume Gregorian calendar
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
         date, &                         ! starting date in yyyymmdd format (can be negative)
         iday                            ! number of days to increment date (can be negative)

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer :: &
        curjday, newjday

      if( iday == 0 ) then
         NEWDATEGRG = date
      else
         NEWDATEGRG = JDAY2GREG( GREG2JDAY( date ) + iday )
      end if

      end function NEWDATEGRG

      real function CALDAYR( idate, isec )
!-----------------------------------------------------------------------
! 	... Calendar day with fractional part.  Returns values in the
!           range [1., 366.)
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: idate, isec

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer :: mon, day

      mon = MOD( idate,10000 ) / 100
      day = MOD( idate,100 )

      CALDAYR = DOY( mon, day ) + REAL( isec )/86400.

      end function CALDAYR

!-----------------------------------------------------------------------
! Currently the routines greg2jday and jday2greg are working in sense of
! j = greg2jday( jday2greg( j ) ) down to jday 1684595 which corresponds
! approx. to the date -1000228.
!-----------------------------------------------------------------------

      integer function GREG2JDAY( date )
!-----------------------------------------------------------------------
! 	... Return Julian day number given Gregorian date.
!
! Algorithm from Hatcher,D.A., Simple Formulae for Julian Day Numbers
! and Calendar Dates, Q.Jl.R.astr.Soc. (1984) v25, pp 53-55.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: date

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer :: yy, mm, dd
      integer :: ap, mp
      integer :: y, d, n, g

!-----------------------------------------------------------------------
!     	... Extract year, month, and day from date
!-----------------------------------------------------------------------
      yy = date / 10000
      mm = MOD( ABS(date),10000 ) / 100
      dd = MOD( ABS(date),100 )

!-----------------------------------------------------------------------
!     	... Modify year and month numbers
!-----------------------------------------------------------------------
      ap = yy - (12 - mm)/10
      mp = MOD( mm-3,12 )
      if( mp < 0 ) then
         mp = mp + 12
      end if

!-----------------------------------------------------------------------
!     	... Julian day
!-----------------------------------------------------------------------
      y = INT( 365.25*( ap + 4712 ) )
      d = INT( 30.6*mp + .5 )
      n = y + d + dd  + 59
      g = INT( .75*INT( ap/100 + 49 ) ) - 38
      GREG2JDAY = n - g

      end function GREG2JDAY

      integer function JDAY2GREG( day )
!-----------------------------------------------------------------------
! 	... Return Gregorian date given Julian day number.
!
! Algorithm from Hatcher,D.A., Simple Formulae for Julian Day Numbers
! and Calendar Dates, Q.Jl.R.astr.Soc. (1984) v25, pp 53-55.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: day

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer :: g, n
      integer :: yy, dp, mm, dd

      g = INT( .75*INT( (day - 4479.5)/36524.25 ) + .5 ) - 37
      n = day + g

      yy = INT( n/365.25 ) - 4712
      dp = INT( MOD( n-59.25, 365.25 ) )
      mm = MOD( INT( (dp+.5)/30.6 )+2, 12 ) + 1
      dd = INT( MOD( dp+.5, 30.6 ) ) + 1

      JDAY2GREG = SIGN( ABS(yy)*10000 + mm*100 +dd,yy )

      end function JDAY2GREG

      integer function ND2ENDM( mon, day )
!-----------------------------------------------------------------------
! 	... Returns the number of days to the end of the month.
!           This number added to the input arguement day gives the
!           last day of month mon.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: mon, day

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer, save :: mdays(12) = &
         (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

      if( mon < 1  .or. mon > 12 ) then
         write(*,*) 'ND2EDMm: illegal input... mon= ', mon
         stop
      end if

      ND2ENDM = mdays(mon) - day

      if( day < 0 ) then
         write(*,*) 'ND2ENDM: illegal input... mon= ',mon,' day= ',day
         stop
      end if

      end function ND2ENDM

      end module MO_CALENDAR
