! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program advance_time


! interface identical to da_advance_cymdh, except for reading the arg line
! from standard input, to be more portable since iargc() is nonstandard across
! different fortran implementations.
!
! i/o sections of file lightly modified from da_advance_cymdh
! time computations all call DART time manager
!
!   - has accuracy down to second,
!   - can use day/hour/minute/second (with/without +/- sign) to advance time,
!   - can digest various input date format if it still has the right order (ie. cc yy mm dd hh nn ss)
!   - can digest flexible time increment 
!   - can output in wrf date format (ccyy-mm-dd_hh:nn:ss)
!   - can specify output date format
!   - can output Julian day
!   - can output Gregorian days and seconds (since year 1601)
!
! e.g:
!  echo 20070730      12         | advance_time    # advance 12 h 
!  echo 2007073012   -1d2h30m30s | advance_time    # back 1 day 2 hours 30 minutes and 30 seconds
!  echo 2007073012    1s-3h30m   | advance_time    # back 3 hours 30 minutes less 1 second
!  echo 200707301200  2d1s -w    | advance_time    # advance 2 days and 1 second, output in wrf date format
!  echo 2007-07-30_12:00:00 2d1s -w  | advance_time  # same as previous example
!  echo 200707301200  2d1s -f ccyy-mm-dd_hh:nn:ss | advance_time # same as previous example
!  echo 2007073006    120 -j     | advance_time    # advance 120 h, and print year and Julian day
!  echo 2007073006    120 -J     | advance_time    # advance 120 h, print year, Julian day, hour, minute and second
!  echo 2007073006    0 -g       | advance_time    # print Gregorian day and second (since year 1601)
!

use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                             increment_time, decrement_time, get_time, &
                             set_date, get_date, julian_day
use    utilities_mod, only : initialize_utilities
use   parse_args_mod, only : get_args_from_string

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

   integer :: ccyy, mm, dd, hh, nn, ss, dday, dh, dn, ds, gday, gsec
   integer :: nargum, i
   character(len=80), dimension(10) :: argum
   character(len=14) :: ccyymmddhhnnss
   character(len=80) :: out_date_format, dtime
   character(len=256) :: in_string
   integer :: datelen
   integer, parameter :: stdout=6
   type(time_type) :: base_time


   ! Initialize modules used that require it, and be silent about it
   call initialize_utilities('advance_time', output_flag = .false.)

   !call register_module(source,revision,revdate)


   call set_calendar_type(GREGORIAN)

   ! this routine reads a line from standard input and parses it up 
   ! into blank-separated words.
   read(*, '(A)') in_string
   call get_args_from_string(in_string, nargum, argum)

   if ( nargum < 2 ) then
      write(unit=stdout, fmt='(a)') &
         'Usage:   echo ccyymmddhh[nnss] [+|-]dt[d|h|m|s] [-w|-W|-wrf|-WRF] [-f|-F date_format] [-j|-J] [-g|-G] | advance_time'
      write(unit=stdout, fmt='(a)') &
         'Option:  -w|-W|-wrf|-WRF  output in wrf date format as ccyy-mm-dd_hh:nn:ss'
      write(unit=stdout, fmt='(a)') &
         '         -f|-F  specify output date format, such as ccyy-mm-dd_hh:nn:ss, or ''ccyy/mm/dd  hh:nn:ss'''
      write(unit=stdout, fmt='(a)') &
         '         -j|-J  print Julian day'
      write(unit=stdout, fmt='(a)') &
         '         -g|-G  print Gregorian days and seconds (since year 1601)'
      write(unit=stdout, fmt='(a)') &
         'Example: echo 20070730      12         | advance_time   # advance 12 h'
      write(unit=stdout, fmt='(a)') &
         '         echo 2007073012   -1d2h30m30s | advance_time   # back 1 day 2 hours 30 min and 30 sec'
      write(unit=stdout, fmt='(a)') &
         '         echo 2007073012    1s-3h30m   | advance_time   # back 3 hours 30 minutes less 1 second'
      write(unit=stdout, fmt='(a)') &
         '         echo 200707301200  1d1s -w    | advance_time   # advance 1 day 1 sec, output in wrf date format'
      write(unit=stdout, fmt='(a)') &
         '         echo 2007-07-30_12:00:00 2d1s -w | advance_time              # same as previous example'
      write(unit=stdout, fmt='(a)') &
         '         echo 200707301200  2d1s -f ccyy-mm-dd_hh:nn:ss | advance_time # same as previous' 
      write(unit=stdout, fmt='(a)') &
         '         echo 2007073006    120 -j     | advance_time    # advance 120 h, and print year and Julian day'
      write(unit=stdout, fmt='(a)') &
         '         echo 2007073006    120 -J     | advance_time    # advance 120 h, print year, Julian day, hour, minute and second'
      write(unit=stdout, fmt='(a)') &
         '         echo 2007073006    0 -g       | advance_time    # print Gregorian day and second (since year 1601)'
      write(unit=stdout, fmt='(a)') ''
      stop 'try again.'
   end if

   ccyymmddhhnnss = parsedate(argum(1))
   datelen = len_trim(ccyymmddhhnnss)

   if (datelen == 8) then
      read(ccyymmddhhnnss(1:10), fmt='(i4, 2i2)')  ccyy, mm, dd
      hh = 0
      nn = 0
      ss = 0
   else if (datelen == 10) then
      read(ccyymmddhhnnss(1:10), fmt='(i4, 3i2)')  ccyy, mm, dd, hh
      nn = 0
      ss = 0
   else if (datelen == 12) then
      read(ccyymmddhhnnss(1:12), fmt='(i4, 4i2)')  ccyy, mm, dd, hh, nn
      ss = 0
   else if (datelen == 14) then
      read(ccyymmddhhnnss(1:14), fmt='(i4, 5i2)')  ccyy, mm, dd, hh, nn, ss
   else
      stop 'wrong input date'
   endif

   base_time = set_date(ccyy, mm, dd, hh, nn, ss)


   dtime = trim(argum(2))
   call parsedt(dtime,dday,dh,dn,ds)


!print*, 'delta t: ', dday, dh, dn, ds

   ! each part can be positive or negative, or 0. 
   if (dday > 0) then
      base_time = increment_time(base_time, 0, dday)
   else if (dday < 0) then
      base_time = decrement_time(base_time, 0, -dday)
   endif
      
   if (dh > 0) then
      base_time = increment_time(base_time, dh*3600)
   else if (dh < 0) then
      base_time = decrement_time(base_time, -dh*3600)
   endif
      
   if (dn > 0) then
      base_time = increment_time(base_time, dn*60)
   else if (dn < 0) then
      base_time = decrement_time(base_time, -dn*60)
   endif
      
   if (ds > 0) then
      base_time = increment_time(base_time, ds)
   else if (ds < 0) then
      base_time = decrement_time(base_time, -ds)
   endif
      

   call get_date(base_time, ccyy, mm, dd, hh, nn, ss)


   write(ccyymmddhhnnss(1:14), fmt='(i4, 5i2.2)')  ccyy, mm, dd, hh, nn, ss
   if ( nargum == 2 ) then
      if (datelen<14) then
         if(nn /= 0) datelen=12
         if(ss /= 0) datelen=14
      endif
      write(unit=stdout, fmt='(a)') ccyymmddhhnnss(1:datelen)
   else if ( nargum > 2 ) then
      i = 3
      do while (i <= nargum)
        select case ( trim(argum(i)) )
           case ('-w', '-W', '-wrf','-WRF')
              out_date_format = 'ccyy-mm-dd_hh:nn:ss'
              write(unit=stdout, fmt='(a)') trim(formatdate(ccyymmddhhnnss, out_date_format))
              i = i+1
           case ('-f', '-F')
              out_date_format = trim(argum(i+1))
              write(unit=stdout, fmt='(a)') trim(formatdate(ccyymmddhhnnss, out_date_format))
              i = i+2
           case ('-j')
              write(unit=stdout, fmt='(I4,I4)') ccyy, julian_day(ccyy,mm,dd)
              i = i+1
           case ('-J')
              write(unit=stdout, fmt='(I4,I4,I3,I3,I3)') ccyy, julian_day(ccyy,mm,dd),hh,nn,ss
              i = i+1
           case ('-g','-G')
              call get_time(base_time, gsec, gday)
              write(unit=stdout, fmt='(I8,I8)') gday, gsec
              i = i+1
           case default
              i = i+1
        end select
      end do
   end if

contains


function parsedate(datein)
   character(len=80) :: datein
   character(len=14) :: parsedate
   character(len=1 ) :: ch
   integer :: n, i
   parsedate = '00000000000000'
   i=0
   do n = 1, len_trim(datein)
      ch = datein(n:n)
      if (ch >= '0' .and. ch <= '9') then
         i=i+1
         parsedate(i:i)=ch
      end if
   end do
   if (parsedate(11:14) == '0000') then
      parsedate(11:14) = ''
   else if(parsedate(13:14) == '00') then
      parsedate(13:14) = ''
   end if
   return 
end function parsedate

subroutine parsedt(dt,dday,dh,dn,ds)
   character(len=80) :: dt
   integer :: dday, dh, dn, ds
   character(len=1 ) :: ch
   integer :: n,i,d,s,nounit
   ! initialize time and sign
   nounit=1
   dday=0
   dh=0
   dn=0
   ds=0
   d=0
   s=1
   do n = 1, len_trim(dt)
      ch = dt(n:n)
      select case (ch)
         case ('0':'9')
           read(ch,fmt='(i1)') i
           d=d*10+i
         case ('-')
           s=-1
         case ('+')
           s=1
         case ('d')
           nounit=0
           dday=dday+d*s
           d=0
         case ('h')
           nounit=0
           dh=dh+d*s
           d=0
         case ('n','m')
           nounit=0
           dn=dn+d*s
           d=0
         case ('s')
           nounit=0
           ds=ds+d*s
           d=0
         case default
      end select
   end do
   if (nounit==1) dh=d*s
end subroutine parsedt

function formatdate(datein,dateform)
   character(len=14) :: datein
   character(len=80) :: dateform
   character(len=80) :: formatdate
   integer :: ic,iy,im,id,ih,in,is
   ic=index(dateform,'cc')
   iy=index(dateform,'yy')
   im=index(dateform,'mm')
   id=index(dateform,'dd')
   ih=index(dateform,'hh')
   in=index(dateform,'nn')
   is=index(dateform,'ss')
   formatdate=trim(dateform)
   if (ic /= 0) formatdate(ic:ic+1) = datein(1:2)
   if (iy /= 0) formatdate(iy:iy+1) = datein(3:4)
   if (im /= 0) formatdate(im:im+1) = datein(5:6)
   if (id /= 0) formatdate(id:id+1) = datein(7:8)
   if (ih /= 0) formatdate(ih:ih+1) = datein(9:10)
   if (in /= 0) formatdate(in:in+1) = datein(11:12)
   if (is /= 0) formatdate(is:is+1) = datein(13:14)
   return
end function formatdate


end program advance_time

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
