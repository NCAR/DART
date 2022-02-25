! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Compute with time quantities
!>
!> The advance_time program computes the resulting time when either
!> adding or subtracting time intervals.  The increments can be
!> expressed in days, hours, minutes or seconds.  The output can be
!> formatted as native WRF, CESM, Julian or Gregorian format.
!>
!> Reads input from standard input to be more portable, since older
!> versions of iargc() weren't standardized.
!> 
!> Based on the WRF da_advance_cymdh utility.
!>
!> All time computations call DART time manager.
!>
!>   - has accuracy down to second,
!>   - can use day/hour/minute/second (with/without +/- sign) to advance time,
!>   - can digest various input date format if it still has the right order (ie. cc yy mm dd hh nn ss)
!>   - can digest flexible time increment 
!>   - can output in wrf date format (ccyy-mm-dd_hh:nn:ss)
!>   - can specify output date format
!>   - can output Julian day
!>   - can output Gregorian days and seconds (since year 1601)
!>   - can output in CESM time format (ccyy-mm-dd-fffff where fffff is seconds of day)
!>
!> Examples
!>
!> - advance 12 h
!>
!>    echo 20070730      12         | advance_time
!>
!> - go back 1 day 2 hours 30 minutes and 30 seconds
!>
!>    echo 2007073012   -1d2h30m30s | advance_time
!>
!> - go back 3 hours 30 minutes less 1 second
!>
!>    echo 2007073012    1s-3h30m   | advance_time
!>
!> - advance 2 days and 1 second, output in wrf date format (three ways)
!>
!>    echo 200707301200        2d1s -w | advance_time
!>
!>    echo 2007-07-30_12:00:00 2d1s -w | advance_time
!>
!>    echo 200707301200        2d1s -f ccyy-mm-dd_hh:nn:ss | advance_time
!>
!> - advance 120 h, and print year and Julian day
!>
!>    echo 2007073006    120 -j     | advance_time
!>
!> - advance 120 h, print year, Julian day, hour, minute and second
!>
!>    echo 2007073006    120 -J     | advance_time
!>
!> - print Gregorian day and second (since year 1601)
!>
!>    echo 2007073006    0 -g       | advance_time
!>
!> - print CESM format time (ccyy-mm-dd-fffff where fffff is sec of day)
!>
!>    echo 2007073006    0 -c       | advance_time
!>
!> @todo if run with no inputs ... it just hangs. Can we make it fail straight away?

program advance_time

use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                             increment_time, decrement_time, get_time, &
                             set_date, get_date, julian_day
use    utilities_mod, only : initialize_utilities, E_ERR, error_handler, &
                             finalize_utilities
use   parse_args_mod, only : get_args_from_string

implicit none

character(len=*), parameter :: source = 'advance_time.f90'

integer :: ccyy, mm, dd, hh, nn, ss, dday, dh, dn, ds, gday, gsec
integer :: nargum, i
character(len=80), dimension(10) :: argum
character(len=14) :: ccyymmddhhnnss
character(len=80) :: out_date_format, dtime
character(len=256) :: in_string
integer :: datelen
integer, parameter :: stdout=6
type(time_type) :: base_time

character(len=512) :: string1, string2

! Initialize modules used that require it, and be silent about it
call initialize_utilities('advance_time', output_flag = .false.)

call set_calendar_type(GREGORIAN)

! this routine reads a line from standard input and parses it up 
! into blank-separated words.
read(*, '(A)') in_string
call get_args_from_string(in_string, nargum, argum)

if ( nargum < 2 ) then
   write(unit=stdout, fmt='(a)') &
   'Usage:   echo ccyymmddhh[nnss] [+|-]dt[d|h|m|s] [-w|-W|-wrf|-WRF] [-f|-F date_format] [-j|-J] [-g|-G] [-c|-C] | advance_time'
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
   write(unit=stdout, fmt='(a)') &
      '         echo 2007073006    0 -c       | advance_time    # print CESM format time (ccyy-mm-dd-fffff where fffff is sec of day)'
   write(unit=stdout, fmt='(a)') ''

   call error_handler(E_ERR,'advance_time','Invalid Usage', source)
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
elseif (datelen == 13) then
   read(ccyymmddhhnnss(1:13), fmt='(i4, 2i2, i5)')  ccyy, mm, dd, ss
   if (ss >= 86400) then
      write(string1,*)'seconds-of-day is ',ss,' as parsed from ',trim(ccyymmddhhnnss)
      write(string2,*)'seconds-of-day  must be less than 86400'
      call error_handler(E_ERR,'advance_time',string1, source, text2=string2)
   endif
   hh = ss / 3600
   ss = ss - hh * 3600
   nn = ss / 60
   ss = ss - nn * 60
else
   write(string1,*)'unsupported format for ',trim(argum(1))
   call error_handler(E_ERR,'advance_time',string1, source)
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
   if (datelen == 13) datelen=10
   if (datelen < 14) then
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
        case ('-c','-C')
           out_date_format = 'ccyy-mm-dd-fffff'
           write(unit=stdout, fmt='(a)') trim(formatCESMdate(ccyy,mm,dd,hh,nn,ss))
           i = i+1
        case default
           i = i+1
     end select
   end do
end if

call finalize_utilities()

contains

!-----------------------------------------------------------------------
!> removes non-numeric characters from a date string
!>
!> @param[in] datein character string containing the date. May include
!>               dashes, colons, etc.
!> @return character string with only numeric characters

function parsedate(datein)
character(len=*), intent(in) :: datein

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

if (i == 13) then
   parsedate(14:14) = ''
   return  ! CESM format
else if (parsedate(11:14) == '0000') then
   parsedate(11:14) = ''
else if(parsedate(13:14) == '00') then
   parsedate(13:14) = ''
end if

end function parsedate

!-----------------------------------------------------------------------
!> extracts the day,hour,minutes and seconds from the second input argument
!>
!> @param[in] dt character string with the temporal offset (the second input argument).
!> @param[out] dday the day
!> @param[out] dh the hour
!> @param[out] dn the minute
!> @param[out] ds the second

subroutine parsedt(dt,dday,dh,dn,ds)
character(len=*), intent(in)  :: dt
integer,          intent(out) :: dday
integer,          intent(out) :: dh
integer,          intent(out) :: dn
integer,          intent(out) :: ds

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
        continue
   end select
end do

if (nounit==1) dh=d*s

end subroutine parsedt


!-----------------------------------------------------------------------
!>
!> @param[in] datein date in the known character string of length 14
!> @param[in] dateform character string containing desired date format
!> @return character string containing desired date format

function formatdate(datein,dateform)
character(len=*), intent(in) :: datein
character(len=*), intent(in) :: dateform

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

end function formatdate

!-----------------------------------------------------------------------
!>
!> @param[in] ccyy the century and year
!> @param[in] mm the month
!> @param[in] dd the day
!> @param[in] hh the hour
!> @param[in] nn the minute
!> @param[in] ss the second
!> @return character string of date in CESM format YYYY-MM-DD-SSSSS

function formatCESMdate(ccyy,mm,dd,hh,nn,ss)
integer, intent(in) :: ccyy
integer, intent(in) :: mm
integer, intent(in) :: dd
integer, intent(in) :: hh
integer, intent(in) :: nn
integer, intent(in) :: ss
character(len=80) :: formatCESMdate

integer :: fffff

fffff = hh*3600 + nn*60 + ss
write(formatCESMdate, '(i4.4,1a,2(i2.2,1a),i5.5)') ccyy, '-', mm, '-', dd, '-', fffff

end function formatCESMdate


end program advance_time

