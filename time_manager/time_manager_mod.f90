module time_manager_mod

! The use of fms utilities below have all been temporarily commented out.
! Should be re-examined.
!!!use fms_mod, only: error_mesg, FATAL, write_version_number, stdout


use types_mod


implicit none
private

!====================================================================
! The time_manager provides a single defined type, time_type, which is 
! used to store time and date quantities. A time_type is a positive 
! definite quantity that represents an interval of time. It can be most 
! easily thought of as representing the number of seconds in some time 
! interval. A time interval can be mapped to a date under a given calendar 
! definition by using it to represent the time that has passed since some 
! base date. A number of interfaces are provided to operate on time_type 
! variables and their associated calendars. Time intervals can be as large 
! as n days where n is the largest number represented by the default integer 
! type on a compiler. This is typically considerably greater than 10 million 
! years which is likely to be adequate for most applications. The 
! description of the interfaces is separated into two sections. The first 
! deals with operations on time intervals while the second deals with 
! operations that convert time intervals to dates for a given calendar.
!
!====================================================================

! Module defines a single type
public time_type

! Operators defined on time_type
public operator(+),  operator(-),   operator(*),   operator(/),  &
       operator(>),  operator(>=),  operator(==),  operator(/=), &
       operator(<),  operator(<=),  operator(//)

! Subroutines and functions operating on time_type
public set_time, increment_time, decrement_time, get_time, interval_alarm
public repeat_alarm, write_time, read_time, nc_append_time

! List of available calendar types
public    THIRTY_DAY_MONTHS,    JULIAN,    GREGORIAN,  NO_LEAP,   NO_CALENDAR

! Subroutines and functions involving relations between time and calendar
public set_calendar_type, get_calendar_type
public set_date,       set_date_gregorian,         set_date_julian, &
                       set_date_thirty,            set_date_no_leap
public get_date,       get_date_gregorian,         get_date_julian, &
                       get_date_thirty,            get_date_no_leap
public increment_date, increment_gregorian,        increment_julian, &
                       increment_thirty,           increment_no_leap
public decrement_date, decrement_gregorian,        decrement_julian, &
                       decrement_thirty,           decrement_no_leap
public days_in_month,  days_in_month_gregorian,    days_in_month_julian, &
                       days_in_month_no_leap,      days_in_month_thirty
public leap_year,      leap_year_gregorian,        leap_year_julian, &
                       leap_year_no_leap,          leap_year_thirty
public length_of_year, length_of_year_thirty,      length_of_year_julian, &
                       length_of_year_gregorian,   length_of_year_no_leap
public days_in_year,   days_in_year_thirty,        days_in_year_julian, &
                       days_in_year_gregorian,     days_in_year_no_leap
public month_name

! Subroutines for printing version number and time type
public :: time_manager_init, print_time, print_date

!====================================================================

! Global data to define calendar type
integer, parameter :: THIRTY_DAY_MONTHS = 1,      JULIAN = 2, &
                      GREGORIAN = 3,              NO_LEAP = 4, &
                      NO_CALENDAR = 0
integer, private :: calendar_type = NO_CALENDAR, max_type = 4

! Define number of days per month
integer, private :: days_per_month(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)

! time_type is implemented as seconds and days to allow for larger intervals
type time_type
   private
   integer:: seconds
   integer:: days
end type time_type

!======================================================================

interface operator (+);   module procedure time_plus;        end interface
interface operator (-);   module procedure time_minus;       end interface
interface operator (*);   module procedure time_scalar_mult 
                          module procedure scalar_time_mult; end interface
interface operator (/);   module procedure time_scalar_divide
                          module procedure time_divide;      end interface
interface operator (>);   module procedure time_gt;          end interface
interface operator (>=);  module procedure time_ge;          end interface
interface operator (<);   module procedure time_lt;          end interface
interface operator (<=);  module procedure time_le;          end interface
interface operator (==);  module procedure time_eq;          end interface
interface operator (/=);  module procedure time_ne;          end interface
interface operator (//);  module procedure time_real_divide; end interface

!======================================================================

character(len=128) :: version='$Id$'
character(len=128) :: tagname='$Id$'
logical :: do_init = .true.

!======================================================================

contains

! First define all operations on time intervals independent of calendar



function set_time(seconds, days)
!---------------------------------------------------------------------------
!
! Returns a time interval corresponding to this number of days and seconds.
! The arguments must not be negative but are otherwise unrestricted.

implicit none

integer, intent(in)           :: seconds
integer, intent(in), optional :: days
type(time_type)               :: set_time

integer :: days_in

days_in = 0;  if (present(days)) days_in = days

! Negative time offset is illegal

if(seconds < 0 .or. days_in < 0) call error_handler('Negative input in set_time')

! Make sure seconds greater than a day are fixed up

set_time%seconds = seconds - seconds / (60*60*24) * (60*60*24)

! Check for overflow on days before doing operation
if(seconds / (60*60*24)  >= huge(days_in) - days_in) &
   call error_handler('Integer overflow in days in set_time')
set_time%days = days_in + seconds / (60*60*24)

end function set_time



subroutine get_time(time, seconds, days)
!---------------------------------------------------------------------------
!
! Returns days and seconds ( < 86400 ) corresponding to a time.

implicit none

type(time_type), intent(in)            :: time
integer,         intent(out)           :: seconds
integer,         intent(out), optional :: days

seconds = time%seconds
if (present(days)) then
  days = time%days
else
  if (time%days > (huge(seconds) - seconds)/(60*60*24)) &
  call error_handler('Integer overflow in seconds in get_time use days')
  seconds = seconds + time%days * (60*60*24)
endif

end subroutine get_time



function increment_time(time, seconds, days)
!-------------------------------------------------------------------------
!
! Increments a time by seconds and days; increments cannot be negative.

implicit none

type(time_type), intent(in)           :: time
integer,         intent(in)           :: seconds
integer,         intent(in), optional :: days
type(time_type)                       :: increment_time

integer :: days_in

days_in = 0;  if (present(days)) days_in = days

! Increment must be positive definite

if(seconds < 0 .or. days_in < 0) &
   call error_handler('Negative increment in increment_time')

! Watch for immediate overflow on days or seconds

if(days_in >= huge(days_in) - time%days) &
   call error_handler('Integer overflow in days in increment_time')

if(seconds >= huge(seconds) - time%seconds) &
   call error_handler('Integer overflow in seconds in increment_time')

increment_time = set_time(time%seconds + seconds, time%days + days_in)

end function increment_time



function decrement_time(time, seconds, days)
!--------------------------------------------------------------------------
!
! Decrements a time by seconds and days; decrements cannot be negative.

implicit none

type(time_type), intent(in)           :: time
integer,         intent(in)           :: seconds
integer,         intent(in), optional :: days
type(time_type)                       :: decrement_time

integer :: cseconds, cdays

cdays = 0;  if (present(days)) cdays = days

! Decrement must be positive definite

if(seconds < 0 .or. cdays < 0) &
   call error_handler('Negative decrement in decrement_time')

cseconds = time%seconds - seconds
cdays    = time%days - cdays

! Borrow if needed

if(cseconds < 0) then
   cdays = cdays - 1 + (cseconds + 1) / (60*60*24)
   cseconds = cseconds - (60*60*24) * (-1 + (cseconds + 1) / (60*60*24))
end if

! Check for illegal negative time

if(cdays < 0) call error_handler('Negative time results in decrement_time')

decrement_time%seconds = cseconds
decrement_time%days    = cdays

end function decrement_time



function time_gt(time1, time2)
!--------------------------------------------------------------------------
!
! Returns true if time1 > time2

implicit none

type(time_type), intent(in) :: time1, time2
logical                     :: time_gt

time_gt = (time1%days > time2%days)

if(time1%days == time2%days) time_gt = (time1%seconds > time2%seconds)

end function time_gt



function time_ge(time1, time2)
!--------------------------------------------------------------------------
!
! Returns true if time1 >= time2

implicit none

type(time_type), intent(in) :: time1, time2
logical                     :: time_ge

time_ge = (time_gt(time1, time2) .or. time_eq(time1, time2))

end function time_ge



function time_lt(time1, time2)
!--------------------------------------------------------------------------
!
! Returns true if time1 < time2

implicit none

type(time_type), intent(in) :: time1, time2
logical                     :: time_lt

time_lt = (time1%days < time2%days)
if(time1%days == time2%days) time_lt = (time1%seconds < time2%seconds)

end function time_lt



function time_le(time1, time2)
!--------------------------------------------------------------------------
!
! Returns true if time1 <= time2

implicit none

type(time_type), intent(in) :: time1, time2
logical                     :: time_le

time_le = (time_lt(time1, time2) .or. time_eq(time1, time2))

end function time_le



function time_eq(time1, time2)
!--------------------------------------------------------------------------
!
! Returns true if time1 == time2

implicit none

type(time_type), intent(in) :: time1, time2
logical                     :: time_eq

time_eq = (time1%seconds == time2%seconds .and. time1%days == time2%days)

end function time_eq



function time_ne(time1, time2)
!--------------------------------------------------------------------------
!
! Returns true if time1 /= time2

implicit none

type(time_type), intent(in) :: time1, time2
logical                     :: time_ne

time_ne = (.not. time_eq(time1, time2))

end function time_ne



function time_plus(time1, time2)
!-------------------------------------------------------------------------
!
! Returns sum of two time_types

implicit none

type(time_type), intent(in) :: time1, time2
type(time_type)             :: time_plus

time_plus = increment_time(time1, time2%seconds, time2%days)

end function time_plus



function time_minus(time1, time2)
!-------------------------------------------------------------------------
!
! Returns difference of two time_types. WARNING: a time type is positive 
! so by definition time1 - time2  is the same as time2 - time1.

implicit none

type(time_type), intent(in) :: time1, time2
type(time_type)             :: time_minus

if(time1 > time2) then
   time_minus = decrement_time(time1, time2%seconds, time2%days)
else 
   time_minus = decrement_time(time2, time1%seconds, time1%days)
endif

end function time_minus



function time_scalar_mult(time, n)
!--------------------------------------------------------------------------
!
! Returns time multiplied by integer factor n

implicit none

type(time_type), intent(in) :: time
integer,         intent(in) :: n
type(time_type)             :: time_scalar_mult

integer  :: days, seconds
real(r8) :: sec_prod 

! Multiplying here in a reasonable fashion to avoid overflow is tricky
! Could multiply by some large factor n, and seconds could be up to 86399
! Need to avoid overflowing integers and wrapping around to negatives

sec_prod = (time%seconds*1.0_r8) * (n*1.0_r8)

! If sec_prod is large compared to precision of double precision, things
! can go bad.  Need to warn and abort on this.

if(sec_prod /= 0.0_r8) then
   if(log10(sec_prod) > precision(sec_prod) - 3) call error_handler( &
      'Insufficient precision to handle scalar product in time_scalar_mult; contact developer')
end if

days    = sec_prod / (24.0_r8 * 60._r8 * 60.0_r8)
seconds = sec_prod - days * (24.0_r8 * 60.0_r8 * 60.0_r8)    ! QJLA nint or round or ...

time_scalar_mult = set_time(seconds, time%days * n + days)

end function time_scalar_mult



function scalar_time_mult(n, time)
!-------------------------------------------------------------------------
!
! Returns time multipled by integer factor n

implicit none

integer,         intent(in) :: n
type(time_type), intent(in) :: time
type(time_type)             :: scalar_time_mult

scalar_time_mult = time_scalar_mult(time, n)

end function scalar_time_mult



function time_divide(time1, time2)
!-------------------------------------------------------------------------
!
! Returns the largest integer, n, for which time1 >= time2 * n.

implicit none

type(time_type), intent(in) :: time1, time2
integer                     :: time_divide

real(r8) :: d1, d2

! Convert time intervals to floating point days; risky for general performance?

d1 = time1%days * (60.0_r8 * 60.0_r8 * 24.0_r8) + (time1%seconds*1.0_r8)
d2 = time2%days * (60.0_r8 * 60.0_r8 * 24.0_r8) + (time2%seconds*1.0_r8) 

! Get integer quotient of this, check carefully to avoid round-off problems.

time_divide = d1 / d2

! Verify time_divide*time2 is <= time1 and (time_divide + 1)*time2 is > time1

if(time_divide * time2 > time1 .or. (time_divide + 1) * time2 <= time1) &
   call error_handler('time_divide quotient error :: notify developer')

end function time_divide



function time_real_divide(time1, time2)
!-------------------------------------------------------------------------
!
! Returns the double precision quotient of two times

implicit none

type(time_type), intent(in) :: time1, time2
real(r8)                    :: time_real_divide

real(r8) :: d1, d2

! Convert time intervals to floating point days; risky for general performance?

d1 = time1%days * (60.0_r8 * 60.0_r8 * 24.0_r8) + (time1%seconds*1.0_r8)
d2 = time2%days * (60.0_r8 * 60.0_r8 * 24.0_r8) + (time2%seconds*1.0_r8) 

time_real_divide = d1 / d2

end function time_real_divide



function time_scalar_divide(time, n)
!-------------------------------------------------------------------------
!
! Returns the largest time, t, for which n * t <= time

implicit none

integer,         intent(in) :: n
type(time_type), intent(in) :: time
type(time_type)             :: time_scalar_divide

real(r8)        :: d, div
integer         :: days, seconds
type(time_type) :: prod1, prod2

! Convert time interval to floating point days; risky for general performance?

d   = time%days * (60.0_r8 * 60.0_r8 * 24.0_r8) + (time%seconds*1.0_r8)
div = d / (1.0_r8 * n)

days    = div / (60.0_r8 * 60.0_r8 * 24.0_r8)
seconds = div - days * (60.0_r8 * 60.0_r8 * 24.0_r8)
time_scalar_divide = set_time(seconds, days)

! Need to make sure that roundoff isn't killing this

prod1 = n * time_scalar_divide
prod2 = n * (increment_time(time_scalar_divide, 1, 0)) 

if(prod1 > time .or. prod2 <= time) &
   call error_handler('time_scalar_divide quotient error :: notify developer')

end function time_scalar_divide



function interval_alarm(time, time_interval, alarm, alarm_interval)
!-------------------------------------------------------------------------
!
! Supports a commonly used type of test on times for models.  Given the
! current time, and a time for an alarm, determines if this is the closest
! time to the alarm time given a time step of time_interval.  If this
! is the closest time (alarm - time <= time_interval/2), the function 
! returns true and the alarm is incremented by the alarm_interval.  Watch
! for problems if the new alarm time is less than time + time_interval

implicit none

type(time_type), intent(in)    :: time, time_interval, alarm_interval
type(time_type), intent(inout) :: alarm
logical                        :: interval_alarm

if((alarm - time) <= (time_interval / 2)) then
   interval_alarm = .TRUE.
   alarm = alarm + alarm_interval
else
   interval_alarm = .FALSE.
end if

end function interval_alarm



function repeat_alarm(time, alarm_frequency, alarm_length)
!--------------------------------------------------------------------------
!
! Repeat_alarm supports an alarm that goes off with alarm_frequency and
! lasts for alarm_length.  If the nearest occurence of an alarm time
! is less than half an alarm_length from the input time, repeat_alarm
! is true.  For instance, if the alarm_frequency is 1 day, and the 
! alarm_length is 2 hours, then repeat_alarm is true from time 2300 on 
! day n to time 0100 on day n + 1 for all n.

implicit none

type(time_type), intent(in) :: time, alarm_frequency, alarm_length
logical                     :: repeat_alarm

type(time_type) :: prev, next

prev = (time / alarm_frequency) * alarm_frequency
next = prev + alarm_frequency

if(time - prev <= alarm_length / 2 .or. next - time <= alarm_length / 2) then
   repeat_alarm = .TRUE.
else
   repeat_alarm = .FALSE.
endif

end function repeat_alarm


!=========================================================================
! CALENDAR OPERATIONS BEGIN HERE
!=========================================================================

subroutine set_calendar_type(type)

! Selects calendar for default mapping from time to date. 

implicit none

integer, intent(in) :: type

if(type <  0 .or. type > max_type) &
   call error_handler('Illegal type in set_calendar_type')
calendar_type = type

if(type == GREGORIAN) &
   call error_handler('set_calendar_type :: GREGORIAN CALENDAR not implemented')

end subroutine set_calendar_type



function get_calendar_type()
!------------------------------------------------------------------------
!
! Returns default calendar type for mapping from time to date.

implicit none

integer :: get_calendar_type

get_calendar_type = calendar_type

end function get_calendar_type



!========================================================================
! START OF get_date BLOCK
!========================================================================

subroutine get_date(time, year, month, day, hour, minute, second)
!------------------------------------------------------------------------
!
! Given a time, computes the corresponding date given the selected calendar

implicit none

type(time_type), intent(in)  :: time
integer,         intent(out) :: second, minute, hour, day, month, year

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   call get_date_thirty(time, year, month, day, hour, minute, second)
case(GREGORIAN)
   call get_date_gregorian(time, year, month, day, hour, minute, second)
case(JULIAN)
   call get_date_julian(time, year, month, day, hour, minute, second)
case(NO_LEAP)
   call get_date_no_leap(time, year, month, day, hour, minute, second)
case default
   call error_handler('Invalid calendar type in get_date')
end select
end subroutine get_date



subroutine get_date_gregorian(time, year, month, day, hour, minute, second)
!------------------------------------------------------------------------
!
! Computes date corresponding to time for gregorian calendar

implicit none

type(time_type), intent(in)  :: time
integer,         intent(out) :: second, minute, hour, day, month, year

integer :: m,t,dyear,dmonth,dday,nday,nleapyr,nfh,nhund,nfour
integer ndiy,nex,ibaseyr
logical :: leap

ibaseyr= 1601
! set nday initially to 109207 (the # of days from 1/1/1601 to 1/1/1900)
! 227 years of 365 days + 72 leap years
nday=109207

! time in seconds from base_year
t = time%seconds
nday=nday+t/(60*60*24)

! find number of four hundred year periods
nfh=nday/146097
nday=modulo(nday,146097)

! find number of hundred year periods
nhund= nday/36524
if(nhund.gt.3) then
  nhund=3
  nday=36524
else
  nday=modulo(nday,36524)
endif

! find number of four year periods
nfour=nday/1461
nday=modulo(nday,1461)
nex=nday/365
if(nex.gt.3) then
 nex=3
 nday=365
else
 nday=modulo(nday,365)
endif

! Is this a leap year? Gregorian calandar assigns each year evenly
! divisible by 4 that is not a century year unevenly divisible by 400
! as a leap-year. (i.e. 1700,1800,1900 are not leap-years, 2000 is)

leap=(nex.eq.3).and.((nfour.ne.24).or.(nhund.eq.3))
 if (leap) then
  ndiy=366
 else
  ndiy=365
 endif
year=ibaseyr+400*nfh+100*nhund+4*nfour+nex
nday=nday+1
! find month 
month=0
do m=1,12
 if (leap.and.(m.eq.2)) then
  if (nday.le. (days_per_month(2)+1)) then
   month = m
   go to 10
  else
   nday = nday - (days_per_month(2)+1)
   month = m
   t = t -  (60*60*24 * (days_per_month(2)+1))
  endif
 else 
  if (nday.le. days_per_month(m)) then
   month = m
   go to 10
  else
   nday = nday - days_per_month(m)
   month = m
   t = t -  (60*60*24 * days_per_month(month))
  endif
 endif
enddo
10 continue

! find day, hour,minute and second

dday   = t / (60*60*24)
day    = nday
t      = t - dday * (60*60*24)
hour   = t / (60 * 60)
t      = t - hour * (60 * 60)
minute = t / 60
second = t - 60 * minute

!if(leap) print*,'1:t,s,m,h,d,m,y=',time,second,minute,hour,day,month,year

end subroutine get_date_gregorian



subroutine get_date_julian(time, year, month, day, hour, minute, second)
!------------------------------------------------------------------------

! Base date for Julian calendar is year 1 with all multiples of 4 
! years being leap years.

implicit none

type(time_type), intent(in)  :: time
integer,         intent(out) :: second, minute, hour, day, month, year

integer :: m, t, nfour, nex, days_this_month
logical :: leap

! find number of four year periods; also get modulo number of days

nfour = time%days / (4 * 365 + 1) 
day   = modulo(time%days, (4 * 365 + 1))

! Find out what year in four year chunk
nex = day / 365
if(nex == 4) then
   nex = 3
   day = 366
else
   day=modulo(day, 365) + 1
endif

! Is this a leap year? 
leap = (nex == 3)

year = 1 + 4 * nfour + nex

! find month and day
do m = 1, 12
   month = m
   days_this_month = days_per_month(m)
   if(leap .and. m == 2) days_this_month = 29
   if(day <= days_this_month) exit
   day = day - days_this_month
end do

! find hour,minute and second

t      = time%seconds
hour   = t / (60 * 60)
t      = t - hour * (60 * 60)
minute = t / 60
second = t - 60 * minute

end subroutine get_date_julian



subroutine get_date_thirty(time, year, month, day, hour, minute, second)
!------------------------------------------------------------------------
!
! Computes date corresponding to time interval for 30 day months, 12
! month years.

implicit none

type(time_type), intent(in)  :: time
integer,         intent(out) :: second, minute, hour, day, month, year

integer :: t, dmonth, dyear

t      = time%days
dyear  = t / (30 * 12)
year   = dyear + 1
t      = t - dyear * (30 * 12)
dmonth = t / 30
month  = 1 + dmonth
day    = t -dmonth * 30 + 1

t      = time%seconds
hour   = t / (60 * 60) 
t      = t - hour * (60 * 60)
minute = t / 60
second = t - 60 * minute

end subroutine get_date_thirty



subroutine get_date_no_leap(time, year, month, day, hour, minute, second)
!------------------------------------------------------------------------
!
! Base date for no_leap calendar is year 1.

implicit none

type(time_type), intent(in)  :: time
integer,         intent(out) :: second, minute, hour, day, month, year

integer :: m, t

! get modulo number of days

year = time%days / 365 + 1
day  = modulo(time%days, 365) + 1

! find month and day

do m = 1, 12
   month = m
   if(day <= days_per_month(m)) exit
   day = day - days_per_month(m)
end do

! find hour,minute and second

t      = time%seconds
hour   = t / (60 * 60)
t      = t - hour * (60 * 60)
minute = t / 60
second = t - 60 * minute

end subroutine get_date_no_leap

!========================================================================
! END OF get_date BLOCK
!========================================================================
! START OF set_date BLOCK
!========================================================================

function set_date(year, month, day, hours, minutes, seconds)

! Given a date, computes the corresponding time given the selected
! date time mapping algorithm.  Note that it is possible to specify
! any number of illegal dates; these should be checked for and generate
! errors as appropriate.

implicit none

integer, intent(in)           :: day, month, year
integer, intent(in), optional :: seconds, minutes, hours
type(time_type)               :: set_date

integer :: oseconds, ominutes, ohours

! Missing optionals are set to 0
oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   set_date =    set_date_thirty(year, month, day, ohours, ominutes, oseconds)
case(GREGORIAN)
   set_date = set_date_gregorian(year, month, day, ohours, ominutes, oseconds)
case(JULIAN)
   set_date =    set_date_julian(year, month, day, ohours, ominutes, oseconds)
case(NO_LEAP)
   set_date =   set_date_no_leap(year, month, day, ohours, ominutes, oseconds)
case default
   call error_handler('Invalid calendar type in set_date')
end select
end function set_date



function set_date_gregorian(year, month, day, hours, minutes, seconds)
!------------------------------------------------------------------------
!
! Computes time corresponding to date for gregorian calendar.

implicit none

integer, intent(in)           :: day, month, year
integer, intent(in), optional :: seconds, minutes, hours
type(time_type)               :: set_date_gregorian

integer :: oseconds, ominutes, ohours
integer :: days, m, nleapyr
integer :: base_year = 1900
logical :: leap

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  ))   ohours = hours

! Need to check for bogus dates

if(    oseconds > 59 .or. oseconds < 0 .or. &
       ominutes > 59 .or. ominutes < 0 .or. &
       ohours   > 23 .or. ohours   < 0 .or. &
       day      > 31 .or. day      < 1 .or. &
       month    > 12 .or. month    < 1 .or. &
                          year     < base_year) &
      call error_handler('Invalid date in set_date_gregorian')

! QJLA ... above check does not work for Feb 30 ... etc.

! Is this a leap year? Gregorian calandar assigns each year evenly
! divisible by 4 that is not a century year unevenly divisible by 400
! as a leap-year. (i.e. 1700,1800,1900 are not leap-years, 2000 is)

  leap=(modulo(year,4).eq.0)
  if((modulo(year,100).eq.0).and.(modulo(year,400).ne.0))then
   leap=.false.
  endif

! compute number of leap years from base_year

nleapyr=((year-1)-base_year)/4-((year-1)-base_year)/100+((year-1)-1600)/400
days = 0
do m=1,month-1
 days = days + days_per_month(m)
 if(leap.and.m.eq.2)days=days+1
enddo

set_date_gregorian%seconds = oseconds + 60*(ominutes + 60*(ohours + 24*((day - 1) + &
        (days + 365*(year - base_year-nleapyr)+366*(nleapyr)))))

end function set_date_gregorian



function set_date_julian(year, month, day, hours, minutes, seconds)
!------------------------------------------------------------------------
!
! Returns time corresponding to date for julian calendar.

implicit none

integer, intent(in)           :: day, month, year
integer, intent(in), optional :: seconds, minutes, hours
type(time_type)               :: set_date_julian

integer :: oseconds, ominutes, ohours
integer :: ndays, m, nleapyr
logical :: leap

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours

! Need to check for bogus dates

if( oseconds > 59 .or. oseconds < 0 .or. &
    ominutes > 59 .or. ominutes < 0 .or. &
    ohours   > 23 .or. ohours   < 0 .or. &
                       day      < 1 .or. &
    month    > 12 .or. month    < 1 .or. &
                       year     < 1) &
   call error_handler('Invalid date in set_date_julian')

if(month /= 2 .and. day > days_per_month(month)) &
   call error_handler('Invalid day in set_date_julian')

! Is this a leap year? 
leap = (modulo(year,4) == 0)

! compute number of complete leap years from year 1
nleapyr = (year - 1) / 4

! Finish checking for day specication errors

if(month == 2 .and. (day > 29 .or. ((.not. leap) .and. day > 28))) &
   call error_handler('Invalid number of days in month 2 in set_date_julian')

ndays = 0
do m = 1, month - 1
   ndays = ndays + days_per_month(m)
   if(leap .and. m == 2) ndays = ndays + 1
enddo

set_date_julian%seconds = oseconds + 60 * (ominutes + 60 * ohours)
set_date_julian%days    = day -1 + ndays + 365*(year - nleapyr - 1) + 366*(nleapyr)

end function set_date_julian



function set_date_thirty(year, month, day, hours, minutes, seconds)
!------------------------------------------------------------------------
!
! Computes time corresponding to date for thirty day months.

implicit none

integer, intent(in)           :: day, month, year
integer, intent(in), optional :: seconds, minutes, hours
type(time_type)               :: set_date_thirty

integer :: oseconds, ominutes, ohours

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours

! Need to check for bogus dates

if( oseconds > 59 .or. oseconds < 0 .or. &
    ominutes > 59 .or. ominutes < 0 .or. &
    ohours   > 23 .or. ohours   < 0 .or. &
    day      > 30 .or. day      < 1 .or. &
    month    > 12 .or. month    < 1 .or. &
                       year     < 1 ) &
      call error_handler('Invalid date in set_date_thirty')

set_date_thirty%days    = (day - 1) + 30 * ((month - 1) + 12 * (year - 1))
set_date_thirty%seconds = oseconds + 60 * (ominutes + 60 * ohours)

end function set_date_thirty



function set_date_no_leap(year, month, day, hours, minutes, seconds)
!------------------------------------------------------------------------
!
! Computes time corresponding to date for fixed 365 day year calendar.

implicit none

integer, intent(in)           :: day, month, year
integer, intent(in), optional :: seconds, minutes, hours
type(time_type)               :: set_date_no_leap

integer :: oseconds, ominutes, ohours
integer :: ndays, m

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours

! Need to check for bogus dates

if( oseconds > 59 .or. oseconds < 0 .or. &
    ominutes > 59 .or. ominutes < 0 .or. &
    ohours   > 23 .or. ohours   < 0 .or. &
    day      > 31 .or. day      < 1 .or. &
    month    > 12 .or. month    < 1 .or. &
                       year     < 1 ) then
   call error_handler('Invalid date in set_date_no_leap')
endif

ndays = 0
do m = 1, month - 1
   ndays = ndays + days_per_month(m)
enddo

set_date_no_leap = set_time(oseconds + 60 * (ominutes + 60 * ohours), &
   day -1 + ndays + 365 * (year - 1))

end function set_date_no_leap

!=========================================================================
! END OF set_date BLOCK
!=========================================================================
! START OF increment_date BLOCK
!=========================================================================


function increment_date(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given a time and some date increment, computes a new time.  Depending
! on the mapping algorithm from date to time, it may be possible to specify
! undefined increments (i.e. if one increments by 68 days and 3 months in
! a Julian calendar, it matters which order these operations are done and
! we don't want to deal with stuff like that, make it an error).

implicit none

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: increment_date

integer :: oseconds, ominutes, ohours, odays, omonths, oyears

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours
odays    = 0; if(present(days   )) odays    = days
omonths  = 0; if(present(months )) omonths  = months
oyears   = 0; if(present(years  )) oyears   = years

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   increment_date = increment_thirty(time, oyears, omonths, odays, &
                                     ohours, ominutes, oseconds)
case(GREGORIAN)
   increment_date = increment_gregorian(time, oyears, omonths, odays, &
                                     ohours, ominutes, oseconds)
case(JULIAN)
   increment_date = increment_julian(time, oyears, omonths, odays, &
                                     ohours, ominutes, oseconds)
case(NO_LEAP)
   increment_date = increment_no_leap(time, oyears, omonths, odays, &
                                     ohours, ominutes, oseconds)
case default
   call error_handler('Invalid calendar type in increment_date')
end select
end function increment_date



function increment_gregorian(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given time and some date increment, computes new time for gregorian calendar.

implicit none

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: increment_gregorian

integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cminute, chour, cday, cmonth, cyear

call error_handler('increment_gregorian not implemented')

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours
odays    = 0; if(present(days   )) odays    = days
omonths  = 0; if(present(months )) omonths  = months
oyears   = 0; if(present(years  )) oyears   = years

! Increment must be positive definite

if( oseconds < 0 .or. ominutes < 0 .or. &
    ohours   < 0 .or. odays    < 0 .or. &
    omonths  < 0 .or. oyears   < 0) &
    call error_handler('Illegal increment in increment_gregorian')

! First convert time into date

call get_date_gregorian(time, cyear, cmonth, cday, chour, cminute, csecond)

! Add on the increments

csecond = csecond + oseconds
cminute = cminute + ominutes
chour   = chour   + ohours
cday    = cday    + odays
cmonth  = cmonth  + omonths
cyear   = cyear   + oyears

! Convert this back into a time

increment_gregorian = set_date_gregorian(cyear, cmonth, cday, chour, cminute, csecond)

end function increment_gregorian



function increment_julian(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given time and some date increment, computes new time for julian calendar.

implicit none

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: increment_julian

integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cminute, chour, cday, cmonth, cyear, dyear
type(time_type) :: t

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours
odays    = 0; if(present(days   )) odays    = days
omonths  = 0; if(present(months )) omonths  = months
oyears   = 0; if(present(years  )) oyears   = years

! Increment must be positive definite

if( oseconds < 0 .or. ominutes < 0 .or. &
    ohours   < 0 .or. odays    < 0 .or. &
    omonths  < 0 .or. oyears   < 0 ) &
    call error_handler('Illegal increment in increment_julian')

!  There are a number of other bad types of increments that should be
!  prohibited here; the addition is not associative
!  Easiest thing is to only let month and year be incremented by themselves
!  This is slight overkill since year is not really a problem.

! QJLA .... see above comment ... will year be a problem?

if( omonths /= 0 .and. &
    (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
     odays    /= 0 .or. oyears   /= 0)  ) call error_handler  &
   ('increment_julian: month must not be incremented with other units')

if( oyears /= 0 .and. &
    (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
     odays    /= 0 .or. omonths  /= 0) ) call error_handler &
   ('increment_julian: year must not be incremented with other units')

!  For non-month and non-year part can just use increment_thirty

t =  increment_thirty(time, 0, 0, odays, ohours, ominutes, oseconds)

!  For month or year increment, first convert to date

call get_date_julian(t, cyear, cmonth, cday, chour, cminute, csecond)

cmonth = cmonth + omonths
cyear  = cyear + oyears

! Check for months larger than 12 and fix

if(cmonth > 12) then
   dyear  = (cmonth - 1) / 12 
   cmonth = cmonth - 12 * dyear
   cyear  = cyear + dyear
end if

! Convert this back into a time

increment_julian = set_date_julian(cyear, cmonth, cday, chour, cminute, csecond)

end function increment_julian



function increment_thirty(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given a time and some date increment, computes new time for thirty-day months.

implicit none

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: increment_thirty

integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cday

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours
odays    = 0; if(present(days   )) odays    = days
omonths  = 0; if(present(months )) omonths  = months
oyears   = 0; if(present(years  )) oyears   = years

! Increment must be positive definite

if( oseconds < 0 .or. ominutes < 0 .or. &
    ohours   < 0 .or. odays    < 0 .or. &
    omonths  < 0 .or. oyears   < 0 ) &
    call error_handler('Illegal increment in increment_thirty')

! Do increment to seconds portion first

csecond = oseconds + 60 * (ominutes + 60 * ohours)
cday    = odays    + 30 * (omonths + 12 * oyears)
increment_thirty = increment_time(time, csecond, cday)

end function increment_thirty



function increment_no_leap(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given time and some date increment, computes new time for julian calendar.

implicit none

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: increment_no_leap

integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cminute, chour, cday, cmonth, cyear, dyear
type(time_type) :: t

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours
odays    = 0; if(present(days   )) odays    = days
omonths  = 0; if(present(months )) omonths  = months
oyears   = 0; if(present(years  )) oyears   = years

! Increment must be positive definite

if( oseconds < 0 .or. ominutes < 0 .or. &
    ohours   < 0 .or. odays    < 0 .or. &
    omonths  < 0 .or. oyears   < 0) &
    call error_handler('Illegal increment in increment_no_leap')

!  There are a number of other bad types of increments that should be
!  prohibited here; the addition is not associative
!  Easiest thing is to only let month and year be incremented by themselves
!  This is slight overkill since year is not really a problem.

if(omonths /= 0 .and. &
   (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
    odays    /= 0 .or. oyears   /= 0)) call error_handler &
   ('increment_no_leap: month must not be incremented with other units')

if(oyears /= 0 .and. &
   (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
    odays    /= 0 .or. omonths  /= 0)) call error_handler &
   ('increment_no_leap: year must not be incremented with other units')

!  For non-month and non-year part can just use increment_thirty

t =  increment_thirty(time, 0, 0, odays, ohours, ominutes, oseconds)

!  For month or year increment, first convert to date

call get_date_no_leap(t, cyear, cmonth, cday, chour, cminute, csecond)
cmonth = cmonth + omonths
cyear  = cyear + oyears

! Check for months larger than 12 and fix
if(cmonth > 12) then
   dyear  = (cmonth - 1) / 12 
   cmonth = cmonth - 12 * dyear
   cyear  = cyear + dyear
end if

! Convert this back into a time

increment_no_leap = set_date_no_leap(cyear, cmonth, cday, chour, cminute, csecond)

end function increment_no_leap


!=========================================================================
! END OF increment_date BLOCK
!=========================================================================
! START OF decrement_date BLOCK
!=========================================================================


function decrement_date(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given a time and some date decrement, computes a new time.  Depending
! on the mapping algorithm from date to time, it may be possible to specify
! undefined decrements (i.e. if one decrements by 68 days and 3 months in
! a Julian calendar, it matters which order these operations are done and
! we don't want to deal with stuff like that, make it an error).

implicit none

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: decrement_date

integer :: oseconds, ominutes, ohours, odays, omonths, oyears

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours
odays    = 0; if(present(days   )) odays    = days
omonths  = 0; if(present(months )) omonths  = months
oyears   = 0; if(present(years  )) oyears   = years

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   decrement_date = decrement_thirty(time, oyears, omonths, odays, &
                       ohours, ominutes, oseconds)
case(GREGORIAN)
   decrement_date = decrement_gregorian(time, oyears, omonths, odays, &
                       ohours, ominutes, oseconds)
case(JULIAN)
   decrement_date = decrement_julian(time, oyears, omonths, odays, &
                       ohours, ominutes, oseconds)
case(NO_LEAP)
   decrement_date = decrement_no_leap(time, oyears, omonths, odays, &
                       ohours, ominutes, oseconds)
case default
   call error_handler('Invalid calendar type in decrement_date')
end select

end function decrement_date



function decrement_gregorian(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given time and some date decrement, computes new time for gregorian calendar.

implicit none

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: decrement_gregorian

integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cminute, chour, cday, cmonth, cyear

call error_handler('decrement_gregorian not implemented')

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours
odays    = 0; if(present(days   )) odays    = days
omonths  = 0; if(present(months )) omonths  = months
oyears   = 0; if(present(years  )) oyears   = years

! Decrement must be positive definite

! this block needs work ... TJH ....

if( oseconds < 0 .or. ominutes < 0 .or. &
    ohours   < 0 .or. odays    < 0 .or. &
    omonths  < 0 .or. oyears < 0) &
    call error_handler('Illegal decrement in decrement_gregorian')

! First convert time into date

call get_date_gregorian(time, cyear, cmonth, cday, chour, cminute, csecond)

! Remove the increments

csecond = csecond - oseconds
cminute = cminute - ominutes
chour   = chour   - ohours
cday    = cday    - odays
cmonth  = cmonth  - omonths
cyear  = cyear    - oyears

! Convert this back into a time

decrement_gregorian =  set_date_gregorian(cyear, cmonth, cday, chour, cminute, csecond)

end function decrement_gregorian



function decrement_julian(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given time and some date decrement, computes new time for julian calendar.

implicit none

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: decrement_julian

integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cminute, chour, cday, cmonth, cyear
type(time_type) :: t

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours
odays    = 0; if(present(days   )) odays    = days
omonths  = 0; if(present(months )) omonths  = months
oyears   = 0; if(present(years  )) oyears   = years

! Increment must be positive definite

if( oseconds < 0 .or. ominutes < 0 .or. ohours < 0 .or. &
    odays    < 0 .or. omonths  < 0 .or. oyears < 0) &
   call error_handler('Illegal increment in decrement_julian')

!  There are a number of other bad types of decrements that should be
!  prohibited here; the subtraction is not associative
!  Easiest thing is to only let month and year be decremented by themselves
!  This is slight overkill since year is not really a problem.

if(omonths /= 0 .and. (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
   odays  /= 0 .or. oyears /= 0)) call error_handler &
   ('decrement_julian: month must not be decremented with other units')

if(oyears /= 0 .and. (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
   odays  /= 0 .or. omonths /= 0)) call error_handler &
   ('decrement_julian: year must not be decremented with other units')

!  For non-month and non-year can just use decrement_thirty

t = decrement_thirty(time, 0, 0, odays, ohours, ominutes, oseconds)

!  For month or year decrement, first convert to date

call get_date_julian(t, cyear, cmonth, cday, chour, cminute, csecond)
cmonth = cmonth - omonths
cyear  = cyear - oyears

! Check for months less than 12 and fix
if(cmonth < 1) then
   cyear = cyear - 1 + (cmonth) / 12
   cmonth = cmonth - 12 * ( -1 + (cmonth) / 12)
end if

! Check for negative years
if(cyear < 1) call error_handler('Illegal date results in decrement_julian')

! Convert this back into a time
decrement_julian = set_date_julian(cyear, cmonth, cday, chour, cminute, csecond)

end function decrement_julian



function decrement_thirty(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given a time and some date decrement, computes new time for thirty day months.

implicit none

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: decrement_thirty

integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cday

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours
odays    = 0; if(present(days   )) odays    = days
omonths  = 0; if(present(months )) omonths  = months
oyears   = 0; if(present(years  )) oyears   = years

! Increment must be positive definite

if( oseconds < 0 .or. ominutes < 0 .or. ohours < 0 .or. &
    odays    < 0 .or. omonths  < 0 .or. oyears < 0) &
    call error_handler('Illegal decrement in decrement_thirty')

csecond = oseconds + 60 * (ominutes + 60 * ohours)
cday    = odays + 30 * (omonths + 12 * oyears)
decrement_thirty = decrement_time(time, csecond, cday)

end function decrement_thirty



function decrement_no_leap(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given time and some date decrement, computes new time for julian calendar.

implicit none

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: decrement_no_leap

integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cminute, chour, cday, cmonth, cyear
type(time_type) :: t

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours
odays    = 0; if(present(days   )) odays    = days
omonths  = 0; if(present(months )) omonths  = months
oyears   = 0; if(present(years  )) oyears   = years

! Increment must be positive definite

if( oseconds < 0 .or. ominutes < 0 .or. ohours < 0 .or. &
    odays    < 0 .or. omonths  < 0 .or. oyears < 0 ) &
    call error_handler('Illegal increment in decrement_no_leap')

!  There are a number of other bad types of decrements that should be
!  prohibited here; the subtraction is not associative
!  Easiest thing is to only let month and year be decremented by themselves
!  This is slight overkill since year is not really a problem.

if(omonths /= 0 .and. (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
   odays   /= 0 .or. oyears /= 0)) call error_handler &
   ('decrement_no_leap: month must not be decremented with other units')

if(oyears /= 0 .and. (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
   odays  /= 0 .or. omonths /= 0)) call error_handler &
   ('decrement_no_leap: year must not be decremented with other units')

!  For non-month and non-year can just use decrement_thirty

t = decrement_thirty(time, 0, 0, odays, ohours, ominutes, oseconds)

!  For month or year decrement, first convert to date

call get_date_no_leap(t, cyear, cmonth, cday, chour, cminute, csecond)
cmonth = cmonth - omonths
cyear  = cyear  - oyears

! Check for months less than 12 and fix
if(cmonth < 1) then
   cyear  = cyear - 1 + (cmonth) / 12
   cmonth = cmonth - 12 * ( -1 + (cmonth) / 12)
end if

! Check for negative years

if(cyear < 1) call error_handler('Illegal date results in decrement_no_leap')

! Convert this back into a time

decrement_no_leap = set_date_no_leap(cyear, cmonth, cday, chour, cminute, csecond)

end function decrement_no_leap


!=========================================================================
! END OF decrement_date BLOCK
!=========================================================================
! START days_in_month BLOCK
!=========================================================================


function days_in_month(time)
!--------------------------------------------------------------------------
!
! Given a time, computes the corresponding date given the selected
! date time mapping algorithm

implicit none

type(time_type), intent(in) :: time
integer                     :: days_in_month

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   days_in_month = days_in_month_thirty(time)
case(GREGORIAN)
   days_in_month = days_in_month_gregorian(time)
case(JULIAN)
   days_in_month = days_in_month_julian(time)
case(NO_LEAP)
   days_in_month = days_in_month_no_leap(time)
case default
   call error_handler('Invalid calendar type in days_in_month')
end select
end function days_in_month



function days_in_month_gregorian(time)
!--------------------------------------------------------------------------
!
! Returns the number of days in a gregorian month.

implicit none

type(time_type), intent(in) :: time
integer                     :: days_in_month_gregorian

call error_handler('days_in_month_gregorian not implemented')
days_in_month_gregorian = -1

end function days_in_month_gregorian



function days_in_month_julian(time)
!--------------------------------------------------------------------------
!
! Returns the number of days in a julian month.

implicit none

type(time_type), intent(in) :: time
integer                     :: days_in_month_julian

integer :: seconds, minutes, hours, day, month, year

call get_date_julian(time, year, month, day, hours, minutes, seconds)
days_in_month_julian = days_per_month(month)
if(leap_year_julian(time) .and. month == 2) days_in_month_julian = 29

end function days_in_month_julian



function days_in_month_thirty(time)
!--------------------------------------------------------------------------
! Returns the number of days in a thirty day month (needed for transparent
! changes to calendar type).

implicit none

type(time_type), intent(in) :: time
integer                     :: days_in_month_thirty

days_in_month_thirty = 30

end function days_in_month_thirty



function days_in_month_no_leap(time)
!--------------------------------------------------------------------------
! Returns the number of days in a 365 day year month.

implicit none

type(time_type), intent(in) :: time
integer                     :: days_in_month_no_leap

integer :: seconds, minutes, hours, day, month, year

call get_date_no_leap(time, year, month, day, hours, minutes, seconds)
days_in_month_no_leap= days_per_month(month)

end function days_in_month_no_leap


!==========================================================================
! END OF days_in_month BLOCK
!==========================================================================
! START OF leap_year BLOCK
!==========================================================================


function leap_year(time)
!--------------------------------------------------------------------------
! Is this date in a leap year for default calendar?

implicit none

type(time_type), intent(in) :: time
logical                     :: leap_year

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   leap_year = leap_year_thirty(time)
case(GREGORIAN)
   leap_year = leap_year_gregorian(time)
case(JULIAN)
   leap_year = leap_year_julian(time)
case(NO_LEAP)
   leap_year = leap_year_no_leap(time)
case default
   call error_handler('Invalid calendar type in leap_year')
end select
end function leap_year



function leap_year_gregorian(time)
!--------------------------------------------------------------------------
!
! Is this a leap year for gregorian calendar?

implicit none

type(time_type), intent(in) :: time
logical                     :: leap_year_gregorian

call error_handler('leap_year_gregorian not implemented')
leap_year_gregorian = .FALSE.

end function leap_year_gregorian



function leap_year_julian(time)
!--------------------------------------------------------------------------
!
! Is this a leap year for julian calendar?

implicit none

type(time_type), intent(in) :: time
logical                     :: leap_year_julian

integer :: seconds, minutes, hours, day, month, year

call get_date(time, year, month, day, hours, minutes, seconds)
leap_year_julian = ((year / 4 * 4) == year)

! QJLA this cannot be right .... unless every fourth year is a leap ...

end function leap_year_julian



function leap_year_thirty(time)
!--------------------------------------------------------------------------
!
! No leap years in thirty day months, included for transparency. 

implicit none

type(time_type), intent(in) :: time
logical                     :: leap_year_thirty

leap_year_thirty = .FALSE.

end function leap_year_thirty



function leap_year_no_leap(time)
!--------------------------------------------------------------------------
!
! Another tough one; no leap year returns false for leap year inquiry.

implicit none

type(time_type), intent(in) :: time
logical                     :: leap_year_no_leap

leap_year_no_leap = .FALSE.

end function leap_year_no_leap


!==========================================================================
! END OF leap_year BLOCK
!==========================================================================
! START OF length_of_year BLOCK
!==========================================================================


function length_of_year()
!--------------------------------------------------------------------------
!
! What is the length of the year for the default calendar type

implicit none

type(time_type) :: length_of_year

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   length_of_year = length_of_year_thirty()
case(GREGORIAN)
   length_of_year = length_of_year_gregorian()
case(JULIAN)
   length_of_year = length_of_year_julian()
case(NO_LEAP)
   length_of_year = length_of_year_no_leap()
case default
   call error_handler('Invalid calendar type in length_of_year')
end select
end function length_of_year



function length_of_year_thirty()
!--------------------------------------------------------------------------
!

implicit none

type(time_type) :: length_of_year_thirty

length_of_year_thirty = set_time(0, 360)

end function length_of_year_thirty



function length_of_year_gregorian()
!---------------------------------------------------------------------------

implicit none

type(time_type) :: length_of_year_gregorian

length_of_year_gregorian = set_time(0, 0)

call error_handler('length_of_year_gregorian not implemented')

end function length_of_year_gregorian



function length_of_year_julian()
!--------------------------------------------------------------------------

implicit none

type(time_type) :: length_of_year_julian

length_of_year_julian = set_time((24 / 4) * 60 * 60, 365)

end function length_of_year_julian



function length_of_year_no_leap()
!--------------------------------------------------------------------------

implicit none

type(time_type) :: length_of_year_no_leap

length_of_year_no_leap = set_time(0, 365)

end function length_of_year_no_leap


!==========================================================================
! END OF length_of_year BLOCK
!==========================================================================
! START OF days_in_year BLOCK
!==========================================================================


function days_in_year(time)
!--------------------------------------------------------------------------
!
! What is the number of days in this year for the default calendar type

implicit none

type(time_type), intent(in) :: time
integer                     :: days_in_year

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   days_in_year = days_in_year_thirty(time)
case(GREGORIAN)
   days_in_year = days_in_year_gregorian(time)
case(JULIAN)
   days_in_year = days_in_year_julian(time)
case(NO_LEAP)
   days_in_year = days_in_year_no_leap(time)
case default
   call error_handler('Invalid calendar type in days_in_year')
end select
end function days_in_year



function days_in_year_thirty(time)
!--------------------------------------------------------------------------

implicit none

type(time_type), intent(in) :: time
integer                     :: days_in_year_thirty

days_in_year_thirty = 360

end function days_in_year_thirty



function days_in_year_gregorian(time)
!---------------------------------------------------------------------------

implicit none

type(time_type), intent(in) :: time
integer                     :: days_in_year_gregorian

days_in_year_gregorian = 0

call error_handler('days_in_year_gregorian not implemented')

end function days_in_year_gregorian



function days_in_year_julian(time)
!--------------------------------------------------------------------------

implicit none

type(time_type), intent(in) :: time
integer                     :: days_in_year_julian

if(leap_year_julian(time)) then
   days_in_year_julian = 366
else
   days_in_year_julian = 365
endif

end function days_in_year_julian



function days_in_year_no_leap(time)
!--------------------------------------------------------------------------

implicit none

type(time_type), intent(in) :: time
integer                     :: days_in_year_no_leap

days_in_year_no_leap = 365

end function days_in_year_no_leap


!==========================================================================
! END OF days_in_year BLOCK
!==========================================================================


function month_name(n)

! Returns character string associated with a month, for now, all calendars
! have 12 months and will return standard names.

integer, intent(in) :: n
character (len=9)   :: month_name

character(len=9), dimension(12) :: months = (/'January  ','February ','March    ',&
                                              'April    ','May      ','June     ',&
                                              'July     ','August   ','September',&
                                              'October  ','November ','December '/) 

if(n < 1 .or. n > 12) call error_handler('Illegal month index')

month_name = months(n)

end function month_name

!==========================================================================



subroutine error_handler(s)

implicit none

character (*), intent(in) :: s

! Stub until module for error_handler available
write(*, *) 'ERROR: In time_manager.f90: ', s
stop

!   call error_mesg ('time_manager', s, FATAL)

end subroutine error_handler



subroutine time_manager_init ( )
!------------------------------------------------------------------------
!
! initialization routine
! this routine does not have to be called, all it does is write
! the version information to the log file

  if (.not.do_init) return  ! silent return if already called

!  call write_version_number (version, tagname)
!  do_init = .false.

end subroutine time_manager_init



subroutine print_time (time,str,unit)
!------------------------------------------------------------------------

type(time_type)  , intent(in)           :: time
character (len=*), intent(in), optional :: str
integer          , intent(in), optional :: unit

integer           :: s,d, ns,nd, unit_in
character(len=13) :: fmt

! prints the time to standard output (or optional unit) as days and seconds
! NOTE: there is no check for PE number

!!!  unit_in = stdout()

  if (present(unit)) unit_in = unit

  call get_time (time,s,d)

! format output
! get number of digits for days and seconds strings

   nd = int(log10(real(max(1,d))))+1
   ns = int(log10(real(max(1,s))))+1
   write (fmt,10) nd, ns
10 format ('(a,i',i2.2,',a,i',i2.2,')')

  if (present(str)) then
     write (unit_in,fmt) trim(str)//' day=', d, ', sec=', s
  else
     write (unit_in,fmt)       'TIME: day=', d, ', sec=', s
  endif

end subroutine print_time



subroutine print_date (time,str,unit)
!------------------------------------------------------------------------

type(time_type)  , intent(in)           :: time
character (len=*), intent(in), optional :: str
integer          , intent(in), optional :: unit

integer          :: y,mo,d,h,m,s, unit_in
character(len=9) :: mon

! prints the time to standard output (or optional unit) as a date
! NOTE: there is no check for PE number

!!!  unit_in = stdout()

  if (present(unit)) unit_in = unit

  call get_date (time,y,mo,d,h,m,s)
  mon = month_name(mo)
  if (present(str)) then
     write (unit_in,10) trim(str)//' ', y,mon(1:3),' ',d,' ',h,':',m,':',s
  else
     write (unit_in,10)       'DATE: ', y,mon(1:3),' ',d,' ',h,':',m,':',s
  endif
10 format (a,i4,1x,a3,4(a1,i2.2))

end subroutine print_date



function read_time(file)
!--------------------------------------------------------------------------------
!

implicit none

integer, intent(in) :: file
type(time_type)     :: read_time

integer :: secs, days

read(file, *) secs, days
read_time = set_time(secs, days)

end function read_time



subroutine write_time(file, time)
!------------------------------------------------------------------------

implicit none

integer, intent(in) :: file

integer :: secs, days

type(time_type), intent(in) :: time

call get_time(time, secs, days)
write(file, *) secs, days

end subroutine write_time



function nc_append_time(ncFileID, time) result(len)
!------------------------------------------------------------------------
! The current time is appended to the "time" coordinate variable.
! The new length of the "time" variable is returned.
! 
! This REQUIRES that "time" is a coordinate variable AND it is the
! unlimited dimension. If not ... bad things happen.
!
! TJH Wed Aug 28 15:40:25 MDT 2002

use typeSizes
use netcdf
implicit none

integer,         intent(in) :: ncFileID
type(time_type), intent(in) :: time
integer                     :: len

integer  :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer  :: TimeVarID, i,ierr
integer  :: secs, days
real(r8) :: r8time         ! gets promoted to nf90_double ...

character(len=NF90_MAX_NAME)          :: varname
integer                               :: xtype, ndims, nAtts
integer, dimension(NF90_MAX_VAR_DIMS) :: dimids

len = -1 ! assume a bad termination

call check(NF90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(NF90_Inq_Varid(ncFileID, "time", TimeVarID))
call check(NF90_Inquire_Variable(ncFileID, TimeVarID, varname, xtype, ndims, dimids, nAtts))

if ( ndims /= 1 ) then
   write(*,*)'Warning:nc_append_time: "time" expected to be rank-1' 
endif

if ( dimids(1) /= unlimitedDimID ) then
   write(*,*)'Warning:nc_append_time: no idea what you are trying to pull here ...'
endif

call check(NF90_Inquire_Dimension(ncFileID, unlimitedDimID, varname, len ))

write(*,*)'Current length of ',trim(adjustl(varname)),' time variable is ',len  ! DEBUG

call get_time(time, secs, days)
r8time = days + secs/86400.0_r8    ! time base is "days since ..."

len = len + 1
call check(nf90_put_var(ncFileID, TimeVarID, r8time, start=(/ len /) ))

contains

  ! Internal subroutine - checks error status after each netcdf, prints
  !                       text message each time an error code is returned.
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
    end if
  end subroutine check

end function nc_append_time



end module time_manager_mod
