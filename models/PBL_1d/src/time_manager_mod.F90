! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module time_manager_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use     types_mod, only : missing_i, digits12
use utilities_mod, only : error_handler, E_DBG, E_MSG, E_WARN, E_ERR, &
                          register_module, dump_unit_attributes

implicit none
private

!====================================================================
! This module works best when the real variables have more than 7 
! significant digits, so this module uses the 'digits12' parameter 
! defined in types_mod ... 
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
public :: time_type

! Operators defined on time_type
public :: operator(+),  operator(-),   operator(*),   operator(/),  &
          operator(>),  operator(>=),  operator(==),  operator(/=), &
          operator(<),  operator(<=),  operator(//)

! Subroutines and functions operating on time_type
public :: set_time, set_time_missing, increment_time, decrement_time, get_time
public :: interval_alarm, repeat_alarm

! List of available calendar types
!!! NO_LEAP changed to NOLEAP for some weird FMS Havana compliance
public :: THIRTY_DAY_MONTHS,    JULIAN,    GREGORIAN,  NOLEAP,   NO_CALENDAR

! Subroutines and functions involving relations between time and calendar
public :: set_calendar_type, get_calendar_type
public :: set_date,       set_date_gregorian,         set_date_julian, &
                          set_date_thirty,            set_date_no_leap
public :: get_date,       get_date_gregorian,         get_date_julian, &
                          get_date_thirty,            get_date_no_leap
public :: increment_date, increment_gregorian,        increment_julian, &
                          increment_thirty,           increment_no_leap
public :: decrement_date, decrement_gregorian,        decrement_julian, &
                          decrement_thirty,           decrement_no_leap
public :: days_in_month,  days_in_month_gregorian,    days_in_month_julian, &
                          days_in_month_no_leap,      days_in_month_thirty
public :: leap_year,      leap_year_gregorian,        leap_year_julian, &
                          leap_year_no_leap,          leap_year_thirty
public :: length_of_year, length_of_year_thirty,      length_of_year_julian, &
                          length_of_year_gregorian,   length_of_year_no_leap
public :: days_in_year,   days_in_year_thirty,        days_in_year_julian, &
                          days_in_year_gregorian,     days_in_year_no_leap
public :: month_name

public :: julian_day

! Subroutines and functions for basic I/O
public :: time_manager_init, print_time, print_date
public :: write_time, read_time, interactive_time

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

! Global data to define calendar type
integer, parameter :: THIRTY_DAY_MONTHS = 1,      JULIAN = 2, &
                      GREGORIAN = 3,              NOLEAP = 4, &
                      NO_CALENDAR = 0
! HAMMERING DEFAULT CALENDAR TYPE --- MUST FIX FOR REAL --- TJH
! integer, private :: calendar_type = GREGORIAN, max_type = 4
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

logical, save :: module_initialized = .false.

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

integer            :: days_in
character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

days_in = 0;  if (present(days)) days_in = days

! Negative time offset is illegal

write(errstring,*)'seconds, days are ',seconds, days_in,' cannot be negative'
if(seconds < 0 .or. days_in < 0) &
   call error_handler(E_ERR,'set_time',errstring,source,revision,revdate)

! Make sure seconds greater than a day are fixed up

set_time%seconds = seconds - seconds / (60*60*24) * (60*60*24)

! Check for overflow on days before doing operation

write(errstring,*)'seconds is ',seconds,' overflowing conversion to days'
if(seconds / (60*60*24)  >= huge(days_in) - days_in) &
   call error_handler(E_ERR,'set_time',errstring,source,revision,revdate)

set_time%days = days_in + seconds / (60*60*24)

end function set_time



function set_time_missing()
!---------------------------------------------------------------------------
!
! Returns a time interval corresponding to this number of days and seconds.
! The arguments must not be negative but are otherwise unrestricted.

type(time_type) :: set_time_missing

set_time_missing%seconds = missing_i
set_time_missing%days    = missing_i

end function set_time_missing



subroutine get_time(time, seconds, days)
!---------------------------------------------------------------------------
!
! Returns days and seconds ( < 86400 ) corresponding to a time.

implicit none

type(time_type), intent(in)            :: time
integer,         intent(out)           :: seconds
integer,         intent(out), optional :: days

character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

seconds = time%seconds
if (present(days)) then
  days = time%days
else

  write(errstring,*)'seconds is ',seconds,' overflowing conversion to days'

  if (time%days > (huge(seconds) - seconds)/(60*60*24)) &
     call error_handler(E_ERR,'get_time',errstring,source,revision,revdate)

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

integer            :: days_in
character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

days_in = 0;  if (present(days)) days_in = days

! Increment must be positive definite

if(seconds < 0 .or. days_in < 0) then
   write(errstring,*)'Negative increment (',seconds, days,') (seconds,days) not allowed'
   call error_handler(E_ERR,'increment_time',errstring,source,revision,revdate)
endif

! Watch for immediate overflow on days or seconds

if(days_in >= huge(days_in) - time%days) then
   write(errstring,*)'integer overflow (',days_in, time%days,') in days'
   call error_handler(E_ERR,'increment_time',errstring,source,revision,revdate)
endif

if(seconds >= huge(seconds) - time%seconds) then
   write(errstring,*)'integer overflow (',seconds, time%seconds,') in days'
   call error_handler(E_ERR,'increment_time',errstring,source,revision,revdate)
endif

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

integer            :: cseconds, cdays
character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

cdays = 0;  if (present(days)) cdays = days

! Decrement must be positive definite

if(seconds < 0 .or. cdays < 0) then
   write(errstring,*)'Negative decrement (',seconds,cdays,') (seconds,days) not allowed.'
   call error_handler(E_ERR,'decrement_time',errstring,source,revision,revdate)
endif

cseconds = time%seconds - seconds
cdays    = time%days - cdays

! Borrow if needed

if(cseconds < 0) then
   cdays = cdays - 1 + (cseconds + 1) / (60*60*24)
   cseconds = cseconds - (60*60*24) * (-1 + (cseconds + 1) / (60*60*24))
end if

! Check for illegal negative time

if(cdays < 0) then
   write(errstring,*)'Negative time result (',cdays,') (days) not allowed.'
   call error_handler(E_ERR,'decrement_time',errstring,source,revision,revdate)
endif

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

if ( .not. module_initialized ) call time_manager_init

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

if ( .not. module_initialized ) call time_manager_init

time_ge = (time_gt(time1, time2) .or. time_eq(time1, time2))

end function time_ge



function time_lt(time1, time2)
!--------------------------------------------------------------------------
!
! Returns true if time1 < time2

implicit none

type(time_type), intent(in) :: time1, time2
logical                     :: time_lt

if ( .not. module_initialized ) call time_manager_init

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

if ( .not. module_initialized ) call time_manager_init

time_le = (time_lt(time1, time2) .or. time_eq(time1, time2))

end function time_le



function time_eq(time1, time2)
!--------------------------------------------------------------------------
!
! Returns true if time1 == time2

implicit none

type(time_type), intent(in) :: time1, time2
logical                     :: time_eq

if ( .not. module_initialized ) call time_manager_init

time_eq = (time1%seconds == time2%seconds .and. time1%days == time2%days)

end function time_eq



function time_ne(time1, time2)
!--------------------------------------------------------------------------
!
! Returns true if time1 /= time2

implicit none

type(time_type), intent(in) :: time1, time2
logical                     :: time_ne

if ( .not. module_initialized ) call time_manager_init

time_ne = (.not. time_eq(time1, time2))

end function time_ne



function time_plus(time1, time2)
!-------------------------------------------------------------------------
!
! Returns sum of two time_types

implicit none

type(time_type), intent(in) :: time1, time2
type(time_type)             :: time_plus

if ( .not. module_initialized ) call time_manager_init

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

if ( .not. module_initialized ) call time_manager_init

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

integer            :: days, seconds
real(digits12)     :: sec_prod
character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

! Multiplying here in a reasonable fashion to avoid overflow is tricky
! Could multiply by some large factor n, and seconds could be up to 86399
! Need to avoid overflowing integers and wrapping around to negatives

sec_prod = real(time%seconds,digits12) * real(n,digits12)

! If sec_prod is large compared to precision of double precision, things
! can go bad.  Need to warn and abort on this.

if(sec_prod /= 0.0_digits12) then
   if(log10(sec_prod) > precision(sec_prod) - 3) then
      write(errstring,*)'Insufficient precision to handle scalar product.'
      call error_handler(E_ERR,'time_scalar_mult',errstring,source,revision,revdate) 
   endif
end if

days    =  floor(sec_prod / (24.0_digits12 * 60.0_digits12 * 60.0_digits12))
seconds = sec_prod - days * (24.0_digits12 * 60.0_digits12 * 60.0_digits12)

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

if ( .not. module_initialized ) call time_manager_init

scalar_time_mult = time_scalar_mult(time, n)

end function scalar_time_mult



function time_divide(time1, time2)
!-------------------------------------------------------------------------
!
! Returns the largest integer, n, for which time1 >= time2 * n.

implicit none

type(time_type), intent(in) :: time1, time2
integer                     :: time_divide

real(digits12)     :: d1, d2
character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

! Convert time intervals to floating point days; risky for general performance?

d1 = time1%days * (60.0_digits12 * 60.0_digits12 * 24.0_digits12) + (time1%seconds*1.0_digits12)
d2 = time2%days * (60.0_digits12 * 60.0_digits12 * 24.0_digits12) + (time2%seconds*1.0_digits12) 

! Get integer quotient of this, check carefully to avoid round-off problems.

time_divide = floor(d1 / d2)

! Verify time_divide*time2 is <= time1 and (time_divide + 1)*time2 is > time1

if(time_divide * time2 > time1 .or. (time_divide + 1) * time2 <= time1) then
   write(errstring,*)'quotient error time1,time2 = ',time1, time2
   call error_handler(E_ERR,'time_divide',errstring,source,revision,revdate)
endif

end function time_divide



function time_real_divide(time1, time2)
!-------------------------------------------------------------------------
!
! Returns the double precision quotient of two times

implicit none

type(time_type), intent(in) :: time1, time2
real(digits12)              :: time_real_divide

real(digits12) :: d1, d2

if ( .not. module_initialized ) call time_manager_init

! Convert time intervals to floating point days; risky for general performance?

d1 = time1%days * (60.0_digits12 * 60.0_digits12 * 24.0_digits12) + (time1%seconds*1.0_digits12)
d2 = time2%days * (60.0_digits12 * 60.0_digits12 * 24.0_digits12) + (time2%seconds*1.0_digits12) 

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

real(digits12)  :: d, div
integer         :: days, seconds
type(time_type) :: prod1, prod2
character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

! Convert time interval to floating point days; risky for general performance?

d   = time%days * (60.0_digits12 * 60.0_digits12 * 24.0_digits12) + real(time%seconds,digits12)
div = d / real(n,digits12)

days    =  floor(div / (60.0_digits12 * 60.0_digits12 * 24.0_digits12))
seconds = div - days * (60.0_digits12 * 60.0_digits12 * 24.0_digits12)
time_scalar_divide = set_time(seconds, days)

! Need to make sure that roundoff isn't killing this

prod1 = n * time_scalar_divide
prod2 = n * (increment_time(time_scalar_divide, 1, 0)) 

if(prod1 > time .or. prod2 <= time) then
   write(errstring,*)'quotient error: ',n,time,prod1,prod2
   call error_handler(E_ERR,'time_scalar_divide',errstring,source,revision,revdate)
endif

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

if ( .not. module_initialized ) call time_manager_init

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

if ( .not. module_initialized ) call time_manager_init

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

character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

if(type <  0 .or. type > max_type) then
   write(errstring,*)'Illegal type ',type,' must be > 0 or < ',max_type
   call error_handler(E_ERR,'set_calendar_type',errstring,source,revision,revdate)
endif
calendar_type = type

! GREGORIAN Calendar only partially implemented; get and set_date work
!if(type == GREGORIAN) &
!   call error_handler('set_calendar_type :: GREGORIAN CALENDAR not implemented')

end subroutine set_calendar_type



function get_calendar_type()
!------------------------------------------------------------------------
!
! Returns default calendar type for mapping from time to date.

implicit none

integer :: get_calendar_type

if ( .not. module_initialized ) call time_manager_init

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

character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   call get_date_thirty(time, year, month, day, hour, minute, second)
case(GREGORIAN)
   call get_date_gregorian(time, year, month, day, hour, minute, second)
case(JULIAN)
   call get_date_julian(time, year, month, day, hour, minute, second)
case(NOLEAP)
   call get_date_no_leap(time, year, month, day, hour, minute, second)
case default
   write(errstring,*)'type is ',calendar_type,' must be one of ', &
                      THIRTY_DAY_MONTHS,GREGORIAN,JULIAN,NOLEAP
   call error_handler(E_ERR,'get_date',errstring,source,revision,revdate)
end select
end subroutine get_date



subroutine get_date_gregorian(time, year, month, day, hour, minute, second)
!------------------------------------------------------------------------
!
! Computes date corresponding to time for gregorian calendar

implicit none


type(time_type), intent(in)  :: time
integer,         intent(out) :: second, minute, hour, day, month, year

integer :: m,t
integer :: num_days, iyear, days_this_month, days_this_year
logical :: leap

integer, parameter :: base_year= 1601

if ( .not. module_initialized ) call time_manager_init


! Do this the inefficient inelegant way for now, top is 10000 years
num_days = time%days
do iyear = base_year, 10000

! Is this a leap year? Gregorian calendar assigns each year evenly
! divisible by 4 that is not a century year unevenly divisible by 400
! as a leap-year. (i.e. 1700,1800,1900 are not leap-years, 2000 is)
   leap=(modulo(iyear,4).eq.0)
   if((modulo(iyear,100).eq.0).and.(modulo(iyear,400).ne.0))then
      leap=.false.
   endif

   if(leap) then
      days_this_year = 366
   else
      days_this_year = 365
   endif

   if(num_days >= days_this_year) then
      num_days = num_days - days_this_year
   else
      year = iyear
      goto 111
   endif

end do

111 continue


! find month and day
do m = 1, 12
   month = m
   days_this_month = days_per_month(m)
   if(leap .and. m == 2) days_this_month = 29
   if(num_days < days_this_month) exit
   num_days = num_days - days_this_month
end do


day = num_days + 1


! Find hour,minute and second
t      = time%seconds
hour   = t / (60 * 60)
t      = t - hour * (60 * 60)
minute = t / 60
second = t - 60 * minute


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

if ( .not. module_initialized ) call time_manager_init

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

if ( .not. module_initialized ) call time_manager_init

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

if ( .not. module_initialized ) call time_manager_init

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

integer            :: oseconds, ominutes, ohours
character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

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
case(NOLEAP)
   set_date =   set_date_no_leap(year, month, day, ohours, ominutes, oseconds)
case default
   write(errstring,*)'type is ',calendar_type,' must be one of ', &
                      THIRTY_DAY_MONTHS,GREGORIAN,JULIAN,NOLEAP
   call error_handler(E_ERR,'set_date',errstring,source,revision,revdate)
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

character(len=129) :: errstring
integer :: oseconds, ominutes, ohours
integer :: ndays, m, nleapyr
integer :: base_year = 1601
logical :: leap

if ( .not. module_initialized ) call time_manager_init

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  ))   ohours = hours

! Need to check for bogus dates

if(    oseconds > 59 .or. oseconds < 0 .or. &
       ominutes > 59 .or. ominutes < 0 .or. &
       ohours   > 23 .or. ohours   < 0 .or. &
                          day      < 1 .or. &
       month    > 12 .or. month    < 1 .or. &
                          year     < base_year) then
   write(errstring,*)'s,m,h,d,mn,y',oseconds,ominutes,ohours,day,month,year,' not a valid date'
   call error_handler(E_ERR,'set_date_gregorian',errstring,source,revision,revdate)
endif

if(month /= 2 .and. day > days_per_month(month)) then
   write(errstring,*)'month (',month,') does not have ',day,' days.'
   call error_handler(E_ERR,'set_date_gregorian',errstring,source,revision,revdate)
endif

! Is this a leap year? Gregorian calandar assigns each year evenly
! divisible by 4 that is not a century year unevenly divisible by 400
! as a leap-year. (i.e. 1700,1800,1900 are not leap-years, 2000 is)
  leap=(modulo(year,4).eq.0)
  if((modulo(year,100).eq.0).and.(modulo(year,400).ne.0))then
   leap=.false.
  endif

! Finish checking for day specication errors
if(month == 2 .and. (day > 29 .or. ((.not. leap) .and. day > 28))) then
   write(errstring,*)'month (',month,') does not have ',day,' days in a lon-leap year.'
   call error_handler(E_ERR,'set_date_gregorian',errstring,source,revision,revdate)
endif

! compute number of leap years fully past since base_year

nleapyr = (year - base_year) / 4 - (year - base_year) / 100 + (year - base_year) / 400

! Count up days in this year
ndays = 0
do m=1,month-1
 ndays = ndays + days_per_month(m)
 if(leap .and. m == 2) ndays = ndays + 1
enddo

set_date_gregorian%seconds = oseconds + 60*(ominutes + 60 * ohours)
set_date_gregorian%days = day - 1 +  ndays + 365*(year - base_year - nleapyr) + 366*(nleapyr)

end function set_date_gregorian



function set_date_julian(year, month, day, hours, minutes, seconds)
!------------------------------------------------------------------------
!
! Returns time corresponding to date for julian calendar.

implicit none

integer, intent(in)           :: day, month, year
integer, intent(in), optional :: seconds, minutes, hours
type(time_type)               :: set_date_julian

character(len=129) :: errstring
integer :: oseconds, ominutes, ohours
integer :: ndays, m, nleapyr
logical :: leap

if ( .not. module_initialized ) call time_manager_init

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
                       year     < 1) then
   write(errstring,*)'s,m,h,d,mn,y',oseconds,ominutes,ohours,day,month,year,' not a valid date'
   call error_handler(E_ERR,'set_date_julian',errstring,source,revision,revdate)
endif

if(month /= 2 .and. day > days_per_month(month)) then
   write(errstring,*)'month (',month,') does not have ',day,' days.'
   call error_handler(E_ERR,'set_date_julian',errstring,source,revision,revdate)
endif

! Is this a leap year? 
leap = (modulo(year,4) == 0)

! compute number of complete leap years from year 1
nleapyr = (year - 1) / 4

! Finish checking for day specication errors

if(month == 2 .and. (day > 29 .or. ((.not. leap) .and. day > 28))) then
   write(errstring,*)'non-leapyear month (',month,') does not have ',day,' days.'
   call error_handler(E_ERR,'set_date_julian',errstring,source,revision,revdate)
endif

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

character(len=129) :: errstring
integer :: oseconds, ominutes, ohours

if ( .not. module_initialized ) call time_manager_init

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
                       year     < 1 ) then
   write(errstring,*)'s,m,h,d,mn,y',oseconds,ominutes,ohours,day,month,year,' not a valid date'
   call error_handler(E_ERR,'set_date_thirty',errstring,source,revision,revdate)
endif

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

character(len=129) :: errstring
integer :: oseconds, ominutes, ohours
integer :: ndays, m

if ( .not. module_initialized ) call time_manager_init

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
   write(errstring,*)'s,m,h,d,mn,y',oseconds,ominutes,ohours,day,month,year,' not a valid date'
   call error_handler(E_ERR,'set_date_no_leap',errstring,source,revision,revdate)
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
character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

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
case(NOLEAP)
   increment_date = increment_no_leap(time, oyears, omonths, odays, &
                                     ohours, ominutes, oseconds)
case default
   write(errstring,*)'calendar type (',calendar_type,') must be one of ', &
                      THIRTY_DAY_MONTHS ,GREGORIAN,JULIAN,NOLEAP
   call error_handler(E_ERR,'increment_date',errstring,source,revision,revdate)
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
character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

call error_handler(E_ERR,'increment_gregorian','not implemented',source,revision,revdate)

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
    omonths  < 0 .or. oyears   < 0) then
   write(errstring,*)'illegal increment s,m,h,d,mn,y', &
                      oseconds,ominutes,ohours,odays,omonths,oyears
    call error_handler(E_ERR,'increment_gregorian',errstring,source,revision,revdate)
endif

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
character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

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
    omonths  < 0 .or. oyears   < 0 ) then
   write(errstring,*)'illegal increment s,m,h,d,mn,y', &
                     oseconds,ominutes,ohours,odays,omonths,oyears
    call error_handler(E_ERR,'increment_julian',errstring,source,revision,revdate)
endif

!  There are a number of other bad types of increments that should be
!  prohibited here; the addition is not associative
!  Easiest thing is to only let month and year be incremented by themselves
!  This is slight overkill since year is not really a problem.

! QJLA .... see above comment ... will year be a problem?

if( omonths /= 0 .and. &
    ( oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
      odays    /= 0 .or. oyears   /= 0)  ) then
   call error_handler(E_ERR,'increment_julian', &
      'month must not be incremented with other units',source,revision,revdate)
endif

if( oyears /= 0 .and. &
    (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
     odays    /= 0 .or. omonths  /= 0) ) then
   call error_handler(E_ERR,'increment_julian', &
      'year must not be incremented with other units',source,revision,revdate)
endif

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

character(len=129) :: errstring
integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cday

if ( .not. module_initialized ) call time_manager_init

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
    omonths  < 0 .or. oyears   < 0 ) then
   write(errstring,*)'illegal increment s,m,h,d,mn,y', &
                     oseconds,ominutes,ohours,odays,omonths,oyears
    call error_handler(E_ERR,'increment_thirty',errstring,source,revision,revdate)
endif

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
character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

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
    omonths  < 0 .or. oyears   < 0) then
    write(errstring,*)'illegal increment s,m,h,d,mn,y', &
                     oseconds,ominutes,ohours,odays,omonths,oyears
    call error_handler(E_ERR,'increment_no_leap',errstring,source,revision,revdate)
endif

!  There are a number of other bad types of increments that should be
!  prohibited here; the addition is not associative
!  Easiest thing is to only let month and year be incremented by themselves
!  This is slight overkill since year is not really a problem.

if(omonths /= 0 .and. &
   (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
    odays    /= 0 .or. oyears   /= 0)) then
   call error_handler(E_ERR,'increment_no_leap', &
      'month must not be incremented with other units',source,revision,revdate)
endif

if(oyears /= 0 .and. &
   (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
    odays    /= 0 .or. omonths  /= 0)) then
   call error_handler(E_ERR,'increment_no_leap', &
      'year must not be incremented with other units',source,revision,revdate)
endif

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
character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

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
case(NOLEAP)
   decrement_date = decrement_no_leap(time, oyears, omonths, odays, &
                       ohours, ominutes, oseconds)
case default
   write(errstring,*)'calendar type (',calendar_type,') not allowed.'
   call error_handler(E_ERR,'decrement_date',errstring,source,revision,revdate)
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

character(len=129) :: errstring
integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cminute, chour, cday, cmonth, cyear

if ( .not. module_initialized ) call time_manager_init

call error_handler(E_ERR,'decrement_gregorian','not implemented',source,revision,revdate)

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
    omonths  < 0 .or. oyears < 0) then
    write(errstring,*)'illegal decrement ',oseconds,ominutes,ohours,odays,omonths,oyears
    call error_handler(E_ERR,'decrement_gregorian',errstring,source,revision,revdate)
endif

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
character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours
odays    = 0; if(present(days   )) odays    = days
omonths  = 0; if(present(months )) omonths  = months
oyears   = 0; if(present(years  )) oyears   = years

! Increment must be positive definite

if( oseconds < 0 .or. ominutes < 0 .or. ohours < 0 .or. &
    odays    < 0 .or. omonths  < 0 .or. oyears < 0) then
    write(errstring,*)'illegal decrement ',oseconds,ominutes,ohours,odays,omonths,oyears
    call error_handler(E_ERR,'decrement_julian',errstring,source,revision,revdate)
endif

!  There are a number of other bad types of decrements that should be
!  prohibited here; the subtraction is not associative
!  Easiest thing is to only let month and year be decremented by themselves
!  This is slight overkill since year is not really a problem.

if(omonths /= 0 .and. (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
   odays  /= 0 .or. oyears /= 0)) then
   call error_handler(E_ERR,'decrement_julian', &
      'month must not be decremented with other units',source,revision,revdate)
endif

if(oyears /= 0 .and. (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
   odays  /= 0 .or. omonths /= 0)) then
   call error_handler(E_ERR,'decrement_julian', &
      'year must not be decremented with other units',source,revision,revdate)
endif

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
if(cyear < 1) then
   write(errstring,*)'Illegal date results ... ',cyear
   call error_handler(E_ERR,'decrement_julian',errstring,source,revision,revdate)
endif

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

character(len=129) :: errstring
integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cday

if ( .not. module_initialized ) call time_manager_init

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours
odays    = 0; if(present(days   )) odays    = days
omonths  = 0; if(present(months )) omonths  = months
oyears   = 0; if(present(years  )) oyears   = years

! Increment must be positive definite

if( oseconds < 0 .or. ominutes < 0 .or. ohours < 0 .or. &
    odays    < 0 .or. omonths  < 0 .or. oyears < 0) then
    write(errstring,*)'illegal decrement ',oseconds,ominutes,ohours,odays,omonths,oyears
    call error_handler(E_ERR,'decrement_thirty',errstring,source,revision,revdate)
endif

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
character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  )) ohours   = hours
odays    = 0; if(present(days   )) odays    = days
omonths  = 0; if(present(months )) omonths  = months
oyears   = 0; if(present(years  )) oyears   = years

! Increment must be positive definite

if( oseconds < 0 .or. ominutes < 0 .or. ohours < 0 .or. &
    odays    < 0 .or. omonths  < 0 .or. oyears < 0 ) then
    write(errstring,*)'illegal decrement ',oseconds,ominutes,ohours,odays,omonths,oyears
    call error_handler(E_ERR,'decrement_no_leap',errstring,source,revision,revdate)
endif

!  There are a number of other bad types of decrements that should be
!  prohibited here; the subtraction is not associative
!  Easiest thing is to only let month and year be decremented by themselves
!  This is slight overkill since year is not really a problem.

if(omonths /= 0 .and. (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
   odays   /= 0 .or. oyears /= 0)) then
   call error_handler(E_ERR,'decrement_no_leap', &
      'month must not be decremented with other units',source,revision,revdate)
endif

if(oyears /= 0 .and. (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
   odays  /= 0 .or. omonths /= 0)) then
   call error_handler(E_ERR,'decrement_no_leap', &
      'year must not be decremented with other units',source,revision,revdate)
endif

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

if(cyear < 1) then
   write(errstring,*)'Illegal date results ... ',cyear
   call error_handler(E_ERR,'decrement_no_leap',errstring,source,revision,revdate)
endif

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

character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   days_in_month = days_in_month_thirty(time)
case(GREGORIAN)
   days_in_month = days_in_month_gregorian(time)
case(JULIAN)
   days_in_month = days_in_month_julian(time)
case(NOLEAP)
   days_in_month = days_in_month_no_leap(time)
case default
   write(errstring,*)'Invalid calendar type (',calendar_type,')'
   call error_handler(E_ERR,'days_in_month',errstring,source,revision,revdate)
end select
end function days_in_month



function days_in_month_gregorian(time)
!--------------------------------------------------------------------------
!
! Returns the number of days in a gregorian month.

implicit none

type(time_type), intent(in) :: time
integer                     :: days_in_month_gregorian

if ( .not. module_initialized ) call time_manager_init

call error_handler(E_ERR,'days_in_month_gregorian','not implemented',source,revision,revdate)
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

if ( .not. module_initialized ) call time_manager_init

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

if ( .not. module_initialized ) call time_manager_init

days_in_month_thirty = 30

end function days_in_month_thirty



function days_in_month_no_leap(time)
!--------------------------------------------------------------------------
! Returns the number of days in a 365 day year month.

implicit none

type(time_type), intent(in) :: time
integer                     :: days_in_month_no_leap

integer :: seconds, minutes, hours, day, month, year

if ( .not. module_initialized ) call time_manager_init

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

character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   leap_year = leap_year_thirty(time)
case(GREGORIAN)
   leap_year = leap_year_gregorian(time)
case(JULIAN)
   leap_year = leap_year_julian(time)
case(NOLEAP)
   leap_year = leap_year_no_leap(time)
case default
   write(errstring,*)'invalid calendar type (',calendar_type,')'
   call error_handler(E_ERR,'leap_year',errstring,source,revision,revdate)
end select
end function leap_year



function leap_year_gregorian(time)
!--------------------------------------------------------------------------
!
! Is this a leap year for gregorian calendar?

implicit none

type(time_type), intent(in) :: time
logical                     :: leap_year_gregorian

if ( .not. module_initialized ) call time_manager_init

call error_handler(E_ERR,'leap_year_gregorian','not implemented',source,revision,revdate)

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

if ( .not. module_initialized ) call time_manager_init

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

if ( .not. module_initialized ) call time_manager_init

leap_year_thirty = .FALSE.

end function leap_year_thirty



function leap_year_no_leap(time)
!--------------------------------------------------------------------------
!
! Another tough one; no leap year returns false for leap year inquiry.

implicit none

type(time_type), intent(in) :: time
logical                     :: leap_year_no_leap

if ( .not. module_initialized ) call time_manager_init

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

character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   length_of_year = length_of_year_thirty()
case(GREGORIAN)
   length_of_year = length_of_year_gregorian()
case(JULIAN)
   length_of_year = length_of_year_julian()
case(NOLEAP)
   length_of_year = length_of_year_no_leap()
case default
   write(errstring,*)'invalid calendar type (',calendar_type,')'
   call error_handler(E_ERR,'length_of_year',errstring,source,revision,revdate)
end select
end function length_of_year



function length_of_year_thirty()
!--------------------------------------------------------------------------
!

implicit none

type(time_type) :: length_of_year_thirty

if ( .not. module_initialized ) call time_manager_init

length_of_year_thirty = set_time(0, 360)

end function length_of_year_thirty



function length_of_year_gregorian()
!---------------------------------------------------------------------------

implicit none

type(time_type) :: length_of_year_gregorian

if ( .not. module_initialized ) call time_manager_init

length_of_year_gregorian = set_time(0, 0)

call error_handler(E_ERR,'length_of_year_gregorian', &
       'not implemented',source,revision,revdate)

end function length_of_year_gregorian



function length_of_year_julian()
!--------------------------------------------------------------------------

implicit none

type(time_type) :: length_of_year_julian

if ( .not. module_initialized ) call time_manager_init

length_of_year_julian = set_time((24 / 4) * 60 * 60, 365)

end function length_of_year_julian



function length_of_year_no_leap()
!--------------------------------------------------------------------------

implicit none

type(time_type) :: length_of_year_no_leap

if ( .not. module_initialized ) call time_manager_init

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

character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

select case(calendar_type)
case(THIRTY_DAY_MONTHS)
   days_in_year = days_in_year_thirty(time)
case(GREGORIAN)
   days_in_year = days_in_year_gregorian(time)
case(JULIAN)
   days_in_year = days_in_year_julian(time)
case(NOLEAP)
   days_in_year = days_in_year_no_leap(time)
case default
   write(errstring,*)'invalid calendar type (',calendar_type,')'
   call error_handler(E_ERR,'days_in_year',errstring,source,revision,revdate)
end select
end function days_in_year



function days_in_year_thirty(time)
!--------------------------------------------------------------------------

implicit none

type(time_type), intent(in) :: time
integer                     :: days_in_year_thirty

if ( .not. module_initialized ) call time_manager_init

days_in_year_thirty = 360

end function days_in_year_thirty



function days_in_year_gregorian(time)
!---------------------------------------------------------------------------

implicit none

type(time_type), intent(in) :: time
integer                     :: days_in_year_gregorian

if ( .not. module_initialized ) call time_manager_init

days_in_year_gregorian = 0

call error_handler(E_ERR,'days_in_year_gregorian', &
      'not implemented',source,revision,revdate)

end function days_in_year_gregorian



function days_in_year_julian(time)
!--------------------------------------------------------------------------

implicit none

type(time_type), intent(in) :: time
integer                     :: days_in_year_julian

if ( .not. module_initialized ) call time_manager_init

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

if ( .not. module_initialized ) call time_manager_init

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

character(len=129) :: errstring

if ( .not. module_initialized ) call time_manager_init

if( n < 1 .or. n > 12 ) then
   write(errstring,*)'Illegal month (',n,')'
   call error_handler(E_ERR,'month_name',errstring,source,revision,revdate)
endif

month_name = months(n)

end function month_name

!==========================================================================


function julian_day(year, month, day)

! Given a date, computes the day from the beginning of the year.

implicit none

integer, intent(in) :: day, month, year
integer             :: julian_day

integer             :: m
logical             :: leap

if ( .not. module_initialized ) call time_manager_init

if (calendar_type /= GREGORIAN) then
   call error_handler(E_ERR,'julian_day', &
        'only implemented for GREGORIAN calendar',source,revision,revdate)
endif

julian_day = 0
leap = (modulo(year,4) == 0)
if((modulo(year,100).eq.0).and.(modulo(year,400).ne.0))then
   leap=.false.
endif

do m = 1, month - 1
   julian_day = julian_day + days_per_month(m)
   if(leap .and. m == 2) julian_day = julian_day + 1
enddo
julian_day = julian_day + day

end function julian_day


subroutine time_manager_init ( )
!------------------------------------------------------------------------
!
! initialization routine
! this routine should be called, even though all it does is write
! the version information to the log file

   if ( module_initialized ) return  ! silent return if already called

   call register_module (source, revision, revdate)
   module_initialized  = .true.

end subroutine time_manager_init



subroutine print_time (time,str,iunit)
!------------------------------------------------------------------------

type(time_type)  , intent(in)           :: time
character (len=*), intent(in), optional :: str
integer          , intent(in), optional :: iunit

integer           :: s,d, ns,nd, unit_in
character(len=13) :: fmt

if ( .not. module_initialized ) call time_manager_init

! prints the time to standard output (or optional iunit) as days and seconds
! NOTE: there is no check for PE number

! NEED TO GET DEFAULT FOR STANDARD OUT< HARD CODED FOR NOW
!!!  unit_in = stdout()
  unit_in = 6

  if (present(iunit)) unit_in = iunit

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



subroutine print_date (time,str,iunit)
!------------------------------------------------------------------------

type(time_type)  , intent(in)           :: time
character (len=*), intent(in), optional :: str
integer          , intent(in), optional :: iunit

integer          :: y,mo,d,h,m,s, unit_in
character(len=9) :: mon

if ( .not. module_initialized ) call time_manager_init

! prints the time to standard output (or optional iunit) as a date
! NOTE: there is no check for PE number

!!!  unit_in = stdout()
  unit_in = 6

  if (present(iunit)) unit_in = iunit

  call get_date (time,y,mo,d,h,m,s)
  mon = month_name(mo)
  if (present(str)) then
     write (unit_in,10) trim(str)//' ', y,mon(1:3),' ',d,' ',h,':',m,':',s
  else
     write (unit_in,10)       'DATE: ', y,mon(1:3),' ',d,' ',h,':',m,':',s
  endif
10 format (a,i4,1x,a3,4(a1,i2.2))

end subroutine print_date



function read_time(file_unit, form, ios_out)
!--------------------------------------------------------------------------------
!

implicit none

integer,          intent(in)            :: file_unit
character(len=*), intent(in), optional  :: form
integer,          intent(out), optional :: ios_out

type(time_type)   :: read_time
integer           :: secs, days, ios
character(len=32) :: fileformat

character(len=128) :: str1

if ( .not. module_initialized ) call time_manager_init

fileformat = "ascii"   ! supply default
if (present(form)) fileformat = trim(adjustl(form))

SELECT CASE (fileformat)
   CASE ("unf","UNF","unformatted","UNFORMATTED")
      read(file_unit, iostat=ios) secs, days
   CASE DEFAULT
      read(file_unit, *, iostat=ios) secs, days
END SELECT

if ( ios /= 0 ) then

   ! If ios_out argument is present, just return a non-zero ios
   if(present(ios_out)) then
      ios_out = ios
      return
   endif

   ! Otherwise, read error is fatal, print message and stop
   call dump_unit_attributes(file_unit)   ! TJH DEBUG statement
   write(str1,*)'read status is',ios
   call error_handler(E_ERR,'read_time',str1,source,revision,revdate)
else
   read_time = set_time(secs, days)
   if(present(ios_out)) ios_out = 0
endif

end function read_time



subroutine write_time(file_unit, time, form)
!------------------------------------------------------------------------
! The time is expected to be written as either 
! an unformatted binary (form == "unformatted")
! or free-format ascii (form /= "unformatted")

implicit none

integer,          intent(in)           :: file_unit
type(time_type),  intent(in)           :: time
character(len=*), intent(in), optional :: form

integer           :: secs, days
character(len=32) :: fileformat

if ( .not. module_initialized ) call time_manager_init

fileformat = "ascii"   ! supply default
if (present(form)) fileformat = trim(adjustl(form))

call get_time(time, secs, days)

SELECT CASE (fileformat)
   CASE ("unf","UNF","unformatted","UNFORMATTED")
      write(file_unit) secs, days
   CASE DEFAULT
      write(file_unit,'(i6,1x,i10)') secs, days
END SELECT

end subroutine write_time


subroutine interactive_time(time)
!------------------------------------------------------------------------

implicit none

type(time_type), intent(inout) :: time
integer                        :: second, minute, hour, day, month, year

if ( .not. module_initialized ) call time_manager_init

if (calendar_type == GREGORIAN) then
   write(*, *) 'input date (as integers): year month day hour minute second'
   read(*, *) year, month, day, hour, minute, second
   time = set_date(year, month, day, hour, minute, second)
else
   write(*, *) 'input time in days and seconds (as integers)'
   read(*, *) day, second
   time = set_time(second, day)
endif

end subroutine interactive_time


end module time_manager_mod
