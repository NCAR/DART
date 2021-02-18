! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module time_manager_mod

use     types_mod, only : missing_i, digits12, i8
use utilities_mod, only : error_handler, E_DBG, E_MSG, E_WARN, E_ERR, &
                          dump_unit_attributes, to_upper, &
                          ascii_file_format

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
public :: interval_alarm, repeat_alarm, generate_seed, time_index_sort

! List of available calendar types
public :: THIRTY_DAY_MONTHS,    JULIAN,    GREGORIAN,  NOLEAP,   NO_CALENDAR, &
          GREGORIAN_MARS,   SOLAR_MARS

! Subroutines and functions involving relations between time and calendar
public :: set_calendar_type, get_calendar_type, get_calendar_string
public :: set_date,       set_date_gregorian,         set_date_julian, &
                          set_date_thirty,            set_date_no_leap, &
                          set_date_gregorian_mars, set_date_solar_mars
public :: get_date,       get_date_gregorian,         get_date_julian, &
                          get_date_thirty,            get_date_no_leap, &
                          get_date_gregorian_mars, get_date_solar_mars
public :: increment_date, increment_gregorian,        increment_julian, &
                          increment_thirty,           increment_no_leap, &
                          increment_gregorian_mars, increment_solar_mars
public :: decrement_date, decrement_gregorian,        decrement_julian, &
                          decrement_thirty,           decrement_no_leap, &
                          decrement_gregorian_mars, decrement_solar_mars
public :: days_in_month,  days_in_month_gregorian,    days_in_month_julian, &
                          days_in_month_no_leap,      days_in_month_thirty, &
                          days_in_month_gregorian_mars, days_in_month_solar_mars
public :: leap_year,      leap_year_gregorian,        leap_year_julian, &
                          leap_year_no_leap,          leap_year_thirty, &
                          leap_year_gregorian_mars, leap_year_solar_mars
public :: length_of_year, length_of_year_thirty,      length_of_year_julian, &
                          length_of_year_gregorian,   length_of_year_no_leap, &
                          length_of_year_gregorian_mars, length_of_year_solar_mars
public :: days_in_year,   days_in_year_thirty,        days_in_year_julian, &
                          days_in_year_gregorian,     days_in_year_no_leap, &
                          days_in_year_gregorian_mars, days_in_year_solar_mars

public :: month_name

public :: julian_day

! Subroutines and functions for basic I/O
public :: time_manager_init, print_time, print_date
public :: write_time, read_time, interactive_time

character(len=*), parameter :: source = 'time_manager_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

! Global data to define calendar type
integer, parameter :: THIRTY_DAY_MONTHS = 1,      JULIAN = 2, &
                      GREGORIAN = 3,              NOLEAP = 4, &
                      NO_CALENDAR = 0,            GREGORIAN_MARS = 5, &
                      SOLAR_MARS = 6
! FIXME: should be a namelist to select default calendar.  make no calendar the
! default for now
integer, private :: calendar_type = NO_CALENDAR, max_type = 6

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

interface set_calendar_type
   module procedure set_calendar_type_integer
   module procedure set_calendar_type_string
end interface

!======================================================================

logical, save :: module_initialized = .false.

character(len=129) :: errstring

!======================================================================

contains

! First define all operations on time intervals independent of calendar



function set_time(seconds, days)
!---------------------------------------------------------------------------
!
! Returns a time interval corresponding to this number of days and seconds.
! The arguments must not be negative but are otherwise unrestricted.

integer, intent(in)           :: seconds
integer, intent(in), optional :: days
type(time_type)               :: set_time

integer            :: days_in

if ( .not. module_initialized ) call time_manager_init

days_in = 0;  if (present(days)) days_in = days

! Negative time offset is illegal

write(errstring,*)'seconds, days are ',seconds, days_in,' cannot be negative'
if(seconds < 0 .or. days_in < 0) &
   call error_handler(E_ERR,'set_time',errstring,source)

! Make sure seconds greater than a day are fixed up.
! Extra parens to force the divide before the multiply are REQUIRED 
! on some compilers to prevent them from combining constants first 
! which makes the expression always return 0.

set_time%seconds = seconds - (seconds / (60*60*24)) * (60*60*24)

! Check for overflow on days before doing operation

write(errstring,*)'seconds is ',seconds,' overflowing conversion to days'
if(seconds / (60*60*24)  >= huge(days_in) - days_in) &
   call error_handler(E_ERR,'set_time',errstring,source)

set_time%days = days_in + seconds / (60*60*24)

end function set_time



function set_time_missing()
!---------------------------------------------------------------------------
!
! Initialize a time derived type to a known value.  It is impossible for
! this to be set by the set_time() function since the missing integers
! are negative values.

type(time_type) :: set_time_missing

set_time_missing%seconds = missing_i
set_time_missing%days    = missing_i

end function set_time_missing



subroutine get_time(time, seconds, days)
!---------------------------------------------------------------------------
!
! Returns days and seconds ( < 86400 ) corresponding to a time.
! If the optional 'days' argument is not given, the days are converted
! to seconds and the total time is returned as seconds.

type(time_type), intent(in)            :: time
integer,         intent(out)           :: seconds
integer,         intent(out), optional :: days

if ( .not. module_initialized ) call time_manager_init

seconds = time%seconds
if (present(days)) then
  days = time%days
else

  write(errstring,*)'seconds is ',seconds,' overflowing conversion to days'

  if (time%days > (huge(seconds) - seconds)/(60*60*24)) &
     call error_handler(E_ERR,'get_time',errstring,source)

  seconds = seconds + time%days * (60*60*24)

endif

end subroutine get_time



function increment_time(time, seconds, days)
!-------------------------------------------------------------------------
!
! Increments a time by seconds and days; increments cannot be negative.

type(time_type), intent(in)           :: time
integer,         intent(in)           :: seconds
integer,         intent(in), optional :: days
type(time_type)                       :: increment_time

integer            :: days_in

if ( .not. module_initialized ) call time_manager_init

days_in = 0;  if (present(days)) days_in = days

! Increment must be positive definite

if(seconds < 0 .or. days_in < 0) then
   write(errstring,*)'Negative increment (',seconds, days,') (seconds,days) not allowed'
   call error_handler(E_ERR,'increment_time',errstring,source)
endif

! Watch for immediate overflow on days or seconds

if(days_in >= huge(days_in) - time%days) then
   write(errstring,*)'integer overflow (',days_in, time%days,') in days'
   call error_handler(E_ERR,'increment_time',errstring,source)
endif

if(seconds >= huge(seconds) - time%seconds) then
   write(errstring,*)'integer overflow (',seconds, time%seconds,') in seconds'
   call error_handler(E_ERR,'increment_time',errstring,source)
endif

increment_time = set_time(time%seconds + seconds, time%days + days_in)

end function increment_time



function decrement_time(time, seconds, days)
!--------------------------------------------------------------------------
!
! Decrements a time by seconds and days; decrements cannot be negative.

type(time_type), intent(in)           :: time
integer,         intent(in)           :: seconds
integer,         intent(in), optional :: days
type(time_type)                       :: decrement_time

integer            :: cseconds, cdays

if ( .not. module_initialized ) call time_manager_init

cdays = 0;  if (present(days)) cdays = days

! Decrement must be positive definite

if(seconds < 0 .or. cdays < 0) then
   write(errstring,*)'Negative decrement (',seconds,cdays,') (seconds,days) not allowed.'
   call error_handler(E_ERR,'decrement_time',errstring,source)
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
   call error_handler(E_ERR,'decrement_time',errstring,source)
endif

decrement_time%seconds = cseconds
decrement_time%days    = cdays

end function decrement_time



function time_gt(time1, time2)
!--------------------------------------------------------------------------
!
! Returns true if time1 > time2

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

type(time_type), intent(in) :: time1, time2
logical                     :: time_ge

if ( .not. module_initialized ) call time_manager_init

time_ge = (time_gt(time1, time2) .or. time_eq(time1, time2))

end function time_ge



function time_lt(time1, time2)
!--------------------------------------------------------------------------
!
! Returns true if time1 < time2

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

type(time_type), intent(in) :: time1, time2
logical                     :: time_le

if ( .not. module_initialized ) call time_manager_init

time_le = (time_lt(time1, time2) .or. time_eq(time1, time2))

end function time_le



function time_eq(time1, time2)
!--------------------------------------------------------------------------
!
! Returns true if time1 == time2

type(time_type), intent(in) :: time1, time2
logical                     :: time_eq

if ( .not. module_initialized ) call time_manager_init

time_eq = (time1%seconds == time2%seconds .and. time1%days == time2%days)

end function time_eq



function time_ne(time1, time2)
!--------------------------------------------------------------------------
!
! Returns true if time1 /= time2

type(time_type), intent(in) :: time1, time2
logical                     :: time_ne

if ( .not. module_initialized ) call time_manager_init

time_ne = (.not. time_eq(time1, time2))

end function time_ne



function time_plus(time1, time2)
!-------------------------------------------------------------------------
!
! Returns sum of two time_types

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

type(time_type), intent(in) :: time
integer,         intent(in) :: n
type(time_type)             :: time_scalar_mult

integer            :: days, seconds
real(digits12)     :: sec_prod

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
      call error_handler(E_ERR,'time_scalar_mult',errstring,source) 
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

type(time_type), intent(in) :: time1, time2
integer                     :: time_divide

real(digits12)     :: d1, d2

if ( .not. module_initialized ) call time_manager_init

! Convert time intervals to floating point days; risky for general performance?

d1 = time1%days * (60.0_digits12 * 60.0_digits12 * 24.0_digits12) + (time1%seconds*1.0_digits12)
d2 = time2%days * (60.0_digits12 * 60.0_digits12 * 24.0_digits12) + (time2%seconds*1.0_digits12) 

! Get integer quotient of this, check carefully to avoid round-off problems.

time_divide = floor(d1 / d2)

! Verify time_divide*time2 is <= time1 and (time_divide + 1)*time2 is > time1

if(time_divide * time2 > time1 .or. (time_divide + 1) * time2 <= time1) then
   write(errstring,*)'quotient error time1,time2 = ',time1, time2
   call error_handler(E_ERR,'time_divide',errstring,source)
endif

end function time_divide



function time_real_divide(time1, time2)
!-------------------------------------------------------------------------
!
! Returns the double precision quotient of two times

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

integer,         intent(in) :: n
type(time_type), intent(in) :: time
type(time_type)             :: time_scalar_divide

real(digits12)  :: d, div
integer         :: days, seconds
type(time_type) :: prod1, prod2

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
   call error_handler(E_ERR,'time_scalar_divide',errstring,source)
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

subroutine set_calendar_type_integer(mytype)

! Selects calendar for default mapping from time to date - if you know 
! the magic integer for the calendar of interest. 

integer, intent(in) :: mytype

if ( .not. module_initialized ) call time_manager_init

if(mytype <  0 .or. mytype > max_type) then
   write(errstring,*)'Illegal calendar type ',mytype,' must be >= 0 or <= ',max_type
   call error_handler(E_ERR,'set_calendar_type_integer',errstring,source)
endif
calendar_type = mytype

end subroutine set_calendar_type_integer



subroutine set_calendar_type_string(calstring)

! Selects calendar for default mapping from time to date - given a string. 

character(len=*), intent(in) :: calstring

integer, parameter :: max_calendar_string_length = len_trim('THIRTY_DAY_MONTHS')

character(len=len(calstring))             :: str1
character(len=max_calendar_string_length) :: cstring
logical :: found_calendar = .false.
integer :: i

if ( .not. module_initialized ) call time_manager_init

str1 = adjustl(calstring)

if ( len_trim(str1) > max_calendar_string_length ) then
   write(errstring,*)'Illegal calendar ',trim(calstring)
   call error_handler(E_ERR,'set_calendar_type_string',errstring,source)
endif

cstring = trim(str1)
call to_upper(cstring)

! Using 'cstring' as the substring (2nd argument), we remove 
! the ambiguity of someone trying to use a calendar string
! of 'no' ... which could either match no_calendar or noleap.
! We must check for the gregorian_mars calendar before
! the gregorian calendar for similar reasons.

WhichCalendar : do i = 0, max_type 

   if     ( cstring == 'NO_CALENDAR' ) then
           calendar_type  = NO_CALENDAR
           found_calendar = .true.
           exit WhichCalendar
   elseif ( cstring == 'NO CALENDAR' ) then   ! allow this as a synonym 
           calendar_type  = NO_CALENDAR
           found_calendar = .true.
           exit WhichCalendar
   elseif ( cstring == 'NONE' ) then          ! also allow this
           calendar_type  = NO_CALENDAR
           found_calendar = .true.
           exit WhichCalendar
   elseif ( cstring == 'THIRTY_DAY_MONTHS' ) then
           calendar_type  = THIRTY_DAY_MONTHS
           found_calendar = .true.
           exit WhichCalendar
   elseif ( cstring == 'JULIAN' ) then
           calendar_type  = JULIAN
           found_calendar = .true.
           exit WhichCalendar
   elseif ( cstring == 'NOLEAP' ) then
           calendar_type  = NOLEAP
           found_calendar = .true.
           exit WhichCalendar
   elseif ( cstring == 'GREGORIAN_MARS' ) then
           calendar_type  = GREGORIAN_MARS
           found_calendar = .true.
           exit WhichCalendar
   elseif ( cstring == 'SOLAR_MARS' ) then
           calendar_type  = SOLAR_MARS
           found_calendar = .true.
           exit WhichCalendar
   elseif ( cstring == 'GREGORIAN' ) then
           calendar_type  = GREGORIAN
           found_calendar = .true.
           exit WhichCalendar
   endif

enddo WhichCalendar

if( .not. found_calendar ) then
   write(errstring,*)'Unknown calendar ',calstring
   call error_handler(E_ERR,'set_calendar_type_string',errstring,source)
endif

end subroutine set_calendar_type_string



function get_calendar_type()
!------------------------------------------------------------------------
!
! Returns default calendar type for mapping from time to date.

integer :: get_calendar_type

if ( .not. module_initialized ) call time_manager_init

get_calendar_type = calendar_type

end function get_calendar_type



subroutine get_calendar_string(mystring)
!------------------------------------------------------------------------
!
! Returns default calendar type for mapping from time to date.

character(len=*), intent(OUT) :: mystring

integer :: i

if ( .not. module_initialized ) call time_manager_init

mystring = '  '

do i = 0,max_type
   if (calendar_type ==            JULIAN) mystring = 'JULIAN'
   if (calendar_type ==            NOLEAP) mystring = 'NOLEAP'
   if (calendar_type ==         GREGORIAN) mystring = 'GREGORIAN'
   if (calendar_type ==       NO_CALENDAR) mystring = 'NO_CALENDAR'
   if (calendar_type ==    GREGORIAN_MARS) mystring = 'GREGORIAN_MARS'
   if (calendar_type ==        SOLAR_MARS) mystring = 'SOLAR_MARS'
   if (calendar_type == THIRTY_DAY_MONTHS) mystring = 'THIRTY_DAY_MONTHS'
enddo

if (len_trim(mystring) < 3) then
   write(errstring,*)'unknown calendar type ', calendar_type
   call error_handler(E_ERR,'get_calendar_string',errstring,source)
endif

end subroutine get_calendar_string



!========================================================================
! START OF get_date BLOCK
!========================================================================

subroutine get_date(time, year, month, day, hour, minute, second)
!------------------------------------------------------------------------
!
! Given a time, computes the corresponding date given the selected calendar

type(time_type), intent(in)  :: time
integer,         intent(out) :: second, minute, hour, day, month, year

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
case(GREGORIAN_MARS)
   ! NOTE: "month" has no meaning in Mars calendar
   call get_date_gregorian_mars(time, year, month, day, hour, minute, second)
case(SOLAR_MARS)
   ! NOTE: "month" has no meaning in Mars calendar
   call get_date_solar_mars(time, year, month, day, hour, minute, second)
case default
   write(errstring,*)'type is ',calendar_type,' must be one of ', &
                      THIRTY_DAY_MONTHS,GREGORIAN,JULIAN,NOLEAP,GREGORIAN_MARS
   call error_handler(E_ERR,'get_date',errstring,source)
end select
end subroutine get_date



subroutine get_date_gregorian(time, year, month, day, hour, minute, second)
!------------------------------------------------------------------------
!
! Computes date corresponding to time for gregorian calendar

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
   leap=greg_leap(iyear)

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


subroutine get_date_gregorian_mars(time, year, month, day, hour, minute, second)
!------------------------------------------------------------------------
!
! Computes date corresponding to time for gregorian MARS calendar
! wgl - we need this to easily map from MarsWRF calendar to time_type
!         NOTE :: "month" has no meaning for the Mars calendar
!         ALSO :: the first day of the year is day 1, not 0
!         PlanetWRF wrfout files follow :: YYYY-DDDDD_HH:MM:SS
!
! To be clear, the first Mars year is:
!   According to MarsWRF         :: 0001-00001_00:00:00  thru  0001-00669_23:59:59
!   According to DART time_type  :: days=0, secs=0       thru  days=668, secs=86399

type(time_type), intent(in)  :: time
integer,         intent(out) :: second, minute, hour, day, month, year

integer :: t
integer :: num_days, iyear, days_this_year

! "base_year" for Mars will be defined as 1 (earliest wrfout file has year = 1)
integer, parameter :: base_year= 1

if ( .not. module_initialized ) call time_manager_init

! Do this the inefficient inelegant way for now, top is 10000 years
num_days = time%days
do iyear = base_year, 10000

   ! No leap years currently defined on Mars -- we generally run with the 
   !   a WRF constant "PLANET_YEAR = 669", even though Mars year really has
   !   668.6 days in its year (where "day" here is really "sol" = 1 solar day)
   days_this_year = 669

   if (num_days >= days_this_year) then
      num_days = num_days - days_this_year
   else
      year = iyear
      goto 111
   endif

end do

111 continue


! No need to deal with months on Mars -- make it 1 by default so that other
!    functions do not break (like print_date and month_name)
month = 1

! Need to add 1 because Mars has NO day 0
day = num_days + 1


! Find hour,minute and second
t      = time%seconds
hour   = t / (60 * 60)
t      = t - hour * (60 * 60)
minute = t / 60
second = t - 60 * minute

end subroutine get_date_gregorian_mars


subroutine get_date_solar_mars(time, year, month, day, hour, minute, second)
!------------------------------------------------------------------------
!
! Computes date corresponding to time for solar MARS calendar
! wgl - we need this to easily map from MarsWRF calendar to time_type
!         NOTE :: "month" has no meaning for the Mars calendar
!         ALSO :: the first day of the year is day 1, not 0
!         PlanetWRF wrfout files follow :: YYYY-DDDDD_HH:MM:SS
!
! CL New method uses Mars Solar Date with no concept of 'year' in the date string
! The year still starts on day 1
! year is always 1

type(time_type), intent(in)  :: time
integer,         intent(out) :: second, minute, hour, day, month, year

integer :: t
integer :: num_days

! "base_year" for Mars will be defined as 1 (earliest wrfout file has year = 1)
integer, parameter :: base_year= 1

if ( .not. module_initialized ) call time_manager_init

year = 1
num_days = time%days

! No need to deal with months on Mars -- make it 1 by default so that other
!    functions do not break (like print_date and month_name)
month = 1

! Need to add 1 because Mars has NO day 0
day = num_days + 1

!cl Time here corresponds to Martian Local Mean Solar Time (LMST), which is what the model 
!cl uses without telling anybody. Essentially, we assume that the Sun 'moves' at a constant
!cl rate across the sky, completing a circle in 24 Mars hours, and we arbitrarily split the 
!cl Mars hour into 3600 Mars seconds, which correspons to roughly 1.02 SI seconds
!cl Inside WRF, Mars delta seconds are converted to SI delta seconds, and the Mars LMST
!cl is used to calculate the position of the sun either as (LMST)*360/24 or (LTST)*360/24. where 
!cl LTST is Local True Solar Time, accounting for variable Solar day length
!cl with orbit (Equation Of Time).
! Find hour,minute and second

t      = time%seconds
hour   = t / (60 * 60)
t      = t - hour * (60 * 60)
minute = t / 60
second = t - 60 * minute

end subroutine get_date_solar_mars


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

integer, intent(in)           :: day, month, year
integer, intent(in), optional :: seconds, minutes, hours
type(time_type)               :: set_date

integer            :: oseconds, ominutes, ohours

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
case(GREGORIAN_MARS)
   set_date =   set_date_gregorian_mars(year, month, day, ohours, ominutes, oseconds)
case(SOLAR_MARS)
   set_date =   set_date_solar_mars(year, month, day, ohours, ominutes, oseconds)
case default
   write(errstring,*)'type is ',calendar_type,' must be one of ', &
                      THIRTY_DAY_MONTHS,GREGORIAN,JULIAN,NOLEAP,GREGORIAN_MARS,SOLAR_MARS
   call error_handler(E_ERR,'set_date',errstring,source)
end select
end function set_date



function set_date_gregorian(year, month, day, hours, minutes, seconds)
!------------------------------------------------------------------------
!
! Computes time corresponding to date for gregorian calendar.

integer, intent(in)           :: year, month, day 
integer, intent(in), optional :: hours, minutes, seconds
type(time_type)               :: set_date_gregorian

integer :: totseconds, totdays
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

   write(errstring,'(''year,mon,day,hour,min,sec'',6(1x,i4),'' not a valid date.'')') &
              year,month,day,ohours,ominutes,oseconds
   call error_handler(E_ERR,'set_date_gregorian',errstring,source)
   call error_handler(E_ERR,'set_date_gregorian',errstring,source)
endif

if(month /= 2 .and. day > days_per_month(month)) then
   write(errstring,*)'month (',month,') does not have ',day,' days.'
   call error_handler(E_ERR,'set_date_gregorian',errstring,source)
endif

! Is this a leap year? Gregorian calendar assigns each year evenly
! divisible by 4 that is not a century year unevenly divisible by 400
! as a leap-year. (i.e. 1700,1800,1900 are not leap-years, 2000 is)
leap=greg_leap(year)

! Finish checking for day specication errors
if(month == 2 .and. (day > 29 .or. ((.not. leap) .and. day > 28))) then
   write(errstring,*)'month (',month,') does not have ',day,' days in a non-leap year.'
   call error_handler(E_ERR,'set_date_gregorian',errstring,source)
endif

! compute number of leap years fully past since base_year

nleapyr = (year - base_year) / 4 - (year - base_year) / 100 + (year - base_year) / 400

! Count up days in this year
ndays = 0
do m=1,month-1
 ndays = ndays + days_per_month(m)
 if(leap .and. m == 2) ndays = ndays + 1
enddo

totseconds = oseconds + 60*(ominutes + 60 * ohours)
totdays    = day - 1 +  ndays + 365*(year - base_year - nleapyr) + 366*(nleapyr)

set_date_gregorian = set_time(totseconds, totdays)

end function set_date_gregorian



function set_date_julian(year, month, day, hours, minutes, seconds)
!------------------------------------------------------------------------
!
! Returns time corresponding to date for julian calendar.

integer, intent(in)           :: year, month, day
integer, intent(in), optional :: hours, minutes, seconds
type(time_type)               :: set_date_julian

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
   write(errstring,'(''year,mon,day,hour,min,sec'',6(1x,i4),'' not a valid date.'')') &
              year,month,day,ohours,ominutes,oseconds
   call error_handler(E_ERR,'set_date_julian',errstring,source)
endif

if(month /= 2 .and. day > days_per_month(month)) then
   write(errstring,*)'month (',month,') does not have ',day,' days.'
   call error_handler(E_ERR,'set_date_julian',errstring,source)
endif

! Is this a leap year? 
leap = (modulo(year,4) == 0)

! compute number of complete leap years from year 1
nleapyr = (year - 1) / 4

! Finish checking for day specication errors

if(month == 2 .and. (day > 29 .or. ((.not. leap) .and. day > 28))) then
   write(errstring,*)'non-leapyear month (',month,') does not have ',day,' days.'
   call error_handler(E_ERR,'set_date_julian',errstring,source)
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

integer, intent(in)           :: year, month, day
integer, intent(in), optional :: hours, minutes, seconds
type(time_type)               :: set_date_thirty

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
   write(errstring,'(''year,mon,day,hour,min,sec'',6(1x,i4),'' not a valid date.'')') &
              year,month,day,ohours,ominutes,oseconds
   call error_handler(E_ERR,'set_date_thirty',errstring,source)
endif

set_date_thirty%days    = (day - 1) + 30 * ((month - 1) + 12 * (year - 1))
set_date_thirty%seconds = oseconds + 60 * (ominutes + 60 * ohours)

end function set_date_thirty



function set_date_no_leap(year, month, day, hours, minutes, seconds)
!------------------------------------------------------------------------
!
! Computes time corresponding to date for fixed 365 day year calendar.

integer, intent(in)           :: year, month, day
integer, intent(in), optional :: hours, minutes, seconds
type(time_type)               :: set_date_no_leap

integer :: totseconds, totdays
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

   write(errstring,'(''year,mon,day,hour,min,sec'',6(1x,i4),'' not a valid date.'')') &
              year,month,day,ohours,ominutes,oseconds
   call error_handler(E_ERR,'set_date_no_leap',errstring,source)
endif

ndays = 0
do m = 1, month - 1
   ndays = ndays + days_per_month(m)
enddo

totseconds = oseconds + 60 * (ominutes + 60 * ohours)
totdays    = day - 1 + ndays + 365 * (year - 1)

set_date_no_leap = set_time(totseconds, totdays)

end function set_date_no_leap


function set_date_gregorian_mars(year, month, day, hours, minutes, seconds)
!------------------------------------------------------------------------
!
! Computes time corresponding to date for gregorian MARS calendar.
! wgl - we need this to easily map from MarsWRF calendar to time_type
!         NOTE :: "month" has no meaning for the Mars calendar
!         ALSO :: the first day of the year is day 1, not 0
!         PlanetWRF wrfout files follow :: YYYY-DDDDD_HH:MM:SS
!
! To be clear, the first Mars year is:
!   According to MarsWRF         :: 0001-00001_00:00:00  thru  0001-00669_23:59:59
!   According to DART time_type  :: days=0, secs=0       thru  days=668, secs=86399

integer, intent(in)           :: year, month, day
integer, intent(in), optional :: hours, minutes, seconds
type(time_type)               :: set_date_gregorian_mars

integer :: oseconds, ominutes, ohours

! "base_year" for Mars will be defined as 1 (earliest wrfout file has year = 1)
integer :: base_year = 1

if ( .not. module_initialized ) call time_manager_init

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  ))   ohours = hours

! Need to check for bogus dates

! Do not bother checking month bounds because it has no meaning (though it should
!   be set to 1)
if(    oseconds > 59 .or. oseconds < 0 .or. &
       ominutes > 59 .or. ominutes < 0 .or. &
       ohours   > 23 .or. ohours   < 0 .or. &
                          day      < 1 .or. &
                          year     < base_year) then
   write(errstring,'(''year,mon,day,hour,min,sec'',6(1x,i4),'' not a valid date.'')') &
              year,month,day,ohours,ominutes,oseconds
   call error_handler(E_ERR,'set_date_gregorian_mars',errstring,source)
endif

! Month has no meaning for Mars calendar

! Currently there is no leap year defined for Mars

set_date_gregorian_mars%seconds = oseconds + 60*(ominutes + 60 * ohours)
set_date_gregorian_mars%days = day - 1 + 669*(year - base_year)

end function set_date_gregorian_mars

function set_date_solar_mars(year, month, day, hours, minutes, seconds)
!------------------------------------------------------------------------
!
! Computes time corresponding to date for gregorian MARS calendar.
! wgl - we need this to easily map from MarsWRF calendar to time_type
!         NOTE :: "month" has no meaning for the Mars calendar
!         ALSO :: the first day of the year is day 1, not 0
!         PlanetWRF wrfout files follow :: YYYY-DDDDD_HH:MM:SS
!
!cl number of days in year is now 100000, and there can only be 1 year.
!cl valid until 2050 Earth date, by which time WRF will not be used.

integer, intent(in)           :: year, month, day
integer, intent(in), optional :: hours, minutes, seconds
type(time_type)               :: set_date_solar_mars

integer :: oseconds, ominutes, ohours

! "base_year" for Mars will be defined as 1 (earliest wrfout file has year = 1)
integer :: base_year = 1

if ( .not. module_initialized ) call time_manager_init

! Missing optionals are set to 0

oseconds = 0; if(present(seconds)) oseconds = seconds
ominutes = 0; if(present(minutes)) ominutes = minutes
ohours   = 0; if(present(hours  ))   ohours = hours

! Need to check for bogus dates

! Do not bother checking month bounds because it has no meaning (though it should
!   be set to 1)
if(    oseconds > 59 .or. oseconds < 0 .or. &
       ominutes > 59 .or. ominutes < 0 .or. &
       ohours   > 23 .or. ohours   < 0 .or. &
                          day      < 1 .or. &
                          year     < base_year) then
   write(errstring,'(''year,mon,day,hour,min,sec'',6(1x,i4),'' not a valid date.'')') &
              year,month,day,ohours,ominutes,oseconds
   call error_handler(E_ERR,'set_date_solar_mars',errstring,source)
endif

! Month has no meaning for Mars calendar

! Currently there is no leap year defined for Mars


set_date_solar_mars%seconds = oseconds + 60*(ominutes + 60 * ohours)

!cl set_date_solar_mars%days = day - 1 + 669*(year - base_year)
!cl number of days in year changed to 100000 to allow MSD orbital calculation
!cl to allow simpler interface between data and model in the assimilation
set_date_solar_mars%days = day - 1 !this is always 0 -> + 100000*(year - base_year)

end function set_date_solar_mars

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

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: increment_date

integer :: oseconds, ominutes, ohours, odays, omonths, oyears

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
case(GREGORIAN_MARS)
   increment_date = increment_gregorian_mars(time, oyears, omonths, odays, &
                                     ohours, ominutes, oseconds)
case(SOLAR_MARS)
   increment_date = increment_solar_mars(time, oyears, omonths, odays, &
                                     ohours, ominutes, oseconds)
case default
   write(errstring,*)'calendar type (',calendar_type,') must be one of ', &
                      THIRTY_DAY_MONTHS ,GREGORIAN,JULIAN,NOLEAP,GREGORIAN_MARS
   call error_handler(E_ERR,'increment_date',errstring,source)
end select
end function increment_date



function increment_gregorian(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given time and some date increment, computes new time for gregorian calendar.

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: years, months, days, hours, minutes, seconds
type(time_type)                       :: increment_gregorian

integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cminute, chour, cday, cmonth, cyear

if ( .not. module_initialized ) call time_manager_init

call error_handler(E_ERR,'increment_gregorian','not implemented',source)

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
    call error_handler(E_ERR,'increment_gregorian',errstring,source)
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

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: increment_julian

integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cminute, chour, cday, cmonth, cyear, dyear
type(time_type) :: t

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
    call error_handler(E_ERR,'increment_julian',errstring,source)
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
      'month must not be incremented with other units',source)
endif

if( oyears /= 0 .and. &
    (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
     odays    /= 0 .or. omonths  /= 0) ) then
   call error_handler(E_ERR,'increment_julian', &
      'year must not be incremented with other units',source)
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

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: increment_thirty

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
    call error_handler(E_ERR,'increment_thirty',errstring,source)
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

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: increment_no_leap

integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cminute, chour, cday, cmonth, cyear, dyear
type(time_type) :: t

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
    call error_handler(E_ERR,'increment_no_leap',errstring,source)
endif

!  There are a number of other bad types of increments that should be
!  prohibited here; the addition is not associative
!  Easiest thing is to only let month and year be incremented by themselves
!  This is slight overkill since year is not really a problem.

if(omonths /= 0 .and. &
   (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
    odays    /= 0 .or. oyears   /= 0)) then
   call error_handler(E_ERR,'increment_no_leap', &
      'month must not be incremented with other units',source)
endif

if(oyears /= 0 .and. &
   (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
    odays    /= 0 .or. omonths  /= 0)) then
   call error_handler(E_ERR,'increment_no_leap', &
      'year must not be incremented with other units',source)
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


function increment_gregorian_mars(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given time and some date increment, computes new time for gregorian MARS calendar.

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: increment_gregorian_mars

!integer :: oseconds, ominutes, ohours, odays, omonths, oyears
!integer :: csecond, cminute, chour, cday, cmonth, cyear

if ( .not. module_initialized ) call time_manager_init

call error_handler(E_ERR,'increment_gregorian_mars','not implemented',source)

! FIXME: set a return value to avoid compiler warnings
! This code is never reached ... the error_handler() terminates before this. 
increment_gregorian_mars = time
if (present(years)) &
write(*,*)'stop complaining about',years,months,days,hours,minutes,seconds

end function increment_gregorian_mars

function increment_solar_mars(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given time and some date increment, computes new time for gregorian MARS calendar.

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: increment_solar_mars

!integer :: oseconds, ominutes, ohours, odays, omonths, oyears
!integer :: csecond, cminute, chour, cday, cmonth, cyear

if ( .not. module_initialized ) call time_manager_init

call error_handler(E_ERR,'increment_solar_mars','not implemented',source)

! FIXME: set a return value to avoid compiler warnings
! This code is never reached ... the error_handler() terminates before this. 
increment_solar_mars = time

end function increment_solar_mars


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

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: decrement_date

integer :: oseconds, ominutes, ohours, odays, omonths, oyears

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
case(GREGORIAN_MARS)
   decrement_date = decrement_gregorian_mars(time, oyears, omonths, odays, &
                       ohours, ominutes, oseconds)
case(SOLAR_MARS)
   decrement_date = decrement_solar_mars(time, oyears, omonths, odays, &
                       ohours, ominutes, oseconds)
case default
   write(errstring,*)'calendar type (',calendar_type,') not allowed.'
   call error_handler(E_ERR,'decrement_date',errstring,source)
end select

end function decrement_date



function decrement_gregorian(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given time and some date decrement, computes new time for gregorian calendar.

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: decrement_gregorian

integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cminute, chour, cday, cmonth, cyear

if ( .not. module_initialized ) call time_manager_init

call error_handler(E_ERR,'decrement_gregorian','not implemented',source)

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
    call error_handler(E_ERR,'decrement_gregorian',errstring,source)
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

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: decrement_julian

integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cminute, chour, cday, cmonth, cyear
type(time_type) :: t

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
    call error_handler(E_ERR,'decrement_julian',errstring,source)
endif

!  There are a number of other bad types of decrements that should be
!  prohibited here; the subtraction is not associative
!  Easiest thing is to only let month and year be decremented by themselves
!  This is slight overkill since year is not really a problem.

if(omonths /= 0 .and. (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
   odays  /= 0 .or. oyears /= 0)) then
   call error_handler(E_ERR,'decrement_julian', &
      'month must not be decremented with other units',source)
endif

if(oyears /= 0 .and. (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
   odays  /= 0 .or. omonths /= 0)) then
   call error_handler(E_ERR,'decrement_julian', &
      'year must not be decremented with other units',source)
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
   call error_handler(E_ERR,'decrement_julian',errstring,source)
endif

! Convert this back into a time
decrement_julian = set_date_julian(cyear, cmonth, cday, chour, cminute, csecond)

end function decrement_julian



function decrement_thirty(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given a time and some date decrement, computes new time for thirty day months.

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: decrement_thirty

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
    call error_handler(E_ERR,'decrement_thirty',errstring,source)
endif

csecond = oseconds + 60 * (ominutes + 60 * ohours)
cday    = odays + 30 * (omonths + 12 * oyears)
decrement_thirty = decrement_time(time, csecond, cday)

end function decrement_thirty



function decrement_no_leap(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given time and some date decrement, computes new time for julian calendar.

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: decrement_no_leap

integer :: oseconds, ominutes, ohours, odays, omonths, oyears
integer :: csecond, cminute, chour, cday, cmonth, cyear
type(time_type) :: t

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
    call error_handler(E_ERR,'decrement_no_leap',errstring,source)
endif

!  There are a number of other bad types of decrements that should be
!  prohibited here; the subtraction is not associative
!  Easiest thing is to only let month and year be decremented by themselves
!  This is slight overkill since year is not really a problem.

if(omonths /= 0 .and. (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
   odays   /= 0 .or. oyears /= 0)) then
   call error_handler(E_ERR,'decrement_no_leap', &
      'month must not be decremented with other units',source)
endif

if(oyears /= 0 .and. (oseconds /= 0 .or. ominutes /= 0 .or. ohours /= 0 .or. &
   odays  /= 0 .or. omonths /= 0)) then
   call error_handler(E_ERR,'decrement_no_leap', &
      'year must not be decremented with other units',source)
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
   call error_handler(E_ERR,'decrement_no_leap',errstring,source)
endif

! Convert this back into a time

decrement_no_leap = set_date_no_leap(cyear, cmonth, cday, chour, cminute, csecond)

end function decrement_no_leap


function decrement_gregorian_mars(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given time and some date decrement, computes new time for gregorian MARS calendar.

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: decrement_gregorian_mars

!integer :: oseconds, ominutes, ohours, odays, omonths, oyears
!integer :: csecond, cminute, chour, cday, cmonth, cyear

if ( .not. module_initialized ) call time_manager_init

call error_handler(E_ERR,'decrement_gregorian','not implemented',source)

! FIXME: set a return value to avoid compiler warnings
! This code is never reached ... the error_handler() terminates before this. 
decrement_gregorian_mars = time
if (present(years)) &
write(*,*)'stop complaining about',years,months,days,hours,minutes,seconds

end function decrement_gregorian_mars

function decrement_solar_mars(time, years, months, days, hours, minutes, seconds)
!-------------------------------------------------------------------------
!
! Given time and some date decrement, computes new time for gregorian MARS calendar.

type(time_type), intent(in)           :: time
integer,         intent(in), optional :: seconds, minutes, hours, days, months, years
type(time_type)                       :: decrement_solar_mars

!integer :: oseconds, ominutes, ohours, odays, omonths, oyears
!integer :: csecond, cminute, chour, cday, cmonth, cyear

if ( .not. module_initialized ) call time_manager_init

call error_handler(E_ERR,'decrement_gregorian','not implemented',source)

! FIXME: set a return value to avoid compiler warnings
! This code is never reached ... the error_handler() terminates before this. 
decrement_solar_mars = time

end function decrement_solar_mars


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

type(time_type), intent(in) :: time
integer                     :: days_in_month


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
case(GREGORIAN_MARS)
   days_in_month = days_in_month_gregorian_mars(time)
case(SOLAR_MARS)
   days_in_month = days_in_month_solar_mars(time)
case default
   write(errstring,*)'Invalid calendar type (',calendar_type,')'
   call error_handler(E_ERR,'days_in_month',errstring,source)
end select
end function days_in_month



function days_in_month_gregorian(time)
!--------------------------------------------------------------------------
!
! Returns the number of days in a gregorian month.

type(time_type), intent(in) :: time
integer                     :: days_in_month_gregorian

if ( .not. module_initialized ) call time_manager_init

call error_handler(E_ERR,'days_in_month_gregorian','not implemented',source)
days_in_month_gregorian = -1

end function days_in_month_gregorian



function days_in_month_julian(time)
!--------------------------------------------------------------------------
!
! Returns the number of days in a julian month.

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

type(time_type), intent(in) :: time
integer                     :: days_in_month_thirty

if ( .not. module_initialized ) call time_manager_init

days_in_month_thirty = 30

end function days_in_month_thirty



function days_in_month_no_leap(time)
!--------------------------------------------------------------------------
! Returns the number of days in a 365 day year month.

type(time_type), intent(in) :: time
integer                     :: days_in_month_no_leap

integer :: seconds, minutes, hours, day, month, year

if ( .not. module_initialized ) call time_manager_init

call get_date_no_leap(time, year, month, day, hours, minutes, seconds)
days_in_month_no_leap= days_per_month(month)

end function days_in_month_no_leap


function days_in_month_gregorian_mars(time)
!--------------------------------------------------------------------------
!
! Returns the number of days in a gregorian MARS month.

type(time_type), intent(in) :: time
integer                     :: days_in_month_gregorian_mars

if ( .not. module_initialized ) call time_manager_init

call error_handler(E_ERR,'days_in_month_gregorian_mars', &
                  'not implemented; oh, and, dude, Mars has no months',source)
days_in_month_gregorian_mars = -1

end function days_in_month_gregorian_mars

function days_in_month_solar_mars(time)
!--------------------------------------------------------------------------
!
! Returns the number of days in a gregorian MARS month.

type(time_type), intent(in) :: time
integer                     :: days_in_month_solar_mars

if ( .not. module_initialized ) call time_manager_init

call error_handler(E_ERR,'days_in_month_solar_mars', &
                  'not implemented; Mars has no months',source)
days_in_month_solar_mars = -1

end function days_in_month_solar_mars


!==========================================================================
! END OF days_in_month BLOCK
!==========================================================================
! START OF leap_year BLOCK
!==========================================================================


function leap_year(time)
!--------------------------------------------------------------------------
! Is this date in a leap year for default calendar?

type(time_type), intent(in) :: time
logical                     :: leap_year

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
case(GREGORIAN_MARS)
   leap_year = leap_year_gregorian_mars(time)
case(SOLAR_MARS)
   leap_year = leap_year_solar_mars(time)
case default
   write(errstring,*)'invalid calendar type (',calendar_type,')'
   call error_handler(E_ERR,'leap_year',errstring,source)
end select
end function leap_year



function leap_year_gregorian(time)
!--------------------------------------------------------------------------
!
! Is this a leap year for gregorian calendar?

type(time_type), intent(in) :: time
logical                     :: leap_year_gregorian

integer :: iyear, imonth, iday, ihour, imin, isec

if ( .not. module_initialized ) call time_manager_init

call get_date_gregorian(time, iyear, imonth, iday, ihour, imin, isec)

leap_year_gregorian=greg_leap(iyear)

end function leap_year_gregorian



function leap_year_julian(time)
!--------------------------------------------------------------------------
!
! Is this a leap year for julian calendar?

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

type(time_type), intent(in) :: time
logical                     :: leap_year_thirty

if ( .not. module_initialized ) call time_manager_init

leap_year_thirty = .FALSE.

end function leap_year_thirty



function leap_year_no_leap(time)
!--------------------------------------------------------------------------
!
! Another tough one; no leap year returns false for leap year inquiry.

type(time_type), intent(in) :: time
logical                     :: leap_year_no_leap

if ( .not. module_initialized ) call time_manager_init

leap_year_no_leap = .FALSE.

end function leap_year_no_leap


function leap_year_gregorian_mars(time)
!--------------------------------------------------------------------------
!
! Is this a leap year for gregorian calendar?
! trick question: answer is always no.

type(time_type), intent(in) :: time
logical                     :: leap_year_gregorian_mars

if ( .not. module_initialized ) call time_manager_init

call error_handler(E_MSG,'leap_year_gregorian_mars', &
                  'not implemented; oh, and, dude, Mars has no leap year',source)

leap_year_gregorian_mars = .FALSE.

end function leap_year_gregorian_mars

function leap_year_solar_mars(time)
!--------------------------------------------------------------------------
!
! Is this a leap year for gregorian calendar?
! trick question: answer is always no.

type(time_type), intent(in) :: time
logical                     :: leap_year_solar_mars

if ( .not. module_initialized ) call time_manager_init

call error_handler(E_MSG,'leap_year_solar_mars', &
                  'not implemented; Mars has no leap year',source)

leap_year_solar_mars = .FALSE.

end function leap_year_solar_mars


!==========================================================================
! END OF leap_year BLOCK
!==========================================================================
! START OF length_of_year BLOCK
!==========================================================================


function length_of_year()
!--------------------------------------------------------------------------
!
! What is the length of the year for the default calendar type

type(time_type) :: length_of_year

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
case(GREGORIAN_MARS)
   length_of_year = length_of_year_gregorian_mars()
case(SOLAR_MARS)
   length_of_year = length_of_year_solar_mars()
case default
   write(errstring,*)'invalid calendar type (',calendar_type,')'
   call error_handler(E_ERR,'length_of_year',errstring,source)
end select
end function length_of_year



function length_of_year_thirty()
!--------------------------------------------------------------------------
!

type(time_type) :: length_of_year_thirty

if ( .not. module_initialized ) call time_manager_init

length_of_year_thirty = set_time(0, 360)

end function length_of_year_thirty



function length_of_year_gregorian()
!---------------------------------------------------------------------------

type(time_type) :: length_of_year_gregorian

if ( .not. module_initialized ) call time_manager_init

length_of_year_gregorian = set_time(0, 0)

call error_handler(E_ERR,'length_of_year_gregorian', &
       'not implemented',source)

end function length_of_year_gregorian



function length_of_year_julian()
!--------------------------------------------------------------------------

type(time_type) :: length_of_year_julian

if ( .not. module_initialized ) call time_manager_init

length_of_year_julian = set_time((24 / 4) * 60 * 60, 365)

end function length_of_year_julian



function length_of_year_no_leap()
!--------------------------------------------------------------------------

type(time_type) :: length_of_year_no_leap

if ( .not. module_initialized ) call time_manager_init

length_of_year_no_leap = set_time(0, 365)

end function length_of_year_no_leap



function length_of_year_gregorian_mars()
!---------------------------------------------------------------------------

type(time_type) :: length_of_year_gregorian_mars

if ( .not. module_initialized ) call time_manager_init

length_of_year_gregorian_mars = set_time(0, 669)

end function length_of_year_gregorian_mars

function length_of_year_solar_mars()
!---------------------------------------------------------------------------

type(time_type) :: length_of_year_solar_mars

if ( .not. module_initialized ) call time_manager_init

length_of_year_solar_mars = set_time(0, 669)

end function length_of_year_solar_mars

!==========================================================================
! END OF length_of_year BLOCK
!==========================================================================
! START OF days_in_year BLOCK
!==========================================================================


function days_in_year(time)
!--------------------------------------------------------------------------
!
! What is the number of days in this year for the default calendar type

type(time_type), intent(in) :: time
integer                     :: days_in_year

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
case(GREGORIAN_MARS)
   days_in_year = days_in_year_gregorian_mars(time)
case(SOLAR_MARS)
   days_in_year = days_in_year_solar_mars(time)
case default
   write(errstring,*)'invalid calendar type (',calendar_type,')'
   call error_handler(E_ERR,'days_in_year',errstring,source)
end select
end function days_in_year



function days_in_year_thirty(time)
!--------------------------------------------------------------------------

type(time_type), intent(in) :: time
integer                     :: days_in_year_thirty

if ( .not. module_initialized ) call time_manager_init

days_in_year_thirty = 360

end function days_in_year_thirty



function days_in_year_gregorian(time)
!---------------------------------------------------------------------------

type(time_type), intent(in) :: time
integer                     :: days_in_year_gregorian

if ( .not. module_initialized ) call time_manager_init

if(leap_year_gregorian(time)) then
   days_in_year_gregorian = 366
else
   days_in_year_gregorian = 365
endif

end function days_in_year_gregorian



function days_in_year_julian(time)
!--------------------------------------------------------------------------

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

type(time_type), intent(in) :: time
integer                     :: days_in_year_no_leap

if ( .not. module_initialized ) call time_manager_init

days_in_year_no_leap = 365

end function days_in_year_no_leap



function days_in_year_gregorian_mars(time)
!---------------------------------------------------------------------------

type(time_type), intent(in) :: time
integer                     :: days_in_year_gregorian_mars

if ( .not. module_initialized ) call time_manager_init

days_in_year_gregorian_mars = 669

end function days_in_year_gregorian_mars

function days_in_year_solar_mars(time)
!---------------------------------------------------------------------------

type(time_type), intent(in) :: time
integer                     :: days_in_year_solar_mars

if ( .not. module_initialized ) call time_manager_init

!CL METHOD, with realistic rotation
days_in_year_solar_mars = 100000 
! 100000 days is the maximum number of days representable in the format used by planetWRF for days
! we ignore the years and assume that 'days' represents the Mars Solar Date (MSD)
! from MSD we can calculate all required astronomical data using Allison and McEwan(2000)

end function days_in_year_solar_mars

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

if ( .not. module_initialized ) call time_manager_init

if( n < 1 .or. n > 12 ) then
   write(errstring,*)'Illegal month (',n,')'
   call error_handler(E_ERR,'month_name',errstring,source)
endif

month_name = months(n)

end function month_name

!==========================================================================

function greg_leap(iyear)

! private function that takes an actual year number, not a time structure,
! and returns true for leap and false for not.  this should not be made
! public, unless we make wrappers for all the other calendar types.

integer, intent(in) :: iyear
logical             :: greg_leap

greg_leap=(modulo(iyear,4).eq.0)
if((modulo(iyear,100).eq.0).and.(modulo(iyear,400).ne.0))then
   greg_leap=.false.
endif

end function greg_leap

!==========================================================================


function julian_day(year, month, day)

! Given a date, computes the day from the beginning of the year.

integer, intent(in) :: day, month, year
integer             :: julian_day

integer             :: m
logical             :: leap

if ( .not. module_initialized ) call time_manager_init

if (calendar_type /= GREGORIAN) then
   call error_handler(E_ERR,'julian_day', &
        'only implemented for GREGORIAN calendar',source)
endif

julian_day = 0
leap = greg_leap(year)

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

! if there is no calendar return without error and without
! printing anything. 
if (calendar_type == NO_CALENDAR) return

! prints the time to standard output (or optional iunit) as a date
! NOTE: there is no check for PE number

  unit_in = 6

  if (present(iunit)) unit_in = iunit

  call get_date (time,y,mo,d,h,m,s)

  ! print_date assumes an Earth calendar -- so check for calendar_type  
  ! mars day numbers can be 3 digits long since there are no months.
  if ( calendar_type == GREGORIAN_MARS ) then
     mon = 'sol'
     if (present(str)) then
        write (unit_in,10) trim(str)//' ', y,mon(1:3),' ',d,' ',h,':',m,':',s
     else
        write (unit_in,10)       'DATE: ', y,mon(1:3),' ',d,' ',h,':',m,':',s
     endif
10   format (a,i4,1x,a3,a1,i3.3,3(a1,i2.2))
  else if (calendar_type == SOLAR_MARS ) then
     mon = 'sol'
     if (present(str)) then
        write (unit_in,12) trim(str)//' ', y,mon(1:3),' ',d,' ',h,':',m,':',s
     else
        write (unit_in,12)       'DATE: ', y,mon(1:3),' ',d,' ',h,':',m,':',s
     endif
12   format (a,i4,1x,a3,a1,i5.5,3(a1,i2.2))

  ! if not Mars, then use Earth calendar
  else
     mon = month_name(mo)
     if (present(str)) then
        write (unit_in,11) trim(str)//' ', y,mon(1:3),' ',d,' ',h,':',m,':',s
     else
        write (unit_in,11)       'DATE: ', y,mon(1:3),' ',d,' ',h,':',m,':',s
     endif
11   format (a,i4,1x,a3,4(a1,i2.2))
  end if

end subroutine print_date



function read_time(file_unit, form, ios_out)
!--------------------------------------------------------------------------------
!

integer,          intent(in)            :: file_unit
character(len=*), intent(in),  optional :: form
integer,          intent(out), optional :: ios_out

type(time_type)   :: read_time
integer           :: secs, days, ios

character(len=128) :: filename
logical :: is_named
integer :: rc

if ( .not. module_initialized ) call time_manager_init

if (ascii_file_format(form)) then
   read(file_unit, *, iostat=ios) secs, days
else
   read(file_unit, iostat=ios) secs, days
endif

if ( ios /= 0 ) then

   ! If ios_out argument is present, just return a non-zero ios
   if(present(ios_out)) then
      ios_out = ios
      return
   endif

   ! try to extract filename associated with unit to give context
   ! for which file we are trying to read from.
   inquire(file_unit, named=is_named, name=filename, iostat=rc)
   if ((rc /= 0) .or. (.not. is_named)) filename = 'unknown'

   ! Otherwise, read error is fatal, print message and stop
   call dump_unit_attributes(file_unit)   ! TJH DEBUG statement
   write(errstring,*)'read returned status ', ios, 'from input file ', trim(filename)
   call error_handler(E_ERR,'read_time',errstring,source)
else
   read_time = set_time(secs, days)
   if(present(ios_out)) ios_out = 0
endif

end function read_time



subroutine write_time(file_unit, time, form, ios_out)
!------------------------------------------------------------------------
! The time is expected to be written as either 
! an unformatted binary (form == "unformatted")
! or free-format ascii (form /= "unformatted")
! If ios_out is specified, do not call error handler here, 
! but return so a better context message can be generated.

integer,          intent(in)            :: file_unit
type(time_type),  intent(in)            :: time
character(len=*), intent(in),  optional :: form
integer,          intent(out), optional :: ios_out

integer            :: secs, days, io
character(len=129) :: filename
logical :: is_named
integer :: rc

if ( .not. module_initialized ) call time_manager_init

call get_time(time, secs, days)
io = 0

if (ascii_file_format(form)) then
   write(file_unit,'(i6,1x,i10)', iostat = io) secs, days
else
   write(file_unit, iostat = io) secs, days
endif

! return code for write.  0 = good, anything else bad.
if (present(ios_out)) then
   ios_out = io
endif

! if the caller is not asking for the return code, error out here.
! otherwise, return so caller can print out a better error message
! about which time it is trying to write.
if ((io /= 0) .and. (.not. present(ios_out))) then

   ! try to extract filename associated with unit to give context
   ! for which file we are trying to write to.
   inquire(file_unit, named=is_named, name=filename, iostat=rc)
   if ((rc /= 0) .or. (.not. is_named)) filename = 'unknown'

   write(errstring,*)'write returned status ', io, 'from file ', trim(filename)
   call error_handler(E_ERR,'write_time',errstring,source)
endif

end subroutine write_time


subroutine interactive_time(time)
!------------------------------------------------------------------------

type(time_type), intent(inout) :: time
integer                        :: second, minute, hour, day, month, year

if ( .not. module_initialized ) call time_manager_init

if (calendar_type == GREGORIAN) then
   write(*, *) 'input date (as integers): year month day hour minute second'
   read(*, *) year, month, day, hour, minute, second
   time = set_date(year, month, day, hour, minute, second)
elseif (calendar_type == GREGORIAN_MARS) then
   write(*, *) 'input date (as integers): year day hour minute second'
   read(*, *) year, day, hour, minute, second
   ! NOTE: "month" has no meaning for Mars calendar
   month = 1
   time = set_date(year, month, day, hour, minute, second)
elseif (calendar_type == SOLAR_MARS) then
   write(*, *) 'input date (as integers): year day hour minute second'
   read(*, *) year, day, hour, minute, second
   ! NOTE: "month" has no meaning for Mars calendar
   month = 1
   time = set_date(year, month, day, hour, minute, second)
else
   write(*, *) 'input time in days and seconds (as integers)'
   read(*, *) day, second
   time = set_time(second, day)
endif

end subroutine interactive_time

!-------------------------------------------------------------------------

function generate_seed(timestamp)

! given a time type, turn it into as unique an integer as possible.
! expected to be used to seed a random number generator in a way
! that you can reproduce the same sequence if seeded again from
! the same time value.
!
! the return value needs to be an i4 since seeds are only i4. 
! compute total number of seconds using a double integer (i8) and 
! return the least significant 32 bits.  a simple assignment could
! overflow an i4, and the seed needs to be as unique as possible
! so preserving the least significant digits is the better choice.

type(time_type), intent(in) :: timestamp
integer                     :: generate_seed

integer                :: days,seconds
integer(i8), parameter :: secs_day = 86400_i8

if ( .not. module_initialized ) call time_manager_init

call get_time(timestamp, seconds, days)

generate_seed = iand((secs_day * days) + seconds, int(z'00000000FFFFFFFF',i8))

end function generate_seed

!-------------------------------------------------------------------------

subroutine time_index_sort(t, index, num)

! Uses a heap sort alogrithm on t, returns array of sorted indices
! Sorts from earliest (smallest value) time to latest (largest value);
! uses time > and < operators to do the compare.  When index = 1
! that's the earliest time. 

integer,  intent(in)         :: num
type(time_type), intent(in)  :: t(num)
integer,  intent(out)        :: index(num)

integer :: ind, i, j, l_val_index, level 
type(time_type) :: l_val


if ( .not. module_initialized ) call time_manager_init

!  INITIALIZE THE INDEX ARRAY TO INPUT ORDER
do i = 1, num
   index(i) = i
end do

! Only one element, just send it back
if(num <= 1) return

level = num / 2 + 1
ind = num

! Keep looping until finished
do
   ! Keep going down levels until bottom
   if(level > 1) then
      level = level - 1
      l_val = t(index(level))
      l_val_index = index(level)
   else
      l_val = t(index(ind))
      l_val_index = index(ind)


      index(ind) = index(1)
      ind = ind - 1
      if(ind == 1) then
         index(1) = l_val_index
         return
      endif
   endif

   i = level
   j = 2 * level

   do while(j <= ind)
      if(j < ind) then
         if(t(index(j)) < t(index(j + 1))) j = j + 1
      endif
      if(l_val < t(index(j))) then
         index(i) = index(j)
         i = j
         j = 2 * j
      else
         j = ind + 1
      endif
      
   end do
   index(i) = l_val_index

end do

end subroutine time_index_sort


!-------------------------------------------------------------------------

end module time_manager_mod

