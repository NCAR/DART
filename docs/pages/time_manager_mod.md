[]{#TOP}

MODULE time\_manager\_mod
=========================

  ------------------------------------------------------------------------ -------------------------------------------------------------------
  ![DART project logo](../../../docs/images/Dartboard7.png){height="70"}   Jump to [DART Documentation Main Index](../../../docs/index.html)
  ------------------------------------------------------------------------ -------------------------------------------------------------------

[INTERFACES](#Interface) / [NAMELIST](#Namelist) / [FILES](#FilesUsed) /
[REFERENCES](#References) / [ERRORS](#Errors) / [PLANS](#FuturePlans) /
[PRIVATE COMPONENTS](#PrivateComponents) / [TERMS OF USE](#Legalese)

Overview
--------

Provides a set of routines to manipulate both time and calendars of
various types.\
\
Time intervals are stored and defined in terms of integer number of days
and integer seconds. The minimum time resolution is 1 second.
Mathematical operations (e.g. addition, subtraction, multiplication) are
defined on these intervals. Seconds which roll over 86400 (the number of
seconds in a day) are converted into days.\
\
Calendars interpret time intervals in terms of years, months, days.
Various calendars commonly in use in the scientific community are
supported.

[]{#OtherModulesUsed}

<div class="top">

\[[](#top)\]

</div>

------------------------------------------------------------------------

OTHER MODULES USED
------------------

    types_mod
    utilities_mod

[](#Interface)

<div class="top">

\[[top](#)\]

</div>

------------------------------------------------------------------------

PUBLIC INTERFACES
-----------------

  | module name                   | interface                              |
  |------------------------------ | ----------------------------------------
  | use time_manager_mod, only :  |     [time\_type](#time_type)      
  | |                                   [operator(+)](#op_type)
  | |                                   [operator(-)](#op_type)
  | |                                   [operator(\*)](#op_type)
  | |                                   [operator(/)](#op_type)
  | |                                   [operator(&gt;)](#op_type)
  | |                                   [operator(&gt;=)](#op_type)
  | |                                   [operator(==)](#op_type)
  | |                                   [operator(/=)](#op_type)
  | |                                   [operator(&lt;)](#op_type)
  | |                                   [operator(&lt;=)](#op_type)
  | |                                   [operator(//)](#op_type)
  | |                                   [set\_time](#set_time)
  | |                                   [set\_time\_missing](#set_time_missing)
  | |                                   [increment\_time](#increment_time)
  | |                                   [decrement\_time](#decrement_time)
  | |                                   [get\_time](#get_time)
  | |                                   [interval\_alarm](#interval_alarm)
  | |                                   [repeat\_alarm](#repeat_alarm)
  | |                                   [THIRTY\_DAY\_MONTHS](#cal_type)
  | |                                   [JULIAN](#cal_type)
  | |                                   [GREGORIAN](#cal_type)
  | |                                   [NOLEAP](#cal_type)
  | |                                   [NO\_CALENDAR](#cal_type)
  | |                                   [GREGORIAN\_MARS](#cal_type)
  | |                                   [set\_calendar\_type](#set_calendar_type)
  | |                                   [get\_calendar\_type](#get_calendar_type)
  | |                                   [get\_calendar\_string](#get_calendar_string)
  | |                                   [set\_date](#set_date)
  | |                                   [get\_date](#get_date)
  | |                                   [increment\_date](#increment_date)
  | |                                   [decrement\_date](#decrement_date)
  | |                                   [days\_in\_month](#days_in_month)
  | |                                   [leap\_year](#leap_year)
  | |                                   [length\_of\_year](#length_of_year)
  | |                                   [days\_in\_year](#days_in_year)
  | |                                   [month\_name](#month_name)
  | |                                   [julian\_day](#julian_day)
  | |                                   [time\_manager\_init](#time_manager_init)
  | |                                   [print\_time](#print_time)
  | |                                   [print\_date](#print_date)
  | |                                   [write\_time](#write_time)
  | |                                   [read\_time](#read_time)
  | |                                   [interactive\_time](#interactive_time)
  | |                                   [generate\_seed](#generate_seed)



  ---------------------------------------------------------------------------------

set_time
--------

[](#set_time)

<div class="routine">

    var = set_time(seconds [, days])

    type(time_type)               :: set_time
    integer,           intent(in) :: seconds
    integer, optional, intent(in) :: days

</div>

<div class="indent1">

Fills a time type. If seconds are &gt; 86400, they are converted into
the appropriate number of days. Note that seconds are specified first.

  | argument  |  description  |
  |-----------|--------------------------------------------------------------------------------------------------
  |*seconds*  | Number of seconds. If larger than 86400, they are converted into the appropriate number of days.
  |*days*     | Number of days. Default is 0.

</div>


set_time_missing
----------------

[](#set_time_missing)

<div class="routine">

    var = set_time_missing()

    type(time_type)  :: set_time_missing

</div>

<div class="indent1">

Set a time type to a missing value. The resulting time value will cause
an error if used for an arithmetic operation or if get\_time() is
called.

</div>


\
[]{#increment_time}\

<div class="routine">

*var = increment\_time(time, seconds *\[, days\]*)*
    type(time_type)               :: increment_time
    type(time_type),   intent(in) :: time
    integer,           intent(in) :: seconds
    integer, optional, intent(in) :: days

</div>

<div class="indent1">

Adds the specified number of seconds and optionally, days, to the given
time and returns the new time. Increments cannot be negative (see
decrement\_time below).

  ----------- -------------------------------------------------------
  *time*      time value to be incremented.
  *seconds*   number of seconds to add to given time.
  *days*      optionally a number of days to add to the given time.
  ----------- -------------------------------------------------------

</div>

\
[]{#decrement_time}\

<div class="routine">

*var = decrement\_time(time, seconds *\[, days\]*)*
    type(time_type)                        :: decrement_time
    type(time_type), intent(in)            :: time
    integer,         intent(in)            :: seconds
    integer,         intent(in), optional  :: days

</div>

<div class="indent1">

Subtract the specified number of seconds and optionally, days, to the
given time and returns the new time. Decrements cannot be negative (see
increment\_time above).

  ----------- --------------------------------------------------------------
  *time*      time value to be decremented.
  *seconds*   number of seconds to subtract from the given time.
  *days*      optionally a number of days to subtract from the given time.
  ----------- --------------------------------------------------------------

</div>

\
[]{#interval_alarm}\

<div class="routine">

*var = interval\_alarm(time, time\_interval, alarm, alarm\_interval)*
    logical                        :: interval_alarm
    type(time_type), intent(in)    :: time
    type(time_type), intent(in)    :: time_interval
    type(time_type), intent(inout) :: alarm
    type(time_type), intent(in)    :: alarm_interval

</div>

<div class="indent1">

Supports a commonly used type of test on times for models. Given the
current time, and a time for an alarm, determines if this is the closest
time to the alarm time given a time step of time\_interval. If this is
the closest time (alarm - time &lt;= time\_interval/2), the function
returns true and the alarm is incremented by the alarm\_interval. Watch
for problems if the new alarm time is less than time + time\_interval.

  ------------------- ----------------------------------------------------------------
  *time*              Current time.
  *time\_interval*    Bin size for determining if alarm time is close enough to now.
  *alarm*             When alarm next goes off next. Updated by this routine.
  *alarm\_interval*   How often alarm goes off.
  ------------------- ----------------------------------------------------------------

</div>

\
[]{#repeat_alarm}\

<div class="routine">

*var = repeat\_alarm(time, alarm\_frequency, alarm\_length)*
    type(time_type)                :: repeat_alarm
    type(time_type), intent(in)    :: time
    type(time_type), intent(in)    :: alarm_frequency
    type(time_type), intent(in)    :: alarm_length

</div>

<div class="indent1">

Repeat\_alarm supports an alarm that goes off with alarm\_frequency and
lasts for alarm\_length. If the nearest occurence of an alarm time is
less than half an alarm\_length from the input time, repeat\_alarm is
true. For instance, if the alarm\_frequency is 1 day, and the
alarm\_length is 2 hours, then repeat\_alarm is true from time 2300 on
day n to time 0100 on day n + 1 for all n.

  -------------------- -------------------------------
  *time*               Current time.
  *alarm\_frequency*   How often the alarm goes off.
  *alarm\_length*      How long the alarm is true.
  -------------------- -------------------------------

</div>

\
[]{#get_calendar_type}\

<div class="routine">

*var = get\_calendar\_type()*
    integer :: get_calendar_type

</div>

<div class="indent1">

Returns default calendar type for mapping from time to date. Calendar
types are public integer parameters that define various calendars. See
elsewhere in this file for the list.

</div>

\
[]{#set_date}\

<div class="routine">

*var = set\_date(year, month, day *\[, hours, minutes, seconds\]*)*
    type(time_type)                :: set_date
    integer, intent(in)            :: year
    integer, intent(in)            :: month
    integer, intent(in)            :: day
    integer, intent(in), optional  :: hours
    integer, intent(in), optional  :: minutes
    integer, intent(in), optional  :: seconds

</div>

<div class="indent1">

Given a date interpreted using the current calendar type, compute the
corresponding time.

  ----------- --------------------------------
  *year*      Integer year.
  *month*     Integer month number.
  *day*       Integer day number.
  *hours*     Integer hour. Default is 0.
  *minutes*   Integer minutes. Default is 0.
  *seconds*   Integer seconds. Default is 0.
  ----------- --------------------------------

</div>

\
[]{#increment_date}\

<div class="routine">

*var = increment\_date(time *\[, years, months, days, hours, minutes,
seconds\]*)*
    type(time_type)                :: increment_date
    type(time_type), intent(in)    :: time
    integer, intent(in), optional  :: years
    integer, intent(in), optional  :: months
    integer, intent(in), optional  :: days
    integer, intent(in), optional  :: hours
    integer, intent(in), optional  :: minutes
    integer, intent(in), optional  :: seconds

</div>

<div class="indent1">

Given a time and some date increment, compute a new time. The
interpretation of the date depends on the currently selected calendar
type.

  ----------- ---------------------------------------
  *time*      Current time.
  *year*      Integer years to add. Default is 0.
  *month*     Integer months to add. Default is 0.
  *day*       Integer days to add. Default is 0.
  *hours*     Integer hours to add. Default is 0.
  *minutes*   Integer minutes to add. Default is 0.
  *seconds*   Integer seconds to add. Default is 0.
  ----------- ---------------------------------------

</div>

\
[]{#decrement_date}\

<div class="routine">

*var = decrement\_date(time *\[, years, months, days, hours, minutes,
seconds\]*)*
    type(time_type)                :: decrement_date
    type(time_type), intent(in)    :: time
    integer, intent(in), optional  :: years
    integer, intent(in), optional  :: months
    integer, intent(in), optional  :: days
    integer, intent(in), optional  :: hours
    integer, intent(in), optional  :: minutes
    integer, intent(in), optional  :: seconds

</div>

<div class="indent1">

Given a time and some date decrement, compute a new time. The
interpretation of the date depends on the currently selected calendar
type.

  ----------- --------------------------------------------
  *time*      Current time.
  *year*      Integer years to subtract. Default is 0.
  *month*     Integer months to subtract. Default is 0.
  *day*       Integer days to subtract. Default is 0.
  *hours*     Integer hours to subtract. Default is 0.
  *minutes*   Integer minutes to subtract. Default is 0.
  *seconds*   Integer seconds to subtract. Default is 0.
  ----------- --------------------------------------------

</div>

\
[]{#days_in_month}\

<div class="routine">

*var = days\_in\_month(time)*
    integer                        :: days_in_month
    type(time_type), intent(in)    :: time

</div>

<div class="indent1">

Given a time, determine the month based on the currently selected
calendar type and return the numbers of days in that month.

  -------- ---------------
  *time*   Current time.
  -------- ---------------

</div>

\
[]{#leap_year}\

<div class="routine">

*var = leap\_year(time)*
    logical                        :: leap_year
    type(time_type),intent(in)     :: time

</div>

<div class="indent1">

Given a time, determine if the current year is a leap year in the
currently selected calendar type.

  -------- ---------------
  *time*   Current time.
  -------- ---------------

</div>

\
[]{#length_of_year}\

<div class="routine">

*var = length\_of\_year()*
    integer                      :: length_of_year

</div>

<div class="indent1">

For the currently selected calendar type, return the number of days in a
year if that value is fixed (e.g. there are not leap years). For other
calendar types, see [days\_in\_year()](#days_in_year) which takes a time
argument to determine the current year.

</div>

\
[]{#days_in_year}\

<div class="routine">

*var = days\_in\_year(time)*
    integer                        :: days_in_year
    type(time_type), intent(in)    :: time

</div>

<div class="indent1">

Given a time, determine the year based on the currently selected
calendar type and return the numbers of days in that year.

  -------- ---------------
  *time*   Current time.
  -------- ---------------

</div>

\
[]{#month_name}\

<div class="routine">

*var = month\_name(n)*
    character(len=9)               :: month_name
    integer,         intent(in)    :: n

</div>

<div class="indent1">

Return a character string containing the month name corresponding to the
given month number.

  ----- ----------------------------------------------------
  *n*   Month number. Must be between 1 and 12, inclusive.
  ----- ----------------------------------------------------

</div>

\
[]{#julian_day}\

<div class="routine">

*var = julian\_day(year, month, day)*
    integer                        :: julian_day
    integer,        intent(in)     :: year
    integer,        intent(in)     :: month
    integer,        intent(in)     :: day

</div>

<div class="indent1">

Given a date in year/month/day format, compute the day number from the
beginning of the year. The currently selected calendar type must be
GREGORIAN.

  --------- -----------------------------------------
  *year*    Year number in the Gregorian calendar.
  *month*   Month number in the Gregorian calendar.
  *day*     Day of month in the Gregorian calendar.
  --------- -----------------------------------------

</div>

\
[]{#read_time}\

<div class="routine">

*var = read\_time(file\_unit *\[, form, ios\_out\]*)*
    type(time_type)                         :: read_time
    integer,          intent(in)            :: file_unit
    character(len=*), intent(in),  optional :: form
    integer,          intent(out), optional :: ios_out

</div>

<div class="indent1">

Read a time from the given file unit number. The unit must already be
open. The default format is ascii/formatted. If an error is encountered
and ios\_out is specified, the error status will be returned to the
caller; otherwise the error is fatal.

  -------------- -----------------------------------------------------------------------------------------------------------------------------------------
  *file\_unit*   Integer file unit number of an already open file.
  *form*         Format to read the time. Options are 'formatted' or 'unformatted'. Default is 'formatted'.
  *ios\_out*     On error, if specified, the error status code is returned here. If not specified, an error calls the standard error\_handler and exits.
  -------------- -----------------------------------------------------------------------------------------------------------------------------------------

</div>

\
[]{#get_time}\

<div class="routine">

*call get\_time(time, seconds *\[, days\]*)*
    type(time_type), intent(in)             :: time
    integer,         intent(out)            :: seconds
    integer,         intent(out), optional  :: days

</div>

<div class="indent1">

Returns days and seconds ( &lt; 86400 ) corresponding to a time. If the
optional 'days' argument is not given, the days are converted to seconds
and the total time is returned as seconds. Note that seconds preceeds
days in the argument list.

  ----------- ---------------------------------------------------------------------------------------------------------
  *time*      Time to convert into seconds and days.
  *seconds*   If days is specified, number of seconds in the current day. Otherwise, total number of seconds in time.
  *days*      If specified, number of days in time.
  ----------- ---------------------------------------------------------------------------------------------------------

</div>

\
[]{#set_calendar_type}\

<div class="routine">

*call set\_calendar\_type(mytype)* or
*call set\_calendar\_type(calstring)*
    integer, intent(in)               :: mytype
     or
    character(len=*), intent(in)      :: calstring

</div>

<div class="indent1">

Selects the current calendar type, for converting between time and
year/month/day. The argument can either be one of the predefined
calendar integer parameter types (see elsewhere in this file for the
list of types), or a string which matches the name of the integer
parameters. The string interface is especially suitable for namelist
use.

  ---------- ------------------------------------------------
  *mytype*   Integer parameter to select the calendar type.
  ---------- ------------------------------------------------

or
  ------------- --------------------------------------------------------------------------------------------------------
  *calstring*   Character string to select the calendar type. Valid strings match the names of the integer parameters.
  ------------- --------------------------------------------------------------------------------------------------------

</div>

\
[]{#get_calendar_string}\

<div class="routine">

*call get\_calendar\_string(mystring)*
    character(len=*), intent(out)     :: mystring

</div>

<div class="indent1">

Return the character string corresponding to the currently selected
calendar type.

  ------------ --------------------------------------------------------------
  *mystring*   Character string corresponding to the current calendar type.
  ------------ --------------------------------------------------------------

</div>

\
[]{#get_date}\

<div class="routine">

*call get\_date(time, year, month, day, hour, minute, second)*
    type(time_type), intent(in)       :: time
    integer, intent(out)              :: year
    integer, intent(out)              :: month
    integer, intent(out)              :: day
    integer, intent(out)              :: hour
    integer, intent(out)              :: minute
    integer, intent(out)              :: second

</div>

<div class="indent1">

Given a time, compute the corresponding date given the currently
selected calendar type.

  ---------- -------------------------------
  *time*     Input time.
  *year*     Corresponding calendar year.
  *month*    Corresponding calendar month.
  *day*      Corresponding calendar day.
  *hour*     Corresponding hour.
  *minute*   Corresponding minute.
  *second*   Corresponding second.
  ---------- -------------------------------

</div>

\
[]{#time_manager_init}\

<div class="routine">

*call time\_manager\_init()*

</div>

<div class="indent1">

Initializes any internal data needed by the time manager code. Does not
need to be called before using any of the time manager routines; it will
be called internally before executing any of the other routines.

</div>

\
[]{#print_time}\

<div class="routine">

*call print\_time(time *\[, str, iunit\]*)*
    type(time_type),  intent(in)           :: time
    character(len=*), intent(in), optional :: str
    integer,          intent(in), optional :: iunit


</div>

<div class="indent1">

Print the time as days and seconds. If the optional str argument is
specified, print that string as a label. If iunit is specified, write
output to that unit; otherwise write to standard output/terminal.

  --------- -------------------------------------------------------------------------------
  *time*    Time to be printed as days/seconds.
  *str*     String label to print before days/seconds. Default: 'TIME: '.
  *iunit*   Unit number to write output on. Default is standard output/terminal (unit 6).
  --------- -------------------------------------------------------------------------------

</div>

\
[]{#print_date}\

<div class="routine">

*call print\_date(time *\[, str, iunit\]*)*
    type(time_type),  intent(in)           :: time
    character(len=*), intent(in), optional :: str
    integer,          intent(in), optional :: iunit


</div>

<div class="indent1">

Print the time as year/month/day/hour/minute/second, as computed from
the currently selected calendar type. If the optional str argument is
specified, print that string as a label. If iunit is specified, write
output to that unit; otherwise write to standard output/terminal.

  --------- -------------------------------------------------------------------------------
  *time*    Time to be printed as a calendar date/time.
  *str*     String label to print before date. Default: 'DATE: '.
  *iunit*   Unit number to write output on. Default is standard output/terminal (unit 6).
  --------- -------------------------------------------------------------------------------

</div>

\
[]{#write_time}\

<div class="routine">

*call write\_time(file\_unit, time *\[, form, ios\_out\]*)*
    integer,          intent(in)               :: file_unit
    type(time_type),  intent(in)               :: time
    character(len=*), intent(in),  optional    :: form
    integer,          intent(out), optional    :: ios_out

</div>

<div class="indent1">

Write a time to an already open file unit. The optional 'form' argument
controls whether it is formatted or unformatted. On error, the optional
'ios\_out' argument returns the error code; otherwise a fatal error is
triggered.

  -------------- -------------------------------------------------------------------------------------------------------------------------------------------
  *file\_unit*   Integer unit number for an already open file.
  *time*         Time to write to the file.
  *form*         String format specifier; either 'unformatted' or 'formatted'. Defaults to 'formatted'.
  *ios\_out*     If specified, on error the i/o status error code is returned here. Otherwise, the standard error handler is called and the program exits.
  -------------- -------------------------------------------------------------------------------------------------------------------------------------------

</div>

\
[]{#interactive_time}\

<div class="routine">

*call interactive\_time(time)*
    type(time_type), intent(inout) :: time

</div>

<div class="indent1">

Prompt the user for a time as a calendar date, based on the currently
selected calendar type. Writes prompt to standard output and reads from
standard input.

  -------- ---------------------------
  *time*   Time type to be returned.
  -------- ---------------------------

</div>

[](#interactive_time)
<div class="routine">

    val = generate_seed(timestamp)

    type(time_type), intent(in) :: timestamp
    integer :: generate_seed

Given a time type, turn it into as unique an integer as possible.
Expected to be used to seed a random number generator in a way
that you can reproduce the same sequence if seeded again from
the same time value.

| argument | description |
-----------|--------------|
| timestamp | a dart time type |
| return | a unique integer seed |

</div>

\
[]{#time_type}\

<div class="type">

    type time_type
       private
       integer :: seconds
       integer :: days
    end type time_type

</div>

<div class="indent1">

This type is used to define a time interval.

</div>

\
[]{#cal_type}\

<div class="type">

     integer :: NO_CALENDAR
     integer :: GREGORIAN
     integer :: GREGORIAN_MARS
     integer :: JULIAN
     integer :: THIRTY_DAY_MONTHS
     integer :: NOLEAP

</div>

<div class="indent1">

The public integer parameters which define different calendar types. The
same names defined as strings can be used to set the calendar type.

</div>

\
[]{#op_type}\

<div class="type">

     operator(+)
     operator(-)
     operator(*)
     operator(/)
     operator(>)
     operator(>=)
     operator(==)
     operator(/=)
     operator(<)
     operator(<=)
     operator(//)

</div>

<div class="indent1">

Arithmetic operations are defined for time types, so expressions like

    t3 = t1 + t2

can be constructed. To use these operators, they must be listed on the
module use statement in the form specified above.\
\
Multiplication is one time and one scalar.\
\
Division with a single slash is integer, and returns the largest integer
for which time1 &gt;= time2 \* n. Division with a double slash returns a
double precision quotient of the two times.

</div>

\
[]{#Namelist}

<div class="top">

\[[top](#)\]

</div>

------------------------------------------------------------------------

NAMELIST
--------

No namelist is currently defined for the time manager code.

[]{#FilesUsed}

<div class="top">

\[[top](#)\]

</div>

------------------------------------------------------------------------

FILES
-----

-   none

[]{#References}

<div class="top">

\[[top](#)\]

</div>

------------------------------------------------------------------------

REFERENCES
----------

1.  none

[]{#Errors}

<div class="top">

\[[top](#)\]

</div>

------------------------------------------------------------------------

ERROR CODES and CONDITIONS
--------------------------

<div class="errors">

Routine
Message
Comment

</div>

KNOWN BUGS
----------

none at this time

[]{#FuturePlans}

<div class="top">

\[[top](#)\]

</div>

------------------------------------------------------------------------

FUTURE PLANS
------------

none at this time

[]{#PrivateComponents}

<div class="top">

\[[top](#)\]

</div>

------------------------------------------------------------------------

PRIVATE COMPONENTS
------------------

N/A

[]{#Legalese}

<div class="top">

\[[top](#)\]

</div>

------------------------------------------------------------------------

Terms of Use
------------

DART software - Copyright UCAR. This open source software is provided by
UCAR, "as is", without charge, subject to all terms of use at
<http://www.image.ucar.edu/DAReS/DART/DART_download>
