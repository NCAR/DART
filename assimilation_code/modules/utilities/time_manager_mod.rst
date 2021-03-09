MODULE time_manager_mod
=======================

Overview
--------

| Provides a set of routines to manipulate both time and calendars of various types.
| Time intervals are stored and defined in terms of integer number of days and integer seconds. The minimum time
  resolution is 1 second. Mathematical operations (e.g. addition, subtraction, multiplication) are defined on these
  intervals. Seconds which roll over 86400 (the number of seconds in a day) are converted into days.
| Calendars interpret time intervals in terms of years, months, days. Various calendars commonly in use in the
  scientific community are supported.

Other modules used
------------------

::

   types_mod
   utilities_mod

Public interfaces
-----------------

============================== ===================
*use time_manager_mod, only :* time_type
\                              operator(+)
\                              operator(-)
\                              operator(*)
\                              operator(/)
\                              operator(>)
\                              operator(>=)
\                              operator(==)
\                              operator(/=)
\                              operator(<)
\                              operator(<=)
\                              operator(//)
\                              set_time
\                              set_time_missing
\                              increment_time
\                              decrement_time
\                              get_time
\                              interval_alarm
\                              repeat_alarm
\                              THIRTY_DAY_MONTHS
\                              JULIAN
\                              GREGORIAN
\                              NOLEAP
\                              NO_CALENDAR
\                              GREGORIAN_MARS
\                              set_calendar_type
\                              get_calendar_type
\                              get_calendar_string
\                              set_date
\                              get_date
\                              increment_date
\                              decrement_date
\                              days_in_month
\                              leap_year
\                              length_of_year
\                              days_in_year
\                              month_name
\                              julian_day
\                              time_manager_init
\                              print_time
\                              print_date
\                              write_time
\                              read_time
\                              interactive_time
============================== ===================

| 

.. container:: routine

   *var = set_time(seconds [, days])*
   ::

      type(time_type)               :: set_time
      integer,           intent(in) :: seconds
      integer, optional, intent(in) :: days

.. container:: indent1

   Fills a time type. If seconds are > 86400, they are converted into the appropriate number of days. Note that seconds
   are specified first.

   =========== ================================================================================================
   ``seconds`` Number of seconds. If larger than 86400, they are converted into the appropriate number of days.
   *days*      Number of days. Default is 0.
   =========== ================================================================================================

| 

.. container:: routine

   *var = set_time_missing()*
   ::

      type(time_type)                        :: set_time_missing

.. container:: indent1

   Set a time type to a missing value. The resulting time value will cause an error if used for an arithmetic operation
   or if get_time() is called.

| 

.. container:: routine

   *var = increment_time(time, seconds [, days])*
   ::

      type(time_type)               :: increment_time
      type(time_type),   intent(in) :: time
      integer,           intent(in) :: seconds
      integer, optional, intent(in) :: days

.. container:: indent1

   Adds the specified number of seconds and optionally, days, to the given time and returns the new time. Increments
   cannot be negative (see decrement_time below).

   =========== =====================================================
   ``time``    time value to be incremented.
   ``seconds`` number of seconds to add to given time.
   *days*      optionally a number of days to add to the given time.
   =========== =====================================================

| 

.. container:: routine

   *var = decrement_time(time, seconds [, days])*
   ::

      type(time_type)                        :: decrement_time
      type(time_type), intent(in)            :: time
      integer,         intent(in)            :: seconds
      integer,         intent(in), optional  :: days

.. container:: indent1

   Subtract the specified number of seconds and optionally, days, to the given time and returns the new time. Decrements
   cannot be negative (see increment_time above).

   =========== ============================================================
   ``time``    time value to be decremented.
   ``seconds`` number of seconds to subtract from the given time.
   *days*      optionally a number of days to subtract from the given time.
   =========== ============================================================

| 

.. container:: routine

   *var = interval_alarm(time, time_interval, alarm, alarm_interval)*
   ::

      logical                        :: interval_alarm
      type(time_type), intent(in)    :: time
      type(time_type), intent(in)    :: time_interval
      type(time_type), intent(inout) :: alarm
      type(time_type), intent(in)    :: alarm_interval

.. container:: indent1

   Supports a commonly used type of test on times for models. Given the current time, and a time for an alarm,
   determines if this is the closest time to the alarm time given a time step of time_interval. If this is the closest
   time (alarm - time <= time_interval/2), the function returns true and the alarm is incremented by the alarm_interval.
   Watch for problems if the new alarm time is less than time + time_interval.

   ================== ==============================================================
   ``time``           Current time.
   ``time_interval``  Bin size for determining if alarm time is close enough to now.
   ``alarm``          When alarm next goes off next. Updated by this routine.
   ``alarm_interval`` How often alarm goes off.
   ================== ==============================================================

| 

.. container:: routine

   *var = repeat_alarm(time, alarm_frequency, alarm_length)*
   ::

      type(time_type)                :: repeat_alarm
      type(time_type), intent(in)    :: time
      type(time_type), intent(in)    :: alarm_frequency
      type(time_type), intent(in)    :: alarm_length

.. container:: indent1

   Repeat_alarm supports an alarm that goes off with alarm_frequency and lasts for alarm_length. If the nearest
   occurence of an alarm time is less than half an alarm_length from the input time, repeat_alarm is true. For instance,
   if the alarm_frequency is 1 day, and the alarm_length is 2 hours, then repeat_alarm is true from time 2300 on day n
   to time 0100 on day n + 1 for all n.

   =================== =============================
   ``time``            Current time.
   ``alarm_frequency`` How often the alarm goes off.
   ``alarm_length``    How long the alarm is true.
   =================== =============================

| 

.. container:: routine

   *var = get_calendar_type()*
   ::

      integer :: get_calendar_type

.. container:: indent1

   Returns default calendar type for mapping from time to date. Calendar types are public integer parameters that define
   various calendars. See elsewhere in this file for the list.

| 

.. container:: routine

   *var = set_date(year, month, day [, hours, minutes, seconds])*
   ::

      type(time_type)                :: set_date
      integer, intent(in)            :: year
      integer, intent(in)            :: month
      integer, intent(in)            :: day
      integer, intent(in), optional  :: hours
      integer, intent(in), optional  :: minutes
      integer, intent(in), optional  :: seconds

.. container:: indent1

   Given a date interpreted using the current calendar type, compute the corresponding time.

   ========= ==============================
   ``year``  Integer year.
   ``month`` Integer month number.
   ``day``   Integer day number.
   *hours*   Integer hour. Default is 0.
   *minutes* Integer minutes. Default is 0.
   *seconds* Integer seconds. Default is 0.
   ========= ==============================

| 

.. container:: routine

   *var = increment_date(time [, years, months, days, hours, minutes, seconds])*
   ::

      type(time_type)                :: increment_date
      type(time_type), intent(in)    :: time
      integer, intent(in), optional  :: years
      integer, intent(in), optional  :: months
      integer, intent(in), optional  :: days
      integer, intent(in), optional  :: hours
      integer, intent(in), optional  :: minutes
      integer, intent(in), optional  :: seconds

.. container:: indent1

   Given a time and some date increment, compute a new time. The interpretation of the date depends on the currently
   selected calendar type.

   ========= =====================================
   ``time``  Current time.
   *year*    Integer years to add. Default is 0.
   *month*   Integer months to add. Default is 0.
   *day*     Integer days to add. Default is 0.
   *hours*   Integer hours to add. Default is 0.
   *minutes* Integer minutes to add. Default is 0.
   *seconds* Integer seconds to add. Default is 0.
   ========= =====================================

| 

.. container:: routine

   *var = decrement_date(time [, years, months, days, hours, minutes, seconds])*
   ::

      type(time_type)                :: decrement_date
      type(time_type), intent(in)    :: time
      integer, intent(in), optional  :: years
      integer, intent(in), optional  :: months
      integer, intent(in), optional  :: days
      integer, intent(in), optional  :: hours
      integer, intent(in), optional  :: minutes
      integer, intent(in), optional  :: seconds

.. container:: indent1

   Given a time and some date decrement, compute a new time. The interpretation of the date depends on the currently
   selected calendar type.

   ========= ==========================================
   ``time``  Current time.
   *year*    Integer years to subtract. Default is 0.
   *month*   Integer months to subtract. Default is 0.
   *day*     Integer days to subtract. Default is 0.
   *hours*   Integer hours to subtract. Default is 0.
   *minutes* Integer minutes to subtract. Default is 0.
   *seconds* Integer seconds to subtract. Default is 0.
   ========= ==========================================

| 

.. container:: routine

   *var = days_in_month(time)*
   ::

      integer                        :: days_in_month
      type(time_type), intent(in)    :: time

.. container:: indent1

   Given a time, determine the month based on the currently selected calendar type and return the numbers of days in
   that month.

   ======== =============
   ``time`` Current time.
   ======== =============

| 

.. container:: routine

   *var = leap_year(time)*
   ::

      logical                        :: leap_year
      type(time_type),intent(in)     :: time

.. container:: indent1

   Given a time, determine if the current year is a leap year in the currently selected calendar type.

   ======== =============
   ``time`` Current time.
   ======== =============

| 

.. container:: routine

   *var = length_of_year()*
   ::

      integer                      :: length_of_year

.. container:: indent1

   For the currently selected calendar type, return the number of days in a year if that value is fixed (e.g. there are
   not leap years). For other calendar types, see days_in_year() which takes a time argument to determine the current
   year.

| 

.. container:: routine

   *var = days_in_year(time)*
   ::

      integer                        :: days_in_year
      type(time_type), intent(in)    :: time

.. container:: indent1

   Given a time, determine the year based on the currently selected calendar type and return the numbers of days in that
   year.

   ======== =============
   ``time`` Current time.
   ======== =============

| 

.. container:: routine

   *var = month_name(n)*
   ::

      character(len=9)               :: month_name
      integer,         intent(in)    :: n

.. container:: indent1

   Return a character string containing the month name corresponding to the given month number.

   ===== ==================================================
   ``n`` Month number. Must be between 1 and 12, inclusive.
   ===== ==================================================

| 

.. container:: routine

   *var = julian_day(year, month, day)*
   ::

      integer                        :: julian_day
      integer,        intent(in)     :: year
      integer,        intent(in)     :: month
      integer,        intent(in)     :: day

.. container:: indent1

   Given a date in year/month/day format, compute the day number from the beginning of the year. The currently selected
   calendar type must be GREGORIAN.

   ========= =======================================
   ``year``  Year number in the Gregorian calendar.
   ``month`` Month number in the Gregorian calendar.
   ``day``   Day of month in the Gregorian calendar.
   ========= =======================================

| 

.. container:: routine

   *var = read_time(file_unit [, form, ios_out])*
   ::

      type(time_type)                         :: read_time
      integer,          intent(in)            :: file_unit
      character(len=*), intent(in),  optional :: form
      integer,          intent(out), optional :: ios_out

.. container:: indent1

   Read a time from the given file unit number. The unit must already be open. The default format is ascii/formatted. If
   an error is encountered and ios_out is specified, the error status will be returned to the caller; otherwise the
   error is fatal.

   +---------------+-----------------------------------------------------------------------------------------------------+
   | ``file_unit`` | Integer file unit number of an already open file.                                                   |
   +---------------+-----------------------------------------------------------------------------------------------------+
   | ``form``      | Format to read the time. Options are 'formatted' or 'unformatted'. Default is 'formatted'.          |
   +---------------+-----------------------------------------------------------------------------------------------------+
   | ``ios_out``   | On error, if specified, the error status code is returned here. If not specified, an error calls    |
   |               | the standard error_handler and exits.                                                               |
   +---------------+-----------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call get_time(time, seconds [, days])*
   ::

      type(time_type), intent(in)             :: time
      integer,         intent(out)            :: seconds
      integer,         intent(out), optional  :: days

.. container:: indent1

   Returns days and seconds ( < 86400 ) corresponding to a time. If the optional 'days' argument is not given, the days
   are converted to seconds and the total time is returned as seconds. Note that seconds preceeds days in the argument
   list.

   =========== =======================================================================================================
   ``time``    Time to convert into seconds and days.
   ``seconds`` If days is specified, number of seconds in the current day. Otherwise, total number of seconds in time.
   ``days``    If specified, number of days in time.
   =========== =======================================================================================================

| 

.. container:: routine

   *call set_calendar_type(mytype)* or *call set_calendar_type(calstring)*
   ::

      integer, intent(in)               :: mytype
       or
      character(len=*), intent(in)      :: calstring

.. container:: indent1

   Selects the current calendar type, for converting between time and year/month/day. The argument can either be one of
   the predefined calendar integer parameter types (see elsewhere in this file for the list of types), or a string which
   matches the name of the integer parameters. The string interface is especially suitable for namelist use.

   ========== ==============================================
   ``mytype`` Integer parameter to select the calendar type.
   ========== ==============================================

   or

   ============= ======================================================================================================
   ``calstring`` Character string to select the calendar type. Valid strings match the names of the integer parameters.
   ============= ======================================================================================================

| 

.. container:: routine

   *call get_calendar_string(mystring)*
   ::

      character(len=*), intent(out)     :: mystring

.. container:: indent1

   Return the character string corresponding to the currently selected calendar type.

   ============ ============================================================
   ``mystring`` Character string corresponding to the current calendar type.
   ============ ============================================================

| 

.. container:: routine

   *call get_date(time, year, month, day, hour, minute, second)*
   ::

      type(time_type), intent(in)       :: time
      integer, intent(out)              :: year
      integer, intent(out)              :: month
      integer, intent(out)              :: day
      integer, intent(out)              :: hour
      integer, intent(out)              :: minute
      integer, intent(out)              :: second

.. container:: indent1

   Given a time, compute the corresponding date given the currently selected calendar type.

   ========== =============================
   ``time``   Input time.
   ``year``   Corresponding calendar year.
   ``month``  Corresponding calendar month.
   ``day``    Corresponding calendar day.
   ``hour``   Corresponding hour.
   ``minute`` Corresponding minute.
   ``second`` Corresponding second.
   ========== =============================

| 

.. container:: routine

   *call time_manager_init()*

.. container:: indent1

   Initializes any internal data needed by the time manager code. Does not need to be called before using any of the
   time manager routines; it will be called internally before executing any of the other routines.

| 

.. container:: routine

   *call print_time(time [, str, iunit])*
   ::

      type(time_type),  intent(in)           :: time
      character(len=*), intent(in), optional :: str
      integer,          intent(in), optional :: iunit
       

.. container:: indent1

   Print the time as days and seconds. If the optional str argument is specified, print that string as a label. If iunit
   is specified, write output to that unit; otherwise write to standard output/terminal.

   ======== =============================================================================
   ``time`` Time to be printed as days/seconds.
   *str*    String label to print before days/seconds. Default: 'TIME: '.
   *iunit*  Unit number to write output on. Default is standard output/terminal (unit 6).
   ======== =============================================================================

| 

.. container:: routine

   *call print_date(time [, str, iunit])*
   ::

      type(time_type),  intent(in)           :: time
      character(len=*), intent(in), optional :: str
      integer,          intent(in), optional :: iunit
       

.. container:: indent1

   Print the time as year/month/day/hour/minute/second, as computed from the currently selected calendar type. If the
   optional str argument is specified, print that string as a label. If iunit is specified, write output to that unit;
   otherwise write to standard output/terminal.

   ======== =============================================================================
   ``time`` Time to be printed as a calendar date/time.
   *str*    String label to print before date. Default: 'DATE: '.
   *iunit*  Unit number to write output on. Default is standard output/terminal (unit 6).
   ======== =============================================================================

| 

.. container:: routine

   *call write_time(file_unit, time [, form, ios_out])*
   ::

      integer,          intent(in)               :: file_unit
      type(time_type),  intent(in)               :: time
      character(len=*), intent(in),  optional    :: form
      integer,          intent(out), optional    :: ios_out

.. container:: indent1

   Write a time to an already open file unit. The optional 'form' argument controls whether it is formatted or
   unformatted. On error, the optional 'ios_out' argument returns the error code; otherwise a fatal error is triggered.

   +---------------+-----------------------------------------------------------------------------------------------------+
   | ``file_unit`` | Integer unit number for an already open file.                                                       |
   +---------------+-----------------------------------------------------------------------------------------------------+
   | ``time``      | Time to write to the file.                                                                          |
   +---------------+-----------------------------------------------------------------------------------------------------+
   | *form*        | String format specifier; either 'unformatted' or 'formatted'. Defaults to 'formatted'.              |
   +---------------+-----------------------------------------------------------------------------------------------------+
   | *ios_out*     | If specified, on error the i/o status error code is returned here. Otherwise, the standard error    |
   |               | handler is called and the program exits.                                                            |
   +---------------+-----------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call interactive_time(time)*
   ::

      type(time_type), intent(inout) :: time

.. container:: indent1

   Prompt the user for a time as a calendar date, based on the currently selected calendar type. Writes prompt to
   standard output and reads from standard input.

   ======== =========================
   ``time`` Time type to be returned.
   ======== =========================

| 

.. container:: type

   ::

      type time_type
         private
         integer :: seconds
         integer :: days
      end type time_type

.. container:: indent1

   This type is used to define a time interval.

| 

.. container:: type

   ::

       integer :: NO_CALENDAR
       integer :: GREGORIAN
       integer :: GREGORIAN_MARS
       integer :: JULIAN
       integer :: THIRTY_DAY_MONTHS
       integer :: NOLEAP

.. container:: indent1

   The public integer parameters which define different calendar types. The same names defined as strings can be used to
   set the calendar type.

| 

.. container:: type

   ::

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

.. container:: indent1

   Arithmetic operations are defined for time types, so expressions like

   ::

      t3 = t1 + t2

   | can be constructed. To use these operators, they must be listed on the module use statement in the form specified
     above.
   | Multiplication is one time and one scalar.
   | Division with a single slash is integer, and returns the largest integer for which time1 >= time2 \* n. Division
     with a double slash returns a double precision quotient of the two times.

| 

Namelist
--------

No namelist is currently defined for the time manager code.

Files
-----

-  none

References
----------

#. none

Private components
------------------

N/A
