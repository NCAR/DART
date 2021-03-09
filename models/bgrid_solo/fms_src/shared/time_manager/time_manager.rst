module time_manager_mod
=======================

Overview
--------

A software package that provides a set of simple interfaces for modelers to perform computations related to time and
dates.

.. container::

   The module defines a type that can be used to represent discrete times (accurate to one second) and to map these
   times into dates using a variety of calendars. A time is mapped to a date by representing the time with respect to an
   arbitrary base date (refer to <B>NOTES</B> section for the base date setting).
   The time_manager provides a single defined type, time_type, which is used to store time and date quantities. A
   time_type is a positive definite quantity that represents an interval of time. It can be most easily thought of as
   representing the number of seconds in some time interval. A time interval can be mapped to a date under a given
   calendar definition by using it to represent the time that has passed since some base date. A number of interfaces
   are provided to operate on time_type variables and their associated calendars. Time intervals can be as large as n
   days where n is the largest number represented by the default integer type on a compiler. This is typically
   considerably greater than 10 million years (assuming 32 bit integer representation) which is likely to be adequate
   for most applications. The description of the interfaces is separated into two sections. The first deals with
   operations on time intervals while the second deals with operations that convert time intervals to dates for a given
   calendar.

| 

Other modules used
------------------

.. container::

   ::

      fms_mod

Public interface
----------------

.. container::

   ::

      use time_manager_mod [, only:  set_time,
                                     get_time,
                                     increment_time,
                                     decrement_time,
                                     time_gt,
                                     time_ge,
                                     time_lt,
                                     time_le,
                                     time_eq,
                                     time_ne,
                                     time_plus,
                                     time_minus,
                                     time_scalar_mult,
                                     scalar_time_mult,
                                     time_divide,
                                     time_real_divide,
                                     time_scalar_divide,
                                     interval_alarm,
                                     repeat_alarm,
                                     set_calendar_type,
                                     get_calendar_type,
                                     get_date,
                                     set_date,
                                     increment_date,
                                     decrement_date,
                                     days_in_month,
                                     leap_year,
                                     length_of_year,
                                     days_in_year,
                                     month_name,
                                     time_manager_init,
                                     print_time,
                                     print_date ]

   set_time:
      Given some number of seconds and days, returns the corresponding time_type.
   get_time:
      Given a time interval, returns the corresponding seconds and days.
   increment_time:
      Given a time and an increment of days and seconds, returns a time that adds this increment to an input time.
   decrement_time:
      Given a time and a decrement of days and seconds, returns a time that subtracts this decrement from an input time.
   time_gt:
      Returns true if time1 > time2.
   time_ge:
      Returns true if time1 >= time2.
   time_lt:
      Returns true if time1 < time2.
   time_le:
      Returns true if time1 <= time2.
   time_eq:
      Returns true if time1 == time2.
   time_ne:
      Returns true if time1 /= time2.
   time_plus:
      Returns sum of two time_types.
   time_minus:
      Returns difference of two time_types.
   time_scalar_mult:
      Returns time multiplied by integer factor n.
   scalar_time_mult:
      Returns time multiplied by integer factor n.
   time_divide:
      Returns the largest integer, n, for which time1 >= time2 \* n.
   time_real_divide:
      Returns the double precision quotient of two times.
   time_scalar_divide:
      Returns the largest time, t, for which n \* t <= time.
   interval_alarm:
      Given a time, and a time interval, this function returns true if this is the closest time step to the alarm time.
   repeat_alarm:
      Repeat_alarm supports an alarm that goes off with alarm_frequency and lasts for alarm_length.
   set_calendar_type:
      Sets the default calendar type for mapping time intervals to dates.
   get_calendar_type:
      Returns the value of the default calendar type for mapping from time to date.
   get_date:
      Given a time_interval, returns the corresponding date under the selected calendar.
   set_date:
      Given an input date in year, month, days, etc., creates a time_type that represents this time interval from the
      internally defined base date.
   increment_date:
      Increments the date represented by a time interval and the default calendar type by a number of seconds, etc.
   decrement_date:
      Decrements the date represented by a time interval and the default calendar type by a number of seconds, etc.
   days_in_month:
      Given a time interval, gives the number of days in the month corresponding to the default calendar.
   leap_year:
      Returns true if the year corresponding to the date for the default calendar is a leap year. Returns false for
      THIRTY_DAY_MONTHS and NO_LEAP.
   length_of_year:
      Returns the mean length of the year in the default calendar setting.
   days_in_year:
      Returns the number of days in the calendar year corresponding to the date represented by time for the default
      calendar.
   month_name:
      Returns a character string containing the name of the month corresponding to month number n.
   time_manager_init:
      Write the version information to the log file.
   print_time:
      Prints the given time_type argument as a time (using days and seconds).
   print_date:
      prints the time to standard output (or optional unit) as a date.

| 

Public data
-----------

.. container::

   +-----------+--------------+-------+-------+-----------------------------------------------------------------------+
   | Name      | Type         | Value | Units | Description                                                           |
   +===========+==============+=======+=======+=======================================================================+
   | time_type | derived type | ---   | ---   | Derived-type data variable used to store time and date quantities. It |
   |           |              |       |       | contains two PRIVATE variables: seconds and days.                     |
   +-----------+--------------+-------+-------+-----------------------------------------------------------------------+

Public routines
---------------

a. .. rubric:: Set_time
      :name: set_time

   ::

      <B> set_time </B>(seconds, days)

   **DESCRIPTION**
      Given some number of seconds and days, returns the corresponding time_type.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``seconds``                                               | A number of seconds (can be greater than 86400), must be  |
      |                                                           | positive.                                                 |
      |                                                           | [integer, dimension(scalar)]                              |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``days``                                                  | A number of days, must be positive.                       |
      |                                                           | [integer, dimension(scalar)]                              |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | A time interval corresponding to this number of days and  |
      |                                                           | seconds.                                                  |
      |                                                           | [, dimension]                                             |
      +-----------------------------------------------------------+-----------------------------------------------------------+

b. .. rubric:: Get_time
      :name: get_time

   ::

      call get_time </B>(time, seconds, days)

   **DESCRIPTION**
      Given a time interval, returns the corresponding seconds and days.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | A time interval.                                          |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``seconds``                                               | A number of seconds (< 86400).                            |
      |                                                           | [integer, dimension(scalar)]                              |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``days``                                                  | A number of days, must be positive.                       |
      |                                                           | [integer, dimension(scalar)]                              |
      +-----------------------------------------------------------+-----------------------------------------------------------+

c. .. rubric:: Increment_time
      :name: increment_time

   ::

       
      increment_time (time, seconds, days)

   **DESCRIPTION**
      Given a time and an increment of days and seconds, returns a time that adds this increment to an input time.
      Increments a time by seconds and days; increments cannot be negative.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``seconds``                                               | Increment of seconds (can be greater than 86400); must be |
      |                                                           | positive.                                                 |
      |                                                           | [integer, dimension(scalar)]                              |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``days``                                                  | Increment of days; must be positive.                      |
      |                                                           | [integer, dimension(scalar)]                              |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | A time that adds this increment to the input time.        |
      |                                                           | [, dimension]                                             |
      +-----------------------------------------------------------+-----------------------------------------------------------+

d. .. rubric:: Decrement_time
      :name: decrement_time

   ::

       
      decrement_time (time, seconds, days)

   **DESCRIPTION**
      Decrements a time by seconds and days; decrements cannot be negative.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``seconds``                                               | Decrement of seconds (can be greater than 86400); must be |
      |                                                           | positive.                                                 |
      |                                                           | [integer, dimension(scalar)]                              |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``days``                                                  | Decrement of days; must be positive.                      |
      |                                                           | [integer, dimension(scalar)]                              |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | A time that subtracts this decrement from an input time.  |
      |                                                           | If the result is negative, it is considered a fatal       |
      |                                                           | error.                                                    |
      |                                                           | [, dimension]                                             |
      +-----------------------------------------------------------+-----------------------------------------------------------+

e. .. rubric:: Time_gt
      :name: time_gt

   ::

      <B> time_gt </B>(time1, time2)

   **DESCRIPTION**
      Returns true if time1 > time2.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time1``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time2``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | Returns true if time1 > time2                             |
      |                                                           | [logical, dimension]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

f. .. rubric:: Time_ge
      :name: time_ge

   ::

      <B> time_ge </B>(time1, time2)

   **DESCRIPTION**
      Returns true if time1 >= time2.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time1``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time2``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | Returns true if time1 >= time2                            |
      |                                                           | [logical, dimension]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

g. .. rubric:: Time_lt
      :name: time_lt

   ::

      <B> time_lt </B>(time1, time2)

   **DESCRIPTION**
      Returns true if time1 < time2.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time1``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time2``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | Returns true if time1 < time2                             |
      |                                                           | [logical, dimension]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

h. .. rubric:: Time_le
      :name: time_le

   ::

      <B> time_le </B>(time1, time2)

   **DESCRIPTION**
      Returns true if time1 <= time2.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time1``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time2``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | Returns true if time1 <= time2                            |
      |                                                           | [logical, dimension]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

i. .. rubric:: Time_eq
      :name: time_eq

   ::

      <B> time_eq </B>(time1, time2)

   **DESCRIPTION**
      Returns true if time1 == time2.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time1``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time2``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | Returns true if time1 == time2                            |
      |                                                           | [logical, dimension]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

j. .. rubric:: Time_ne
      :name: time_ne

   ::

      <B> time_ne </B>(time1, time2)

   **DESCRIPTION**
      Returns true if time1 /= time2.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time1``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time2``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | Returns true if time1 /= time2                            |
      |                                                           | [logical, dimension]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

k. .. rubric:: Time_plus
      :name: time_plus

   ::

      <B> time_plus </B>(time1, time2)

   **DESCRIPTION**
      Returns sum of two time_types.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time1``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time2``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | Returns sum of two time_types.                            |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

l. .. rubric:: Time_minus
      :name: time_minus

   ::

      <B> time_minus </B>(time1, time2)

   **DESCRIPTION**
      Returns difference of two time_types. WARNING: a time type is positive so by definition time1 - time2 is the same
      as time2 - time1.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time1``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time2``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | Returns difference of two time_types.                     |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

m. .. rubric:: Time_scalar_mult
      :name: time_scalar_mult

   ::

      <B> time_scalar_mult </B>(time, n)

   **DESCRIPTION**
      Returns time multiplied by integer factor n.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | A time interval.                                          |
      |                                                           | [integer, dimension]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | Returns time multiplied by integer factor n.              |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

n. .. rubric:: Scalar_time_mult
      :name: scalar_time_mult

   ::

      <B> scalar_time_mult </B>(n, time)

   **DESCRIPTION**
      Returns time multiplied by integer factor n.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | An integer.                                               |
      |                                                           | [integer, dimension]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | Returns time multiplied by integer factor n.              |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

o. .. rubric:: Time_divide
      :name: time_divide

   ::

      <B> time_divide </B>(time1, time2)

   **DESCRIPTION**
      Returns the largest integer, n, for which time1 >= time2 \* n.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time1``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time2``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | Returns the largest integer, n, for which time1 >= time2  |
      |                                                           | \* n.                                                     |
      |                                                           | [integer, dimension]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

p. .. rubric:: Time_real_divide
      :name: time_real_divide

   ::

      <B> time_real_divide </B>(time1, time2)

   **DESCRIPTION**
      Returns the double precision quotient of two times.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time1``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time2``                                                 | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | Returns the double precision quotient of two times        |
      |                                                           | [integer, dimensiondouble precision]                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

q. .. rubric:: Time_scalar_divide
      :name: time_scalar_divide

   ::

      <B> time_scalar_divide </B>(time, n)

   **DESCRIPTION**
      Returns the largest time, t, for which n \* t <= time.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | An integer factor.                                        |
      |                                                           | [integer, dimension]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | Returns the largest time, t, for which n \* t <= time.    |
      |                                                           | [integer, dimensiondouble precision]                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

r. .. rubric:: Interval_alarm
      :name: interval_alarm

   ::

       
      interval_alarm (time, time_interval, alarm, alarm_interval)

   **DESCRIPTION**
      This is a specialized operation that is frequently performed in models. Given a time, and a time interval, this
      function is true if this is the closest time step to the alarm time. The actual computation is:
      if((alarm_time - time) <= (time_interval / 2))
      If the function is true, the alarm time is incremented by the alarm_interval; WARNING, this is a featured side
      effect. Otherwise, the function is false and there are no other effects. CAUTION: if the alarm_interval is smaller
      than the time_interval, the alarm may fail to return true ever again. Watch for problems if the new alarm time is
      less than time + time_interval
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | Current time.                                             |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time_interval``                                         | A time interval.                                          |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``alarm_interval``                                        | A time interval.                                          |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **INPUT/OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``alarm``                                                 | An alarm time, which is incremented by the alarm_interval |
      |                                                           | if the function is true.                                  |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``interval_alarm``                                        | Returns either True or false.                             |
      |                                                           | [logical]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

s. .. rubric:: Repeat_alarm
      :name: repeat_alarm

   ::

       
      repeat_alarm 

   **DESCRIPTION**
      Repeat_alarm supports an alarm that goes off with alarm_frequency and lasts for alarm_length. If the nearest
      occurence of an alarm time is less than half an alarm_length from the input time, repeat_alarm is true. For
      instance, if the alarm_frequency is 1 day, and the alarm_length is 2 hours, then repeat_alarm is true from time
      2300 on day n to time 0100 on day n + 1 for all n.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | Current time.                                             |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``alarm_frequency``                                       | A time interval for alarm_frequency.                      |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``alarm_length``                                          | A time interval for alarm_length.                         |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``repeat_alarm``                                          | Returns either True or false.                             |
      |                                                           | [logical]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

t. .. rubric:: Set_calendar_type
      :name: set_calendar_type

   ::

      call set_calendar_type (type)

   **DESCRIPTION**
      A constant number for setting the calendar type.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``type``                                                  | A constant number for setting the calendar type.          |
      |                                                           | [integer, dimension]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``calendar_type``                                         | A constant number for default calendar type.              |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **NOTE**
      At present, four integer constants are defined for setting the calendar type: THIRTY_DAY_MONTHS, JULIAN, NO_LEAP,
      and GREGORIAN. However, GREGORIAN CALENDAR is not completely implemented. Selection of this type will result in
      illegal type error.

u. .. rubric:: Get_calendar_type
      :name: get_calendar_type

   ::

       
      get_calendar_type ()

   **DESCRIPTION**
      There are no arguments in this function. It returns the value of the default calendar type for mapping from time
      to date.

v. .. rubric:: Get_date
      :name: get_date

   ::

      call get_date (time, year, month, day, hour, minute, second)

   **DESCRIPTION**
      Given a time_interval, returns the corresponding date under the selected calendar.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | A time interval.                                          |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``day``                                                   | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``month``                                                 | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``year``                                                  | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``second``                                                | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``minute``                                                | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``hour``                                                  | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **NOTE**
      For all but the thirty_day_months calendar, increments to months and years must be made separately from other
      units because of the non-associative nature of the addition. All the input increments must be positive.

w. .. rubric:: Set_date
      :name: set_date

   ::

       
      set_date (year, month, day, hours, minutes, seconds)

   **DESCRIPTION**
      Given a date, computes the corresponding time given the selected date time mapping algorithm. Note that it is
      possible to specify any number of illegal dates; these should be checked for and generate errors as appropriate.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | A time interval.                                          |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``day``                                                   | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``month``                                                 | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``year``                                                  | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``second``                                                | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``minute``                                                | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``hour``                                                  | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``set_date``                                              | A time interval.                                          |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

x. .. rubric:: Increment_date
      :name: increment_date

   ::

       
      increment_date (time, years, months, days, hours, minutes, seconds)

   **DESCRIPTION**
      Given a time and some date increment, computes a new time. Depending on the mapping algorithm from date to time,
      it may be possible to specify undefined increments (i.e. if one increments by 68 days and 3 months in a Julian
      calendar, it matters which order these operations are done and we don't want to deal with stuff like that, make it
      an error).
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | A time interval.                                          |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``day``                                                   | An increment of days.                                     |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``month``                                                 | An increment of months.                                   |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``year``                                                  | An increment of years.                                    |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``second``                                                | An increment of seconds.                                  |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``minute``                                                | An increment of minutes.                                  |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``hour``                                                  | An increment of hours.                                    |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``increment_date``                                        | A new time based on the input time interval and the       |
      |                                                           | default calendar type.                                    |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

y. .. rubric:: Decrement_date
      :name: decrement_date

   ::

       
      decrement_date (time, years, months, days, hours, minutes, seconds)

   **DESCRIPTION**
      Given a time and some date decrement, computes a new time. Depending on the mapping algorithm from date to time,
      it may be possible to specify undefined decrements (i.e. if one decrements by 68 days and 3 months in a Julian
      calendar, it matters which order these operations are done and we don't want to deal with stuff like that, make it
      an error).
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | A time interval.                                          |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``day``                                                   | A decrement of days.                                      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``month``                                                 | A deincrement of months.                                  |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``year``                                                  | A deincrement of years.                                   |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``second``                                                | A deincrement of seconds.                                 |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``minute``                                                | A deincrement of minutes.                                 |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``hour``                                                  | A deincrement of hours.                                   |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``decrement_date``                                        | A new time based on the input time interval and the       |
      |                                                           | default calendar type.                                    |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **NOTE**
      For all but the thirty_day_months calendar, decrements to months and years must be made separately from other
      units because of the non-associative nature of addition. All the input decrements must be positive. If the result
      is a negative time (i.e. date before the base date) it is considered a fatal error.

z. .. rubric:: Days_in_month
      :name: days_in_month

   ::

      <B> days_in_month (time)

   **DESCRIPTION**
      Given a time, computes the corresponding date given the selected date time mapping algorithm.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``days_in_month``                                         | The number of days in the month given the selected time   |
      |                                                           | mapping algorithm.                                        |
      |                                                           | [integer, dimension]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

a. .. rubric:: Leap_year
      :name: leap_year

   ::

       
      leap_year (time)

   **DESCRIPTION**
      Is this date in a leap year for default calendar? Returns true if the year corresponding to the date for the
      default calendar is a leap year. Returns false for THIRTY_DAY_MONTHS and NO_LEAP.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | A time interval.                                          |
      |                                                           | [time_type, dimension]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``leap_year``                                             | True if the year corresponding to the date for the        |
      |                                                           | default calendar is a leap year. False for                |
      |                                                           | THIRTY_DAY_MONTHS and NO_LEAP and otherwise.              |
      |                                                           | [calendar_type, dimension]                                |
      +-----------------------------------------------------------+-----------------------------------------------------------+

b. .. rubric:: Length_of_year
      :name: length_of_year

   ::

       
      length_of_year ()

   **DESCRIPTION**
      There are no arguments in this function. It returns the mean length of the year in the default calendar setting.

c. .. rubric:: Days_in_year
      :name: days_in_year

   ::

       
      days_in_year ()

   **DESCRIPTION**
      Returns the number of days in the calendar year corresponding to the date represented by time for the default
      calendar.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | A time interval.                                          |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      == ==============================================================
      \  The number of days in this year for the default calendar type.
      == ==============================================================

d. .. rubric:: Month_name
      :name: month_name

   ::

       
      month_name (n)

   **DESCRIPTION**
      Returns a character string containing the name of the month corresponding to month number n. Definition is the
      same for all calendar types.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | Month number.                                             |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``month_name``                                            | The character string associated with a month. For now all |
      |                                                           | calendars have 12 months and will return standard names.  |
      |                                                           | [character]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

e. .. rubric:: Time_manager_init
      :name: time_manager_init

   ::

       
      time_manager_init ()

   **DESCRIPTION**
      Initialization routine. This routine does not have to be called, all it does is write the version information to
      the log file.

f. .. rubric:: Print_time
      :name: print_time

   ::

       
      print_time (time,str,unit)

   **DESCRIPTION**
      Prints the given time_type argument either as a time (using days and seconds). NOTE: there is no check for PE
      number.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | Time that will be printed.                                |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``str``                                                   | Character string that precedes the printed time or date.  |
      |                                                           | [character (len=*)]                                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``unit``                                                  | Unit number for printed output. The default unit is       |
      |                                                           | stdout.                                                   |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

g. .. rubric:: Print_date
      :name: print_date

   ::

       
      print_date (time,str,unit)

   **DESCRIPTION**
      Prints the given time_type argument as a date (using year,month,day, hour,minutes and seconds). NOTE: there is no
      check for PE number.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time``                                                  | Time that will be printed.                                |
      |                                                           | [time_type]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``str``                                                   | Character string that precedes the printed time or date.  |
      |                                                           | [character (len=*)]                                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``unit``                                                  | Unit number for printed output. The default unit is       |
      |                                                           | stdout.                                                   |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

Data sets
---------

.. container::

   None.

Error messages
--------------

.. container::

   None.

References
----------

.. container::

   None.

| 

Compiler specifics
------------------

.. container::

   None.

| 

Precompiler options
-------------------

.. container::

   None.

| 

Loader options
--------------

.. container::

   None.

Test PROGRAM
------------

.. container::

   time_main2
      ::

                 use time_manager_mod
                 implicit none
                 type(time_type) :: dt, init_date, astro_base_date, time, final_date
                 type(time_type) :: next_rad_time, mid_date
                 type(time_type) :: repeat_alarm_freq, repeat_alarm_length
                 integer :: num_steps, i, days, months, years, seconds, minutes, hours
                 integer :: months2, length
                 real :: astro_days
            
         Set calendar type
             call set_calendar_type(THIRTY_DAY_MONTHS)
                 call set_calendar_type(JULIAN)
             call set_calendar_type(NO_LEAP)
            
          Set timestep
                 dt = set_time(1100, 0)
            
          Set initial date
                 init_date = set_date(1992, 1, 1)
            
          Set date for astronomy delta calculation
                 astro_base_date = set_date(1970, 1, 1, 12, 0, 0)
            
          Copy initial time to model current time
                 time = init_date
            
          Determine how many steps to do to run one year
                 final_date = increment_date(init_date, years = 1)
                 num_steps = (final_date - init_date) / dt
                 write(*, *) 'Number of steps is' , num_steps
            
          Want to compute radiation at initial step, then every two hours
                 next_rad_time = time + set_time(7200, 0)
            
          Test repeat alarm
                 repeat_alarm_freq = set_time(0, 1)
                 repeat_alarm_length = set_time(7200, 0)
            
          Loop through a year
                 do i = 1, num_steps
            
          Increment time
                 time = time + dt
            
          Test repeat alarm
                 if(repeat_alarm(time, repeat_alarm_freq, repeat_alarm_length)) &
                 write(*, *) 'REPEAT ALARM IS TRUE'
            
          Should radiation be computed? Three possible tests.
          First test assumes exact interval; just ask if times are equal
              if(time == next_rad_time) then
          Second test computes rad on last time step that is <= radiation time
              if((next_rad_time - time) < dt .and. time < next_rad) then
          Third test computes rad on time step closest to radiation time
                  if(interval_alarm(time, dt, next_rad_time, set_time(7200, 0))) then
                    call get_date(time, years, months, days, hours, minutes, seconds)
                    write(*, *) days, month_name(months), years, hours, minutes, seconds
            
          Need to compute real number of days between current time and astro_base
                    call get_time(time - astro_base_date, seconds, days)
                    astro_days = days + seconds / 86400.
                write(*, *) 'astro offset ', astro_days
                 end if
            
          Can compute daily, monthly, yearly, hourly, etc. diagnostics as for rad
            
          Example: do diagnostics on last time step of this month
                 call get_date(time + dt, years, months2, days, hours, minutes, seconds)
                 call get_date(time, years, months, days, hours, minutes, seconds)
                 if(months /= months2) then
                    write(*, *) 'last timestep of month'
                    write(*, *) days, months, years, hours, minutes, seconds
                 endif
            
          Example: mid-month diagnostics; inefficient to make things clear
                 length = days_in_month(time)
                 call get_date(time, years, months, days, hours, minutes, seconds)
                 mid_date = set_date(years, months, 1) + set_time(0, length) / 2
            
                 if(time < mid_date .and. (mid_date - time) < dt) then
                    write(*, *) 'mid-month time'
                    write(*, *) days, months, years, hours, minutes, seconds
                 endif
            
                 end do

      end program time_main2

| 

Notes
-----

.. container::

   The Gregorian calendar type is not completely implemented, and currently no effort is put on it since it doesn't
   differ from Julian in use between 1901 and 2099.
   The <a name="base date">base date</a> is implicitly defined so users don't need to be concerned with it. For the
   curious, the base date is defined as 0 seconds, 0 minutes, 0 hours, day 1, month 1, year 1 for the Julian and
   thirty_day_months calendars, and 1 January, 1900, 0 seconds, 0 minutes, 0 hour for the Gregorian calendar.
   Please note that a time is a positive definite quantity.
   See the `Test Program <TEST%20PROGRAM>`__ for a simple program that shows some of the capabilities of the time
   manager.

| 
