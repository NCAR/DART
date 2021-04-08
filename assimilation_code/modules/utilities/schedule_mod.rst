MODULE schedule_mod
===================

Overview
--------

Provides a set of routines to generate a regular pattern of time windows. This module is only used for converting
observation sequences files to netCDF format. If it stands the test of time, it will likely be used to create an
assimilation schedule independent of the observation sequence file. 

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &schedule_nml
      first_bin_start      =  1601,  1,  1,  0,  0,  0
      first_bin_end        =  2999,  1,  1,  0,  0,  0
      last_bin_end         =  2999,  1,  1,  0,  0,  0
      bin_interval_days    = 1000000
      bin_interval_seconds = 0
      max_num_bins         = 1000
      calendar             = 'Gregorian'
      print_table          = .true.
     /

| 

The default values will cause (pretty much) all possible observations to be put into one output file.

.. container::

   +----------------------+-----------------------+---------------------------------------------------------------------+
   | Item                 | Type                  | Description                                                         |
   +======================+=======================+=====================================================================+
   | first_bin_start      | integer, dimension(6) | Date/time specification for starting time of first bin.             |
   +----------------------+-----------------------+---------------------------------------------------------------------+
   | first_bin_end        | integer, dimension(6) | Date/time specification for ending time of first bin. Sets the bin  |
   |                      |                       | width.                                                              |
   +----------------------+-----------------------+---------------------------------------------------------------------+
   | last_bin_end         | integer, dimension(6) | Date/time specification for ending time of last bin. Sets the       |
   |                      |                       | length of the overall time of the schedule.                         |
   +----------------------+-----------------------+---------------------------------------------------------------------+
   | bin_interval_days    | integer               | Sets the time between bins. Must be larger or equal to the bin      |
   |                      |                       | width.                                                              |
   +----------------------+-----------------------+---------------------------------------------------------------------+
   | bin_interval_seconds | integer               | Sets the time between bins. Must be larger or equal to the bin      |
   |                      |                       | width.                                                              |
   +----------------------+-----------------------+---------------------------------------------------------------------+
   | max_num_bins         | integer               | Upper limit on the number of bins.                                  |
   +----------------------+-----------------------+---------------------------------------------------------------------+
   | calendar             | character(len=32)     | String calendar type. Valid types are listed in the                 |
   |                      |                       | `time_manager_mod <time_manager_mod.html#cal_type>`__ file.         |
   +----------------------+-----------------------+---------------------------------------------------------------------+
   | print_table          | logical               | If .TRUE., print out information about the schedule each time       |
   |                      |                       | set_regular_schedule() is called.                                   |
   +----------------------+-----------------------+---------------------------------------------------------------------+

| 

Other modules used
------------------

::

   types_mod
   utilities_mod
   time_manager_mod

Public interfaces
-----------------

========================== ======================
*use schedule_mod, only :* schedule_type
\                          set_regular_schedule
\                          get_time_from_schedule
\                          get_schedule_length
========================== ======================

Namelist ``&schedule_mod_nml`` may be read from file ``input.nml``.

| 

.. container:: routine

   *call set_regular_schedule(schedule)*
   ::

      type(schedule_type), intent(out) :: schedule

.. container:: indent1

   Uses the namelist information to compute and fill a schedule_type variable.

   ============ ==========================================================================================================
   ``schedule`` Fills this derived type with the information needed to generate a series of regularly spaced time windows.
   ============ ==========================================================================================================

| 

.. container:: routine

   *call get_time_from_schedule(mytime, schedule, iepoch [, edge])*
   ::

      type(time_type),     intent(out) :: mytime
       or
      real(digits12),      intent(out) :: mytime
      type(schedule_type), intent(in)  :: schedule
      integer,             intent(in)  :: iepoch
      integer, optional,   intent(in)  :: edge

.. container:: indent1

   Returns either the leading or trailing time for the specified bin/epoch number for the given schedule. The time can
   be returned in one of two formats, depending on the variable type specified for the first argument: either a DART
   derived time_type, or a real of kind digits12 (defined in the types_mod).

   +--------------+------------------------------------------------------------------------------------------------------+
   | ``mytime``   | Return value with the leading or trailing edge time for the requested bin. There are two supported   |
   |              | return formats, either as a standard DART time_type, or as a real value which will contain the       |
   |              | number of days plus any fraction.                                                                    |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``schedule`` | Schedule type to extract information from.                                                           |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``iepoch``   | The bin number, or epoch number, to return a time for. Unless edge is specified and requests the     |
   |              | ending time, the time returned is the starting time for this bin.                                    |
   +--------------+------------------------------------------------------------------------------------------------------+
   | *edge*       | If specified, and if edge is larger than 1, the trailing edge time of the bin is returned. Any other |
   |              | value, or if this argument is not specified, returns the leading edge time of the bin.               |
   +--------------+------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *var = get_schedule_length()*
   ::

      integer                             :: get_schedule_length
      type(schedule_type), intent(in)     :: schedule

.. container:: indent1

   Return the total number of intervals/bins/epochs defined by this schedule.

   ============ =================================================
   ``schedule`` Return number of time intervals in this schedule.
   ============ =================================================

| 

.. container:: type

   ::

      type schedule_type
         private
         integer :: num_bins
         integer :: current_bin
         logical :: last_bin
         integer :: calendar
         character(len=32) :: calendarstring
         type(time_type)          :: binwidth
         type(time_type)          :: bininterval
         type(time_type), pointer :: binstart(   :) => NULL()
         type(time_type), pointer :: binend(     :) => NULL()
         real(digits12),  pointer :: epoch_start(:) => NULL()
         real(digits12),  pointer :: epoch_end(  :) => NULL()
      end type schedule_type

.. container:: indent1

   This type is used to define a schedule.

| 

Files
-----

========= =================================
filename  purpose
========= =================================
input.nml to read the schedule_mod namelist
========= =================================

References
----------

-  none

Private components
------------------

N/A
