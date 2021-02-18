! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module schedule_mod

use        types_mod, only : missing_i, digits12

use    utilities_mod, only : error_handler, E_ERR, nmlfileunit, do_output, &
                             check_namelist_read, find_namelist_in_file, &
                             do_nml_file, do_nml_term

use time_manager_mod, only : time_type, set_calendar_type, get_calendar_type, &
                             set_time, set_date, get_time, print_time, print_date, &
                             operator(*), operator(+), operator(-),           &
                             operator(>), operator(<), operator(/),           &
                             operator(/=), operator(<=)

implicit none
private

!=======================================================================
!
!
!
!
!=======================================================================

character(len=*), parameter :: source = 'schedule_mod.f90'

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

!-----------------------------------------------------------------------
! Namelist with default values (put all possible obs into one file)
!-----------------------------------------------------------------------

integer, dimension(6) :: first_bin_start = (/ 1601, 1, 1, 0, 0, 0 /)
integer, dimension(6) :: first_bin_end   = (/ 2999, 1, 1, 0, 0, 0 /)
integer, dimension(6) :: last_bin_end    = (/ 2999, 1, 1, 0, 0, 0 /)
integer               :: bin_interval_days    = 1000000
integer               :: bin_interval_seconds = 0
integer               :: max_num_bins         = 1000
character(len=32)     :: calendar             = 'Gregorian'
logical               :: print_table          = .false.

namelist /schedule_nml/ first_bin_start, first_bin_end, last_bin_end, &
                        bin_interval_days, bin_interval_seconds, &
                        max_num_bins, calendar, print_table

!-----------------------------------------------------------------------

character(len = 129) :: msgstring

logical, save :: module_initialized = .false.


public :: schedule_type ! Module defines a single type

! Subroutines and functions operating on time_type
public :: set_regular_schedule, get_time_from_schedule, &
          get_schedule_length

interface get_time_from_schedule
   module procedure get_timetype_from_schedule
   module procedure get_realtime_from_schedule
end interface

!=======================================================================
contains
!=======================================================================


subroutine schedule_init()

integer :: iunit, io


! Read the namelist entry
call find_namelist_in_file("input.nml", "schedule_nml", iunit)
read(iunit, nml = schedule_nml, iostat = io)
call check_namelist_read(iunit, io, "schedule_nml")

! Write the namelist values to the log file
if (do_nml_file()) write(nmlfileunit, nml=schedule_nml)
if (do_nml_term()) write(     *     , nml=schedule_nml)

call set_calendar_type(calendar)

module_initialized = .true.

end subroutine schedule_init




subroutine set_regular_schedule(schedule)
! all of the necessary information to set a schedule is gotten from the namelist
type(schedule_type), intent(out) :: schedule


type(time_type) :: beg_time, end_time
type(time_type) :: bininterval, binwidth
type(time_type) :: TimeMin, TimeMax

integer :: iepoch, Nepochs, seconds, days

character(len=32) :: str1, str2, str3

if ( .not. module_initialized ) call schedule_init   ! reads the namelist

beg_time    = set_date(first_bin_start(1), first_bin_start(2), &
                       first_bin_start(3), first_bin_start(4), &
                       first_bin_start(5), first_bin_start(6) )
end_time    = set_date(first_bin_end(1), first_bin_end(2), &
                       first_bin_end(3), first_bin_end(4), &
                       first_bin_end(5), first_bin_end(6) )
TimeMax     = set_date(last_bin_end(1), last_bin_end(2), &
                       last_bin_end(3), last_bin_end(4), &
                       last_bin_end(5), last_bin_end(6) )

bininterval = set_time(bin_interval_seconds, bin_interval_days)

! do some error-checking first

if (end_time < beg_time) then
   write(msgstring,*)'schedule_nml:first_bin_end must be at or after first_bin_start'
   call error_handler(E_ERR,'set_regular_schedule',msgstring,source)
endif

if (TimeMax < end_time) then
   write(msgstring,*)'schedule_nml:last_bin_end must be at or after first_bin_end'
   call error_handler(E_ERR,'set_regular_schedule',msgstring,source)
endif

binwidth = beg_time - end_time
call get_time(binwidth, seconds, days)
binwidth = set_time(seconds, days)

if (bininterval < binwidth) then
   write(msgstring,*)'schedule_nml:bin interval must be >= bin width'
   call error_handler(E_ERR,'set_regular_schedule',msgstring,source)
endif

! FIXME
! Need to check that bin_interval is a multiple of the model advance step
! for real assimilation experiments. This is not needed for observation-space
! diagnostics or ...

! Determine temporal bin characteristics.
! The user input is not guaranteed to align on bin centers. 
! So -- we will assume the start time is correct and take strides till we
! get past the last time of interest. 
! Nepochs will be the total number of time intervals of the period requested.

TimeMin  = end_time
Nepochs  = 0
NepochLoop : do iepoch = 1,max_num_bins
   if ( TimeMin > TimeMax ) exit NepochLoop
   Nepochs = iepoch
   TimeMin = TimeMin + bininterval
enddo NepochLoop

if (do_output()) write(*,*)'Requesting ',Nepochs,' assimilation periods.'

if (Nepochs < 1) then
   write(msgstring,*)'schedule_nml:Requesting ZERO assimilation periods.'
   call error_handler(E_ERR,'set_regular_schedule',msgstring,source)
endif

! Now that we know the number of assimilation epochs, allocate and fill.
! Our assimilation bins start 1 second AFTER the stated time, and end ON
! the stated time. If you specify a bin from  00Z to 03Z ... and again
! from 03Z to 06Z ... the observation at 03Z is considered to be part
! of the bin from 00Z to 03Z ONLY. Mathematically, ( 00Z, 03Z ]

allocate(schedule%binstart(   Nepochs), schedule%binend(   Nepochs), &
         schedule%epoch_start(Nepochs), schedule%epoch_end(Nepochs))

schedule%binstart(1)    = beg_time + set_time(1,0)
schedule%binend(  1)    = end_time
schedule%binwidth       = binwidth
schedule%bininterval  = bininterval
schedule%current_bin    = 0
schedule%last_bin       = .false.
schedule%calendar       = get_calendar_type()
schedule%calendarstring = calendar
schedule%num_bins       = Nepochs

call get_time(schedule%binstart(1), seconds, days)
schedule%epoch_start(1) = days + seconds/86400.0_digits12

call get_time(schedule%binend(  1), seconds, days)
schedule%epoch_end(  1) = days + seconds/86400.0_digits12

BinLoop : do iepoch = 2,Nepochs

   schedule%binstart(iepoch) = schedule%binstart(iepoch-1) + bininterval
   schedule%binend(  iepoch) = schedule%binend(  iepoch-1) + bininterval

   call get_time(schedule%binstart(iepoch), seconds, days)
   schedule%epoch_start(iepoch) = days + seconds/86400.0_digits12

   call get_time(schedule%binend(  iepoch), seconds, days)
   schedule%epoch_end(  iepoch) = days + seconds/86400.0_digits12

enddo BinLoop

if ( print_table ) then
do iepoch = 1,Nepochs
   write(     *     ,*)
   write(str1,'(''epoch '',i6,''  start'')')iepoch
   write(str2,'(''epoch '',i6,'' center'')')iepoch
   write(str3,'(''epoch '',i6,''    end'')')iepoch

   call print_time( schedule%binstart(iepoch), str1)
   call print_time( schedule%binend(  iepoch), str3)

   call print_date( schedule%binstart(iepoch), str1)
   call print_date( schedule%binend(  iepoch), str3)
enddo
write(     *     ,*)
endif

end subroutine set_regular_schedule



subroutine get_timetype_from_schedule(mytime, schedule, iepoch, edge)

type(time_type),     intent(OUT) :: mytime
type(schedule_type), intent(IN)  :: schedule
integer,             intent(IN)  :: iepoch
integer, optional,   intent(IN)  :: edge

if (iepoch > schedule%num_bins) then
   write(msgstring,*)'schedule has ',schedule%num_bins,' bins; you wanted bin',iepoch
   call error_handler(E_ERR,'get_timetype_from_schedule',msgstring, source)
endif

if (present(edge)) then
   if (edge > 1) then
      mytime = schedule%binend(iepoch)
   else
      mytime = schedule%binstart(iepoch)
   endif
else
   mytime = schedule%binstart(iepoch)
endif

end subroutine



subroutine get_realtime_from_schedule(mytime, schedule, iepoch, edge)

real(digits12),      intent(OUT) :: mytime
type(schedule_type), intent(IN)  :: schedule
integer,             intent(IN)  :: iepoch
integer, optional,   intent(IN)  :: edge

type(time_type) :: yourtime
integer :: seconds, days

if (present(edge)) then
   call get_timetype_from_schedule(yourtime, schedule, iepoch, edge)
else
   call get_timetype_from_schedule(yourtime, schedule, iepoch)
endif

call get_time(yourtime, seconds, days)

mytime = days + seconds/86400.0_digits12

end subroutine


function get_schedule_length(schedule)
   type(schedule_type), intent(IN) :: schedule
   integer :: get_schedule_length

   get_schedule_length = schedule%num_bins
end function

end module schedule_mod

