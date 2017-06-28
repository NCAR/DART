! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module dart_time_io_mod

!> \defgroup dart_time_io_mod dart_time_io_mod
!> Netcdf reading and writing dart model time.
!> Temporary module for dart time.
!>@todo should this go in state_vector_io_mod or io_filename_mod?
!> @{

use types_mod,        only : r8, digits12
use time_manager_mod, only : time_type, set_time, get_time, print_time, &
                             set_calendar_type, set_date, get_calendar_string, &
                             operator(+)

use utilities_mod,    only : nc_check, E_MSG, E_ERR, error_handler, to_upper

use typeSizes
use netcdf

implicit none
private

public :: read_model_time, write_model_time

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2, string3

contains


!--------------------------------------------------------------------
!> Make a stab at reading the time from the input file
!> and setting the calendar.  This routine assumes the file hasn't
!> been opened yet, and closes it when it's done.

function read_model_time(filename)

character(len=*), intent(in) :: filename
type(time_type)              :: read_model_time

integer :: ncid, ios, numdims, xtype, VarID
integer :: ntimes, seconds, days
integer :: year, month, day, hour, minute, second
type(time_type) :: base_time, delta_time

real(digits12) :: model_time, time_array(1)

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=256) :: file_calendar, dart_calendar
character(len=256) :: unitstring

! this is used in many error messages below.  set it here, and 
! don't reuse string3 here, please.
write(string3,*)'You may need to supply a model-specific "read_model_time()" to read the time.'


call nc_check( nf90_open(filename, NF90_NOWRITE, ncid), &
               'read_model_time',  'opening : "'//trim(filename)//'"')

ios = nf90_inq_varid(ncid, "time", VarID)
if (ios /= NF90_NOERR) then
   write(string1,*)'Expecting to be able to read the model time from a "time" variable'
   write(string2,*)'in file "'//trim(filename)//'"'
   call error_handler(E_ERR,'read_model_time', string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

ios = nf90_inquire_variable(ncid, VarID, xtype=xtype, dimids=dimIDs, ndims=numdims)
call nc_check(ios, 'read_model_time', 'inquire_variable "time" from "'//trim(filename)//'"')

if (numdims > 1) then
   write(string1,*)'Expecting the "time" variable to be a single dimension.'
   write(string2,*)'in file "'//trim(filename)//'"'
   call error_handler(E_ERR,'read_model_time', string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

! Since the time variable is known to have only 1 dimension, we know it is the first one.

ios = nf90_inquire_dimension(ncid, dimids(1), len=ntimes)
call nc_check(ios, 'read_model_time', 'inquire_dimension for time dimension from "'//trim(filename) )

! read the last one
ios = nf90_get_var(ncid, VarID, time_array, start=(/ntimes/), count=(/1/))
call nc_check(ios, 'read_model_time','get_var time' )
model_time = time_array(1)

! try to handle the calendar in a generic way

call get_calendar_string(dart_calendar)

ios = nf90_get_att(ncid, VarID, 'calendar', file_calendar)
if (ios /= NF90_NOERR) file_calendar = 'NO_CALENDAR'
call to_upper(file_calendar)

if (dart_calendar /= file_calendar ) then
   ! allow NO_CALENDAR, NO CALENDAR, and NONE to be synonyms.
   ! replace file_calendar with what dart uses to simplify the tests below.
   if (dart_calendar == 'NO_CALENDAR' .and. (file_calendar == 'NONE'        .or. &
                                             file_calendar == 'NO CALENDAR')) then
      file_calendar = 'NO_CALENDAR'
   else
      write(string1,*)'inconsistent calendar types between DART program and input file.'
      write(string2,*)'DART initialized with: ', trim(dart_calendar), ' File uses: ', trim(file_calendar)
      call error_handler(E_ERR, 'read_model_time:', string1, source,revision, revdate, &
                      text2=string2, text3=string3)
   endif
endif

if ( dart_calendar == 'NO_CALENDAR' ) then

   !> assumes time variable is real and fractional days.  if this isn't true,
   !> user has to supply their own read time routine.
   
   days    = floor(model_time)
   seconds = (model_time-real(days,digits12))*86400
   read_model_time = set_time(seconds, days)

else if ( dart_calendar == 'GREGORIAN' ) then

   ios = nf90_get_att(ncid, VarID, 'units', unitstring)
   call nc_check(ios, 'read_model_time', 'get_att time units "'//trim(filename)//'"')
   
   if (unitstring(1:10) == 'days since') then
   
      read(unitstring,'(11x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
      if (ios /= 0) then
         write(string1,*)'Unable to interpret time unit attribute. Error status was ',ios
         write(string2,*)'expected "days since YYYY-MM-DD HH:MM:SS", got "'//trim(unitstring)//'"'
         call error_handler(E_ERR, 'read_model_time:', string1, &
                source, revision, revdate, text2=string2, text3=string3)
      endif
   
      ! This is the start of their calendar
      base_time = set_date(year, month, day, hour, minute, second)
   
      days    = floor(model_time)
      seconds = (model_time-real(days,digits12))*86400

      delta_time = set_time(seconds, days)
      read_model_time = base_time + delta_time
   
   else if (unitstring(1:13) == 'seconds since') then
            
      read(unitstring,'(14x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
      if (ios /= 0) then
         write(string1,*)'Unable to interpret time unit attribute. Error status was ',ios
         write(string2,*)'expected "seconds since YYYY-MM-DD HH:MM:SS", got "'//trim(unitstring)//'"'
         call error_handler(E_ERR, 'read_model_time:', string1, &
                source, revision, revdate, text2=string2, text3=string3)
      endif
   
      ! This is the start of their calendar
      base_time  = set_date(year, month, day, hour, minute, second)
      delta_time = set_time(seconds)

      read_model_time = base_time + delta_time
   
   else
   
      write(string1, *) 'looking for "days since" or "seconds since" in the "units" attribute'
      call error_handler(E_ERR, 'read_model_time:', 'unable to set base time for gregorian calendar', &
                         source, revision, revdate, text2=string1, text3=string3)
   
   endif
else
   call error_handler(E_ERR, 'read_model_time:', &
    'calendar type "'//trim(dart_calendar)//' unsupported by default read_model_time() routine', &
                      source, revision, revdate, text2=string3)
endif

call print_time(read_model_time,'read_model_time')
!>@todo FIXME:
! make print_date() return without error if calendar is no_calendar,
! and then add a call to print_date() here.

call nc_check( nf90_close(ncid) , 'read_model_time closing : ', filename)

end function read_model_time


!--------------------------------------------------------------------
!> Write time to a netcdf file.
!> This routine assumes the file comes in already opened, and doesn't close it.

subroutine write_model_time(ncid, dart_time)

integer,             intent(in) :: ncid
type(time_type),     intent(in) :: dart_time

integer  :: ios
integer  :: xtype, numdims, ntimes
integer  :: VarID
integer  :: dart_days, dart_seconds

real(digits12) :: model_time

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character     (len=NF90_MAX_NAME)     :: dart_calendar, file_calendar

! this is used in many error messages below.  set it here, and
! don't reuse string3 here, please.
write(string3,*)'You may need to supply a model-specific "write_model_time()" to write the time.'

ios = nf90_inq_varid(ncid, "time", VarID)

! if the file doesn't already have a "time" variable, we make one
if (ios /= NF90_NOERR) then

   call error_handler(E_MSG,'write_model_time','no time variable found in file', &
              source, revision, revdate, text2='creating one')

   ! begin define mode
   ios = nf90_Redef(ncid)
   call nc_check(ios, "write_model_time", "redef")

   ! check to see if there is a time dimension
   ios = nf90_inq_dimid(ncid, "time", dimIds(1))

   ! if time dimension does not exist create it
   if (ios /= NF90_NOERR) then
      call nc_check(nf90_def_dim(ncid, "time", nf90_unlimited, dimIds(1)), &
        "write_model_time def_var dimension time")
   endif

   !>@todo NF90_UNLIMITED
   ios = nf90_def_var(ncid, name="time", xtype=nf90_int, varid=VarID)
   call nc_check(ios, "write_model_time", "time def_var")

   ! define time attributes consistent with CF convention
   ios = nf90_put_att(ncid, VarID, "long_name", "valid time of the model state")
   call nc_check(ios, "write_model_time", "time long_name")

   call get_calendar_string(dart_calendar)
   if (dart_calendar == 'NO_CALENDAR') then
      ios = nf90_put_att(ncid, VarID, "calendar", "none")
      call nc_check(ios, "write_model_time", "calendar long_name")

      ! ncview (actually, probably udunits2) crashes or errors out or 
      ! displays misleading plot axes if you use 'days since ...' as the units.
      ! if you simply use 'days' it works much better.

      ios = nf90_put_att(ncid, VarID, "units", "days")
      call nc_check(ios, "write_model_time", "units long_name")

   else if (dart_calendar == 'GREGORIAN') then
      ios = nf90_put_att(ncid, VarID, "calendar", "gregorian")
      call nc_check(ios, "write_model_time", "calendar long_name")

      ios = nf90_put_att(ncid, VarID, "units", "days since 1601-01-01 00:00:00")
      call nc_check(ios, "write_model_time", "units long_name")
   else
      call error_handler(E_ERR, 'write_model_time:', &
      'calendar type "'//trim(dart_calendar)//' unsupported by default write_model_time() routine', &
                      source, revision, revdate, text2=string3)
   endif

   ! end define mode
   call nc_check( nf90_Enddef(ncid),"write_model_time", "Enddef" )
endif

! See if the existing time dimension has a calendar and start date to consider

ios = nf90_get_att(ncid, VarID, 'calendar', file_calendar)
if (ios /= NF90_NOERR) file_calendar = 'NO_CALENDAR'
call to_upper(file_calendar)

if (dart_calendar /= file_calendar ) then
   ! allow NO_CALENDAR, NO CALENDAR, and NONE to be synonyms.
   ! replace file_calendar with what dart uses to simplify the tests below.
   if (dart_calendar == 'NO_CALENDAR' .and. (file_calendar == 'NONE'        .or. &
                                             file_calendar == 'NO CALENDAR')) then
      file_calendar = 'NO_CALENDAR'
   else
      write(string1,*)'inconsistent calendar types between DART program and input file.'
      write(string2,*)'DART initialized with: ', trim(dart_calendar), ' File uses: ', trim(file_calendar)
      call error_handler(E_ERR, 'write_model_time:', string1, source,revision, revdate, &
                      text2=string2, text3=string3)
   endif
endif

! need to know how long the time variable is and hammer the last time
ios = nf90_inquire_variable(ncid, VarID, xtype=xtype, dimids=dimIDs, ndims=numdims)
call nc_check(ios, 'write_model_time', 'inquire_variable "time"')

if (numdims > 1) then
   write(string1,*)'Expecting the "time" variable to be a single dimension.'
   call error_handler(E_ERR,'write_model_time', string1, &
              source, revision, revdate, text2=string3)
endif

! Since the time variable is known to have only 1 dimension, we know it is the first one.

ios = nf90_inquire_dimension(ncid, dimIds(1), len=ntimes)
call nc_check(ios, 'write_model_time', 'inquire_dimension for time dimension')

! convert time to something that netCDF can store, fractional days
call get_time(dart_time, dart_seconds, dart_days)
model_time = real(dart_days,digits12) + real(dart_seconds,digits12)/86400.0_digits12

! write dart days and seconds files to netcdf file
ios = nf90_put_var(ncid, VarID, model_time, start=(/ ntimes /))
call nc_check( ios, "write_model_time", "put_var model_time")

end subroutine write_model_time

!--------------------------------------------------------------------
!> @}
end module dart_time_io_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
