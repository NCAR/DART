! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module dart_time_io_mod

!> \defgroup dart_time_io_mod dart_time_io_mod
!> Default routines for netCDF reading and writing dart model time.
!> If your model uses a different name for the time dimension
!> or has a different way of handing/storing time, it must provide 
!> a custom read_model_time() and write_model_time() routine.
!> @{

use types_mod,            only : r8, digits12
use time_manager_mod,     only : time_type, set_time, get_time, print_time, &
                                 set_calendar_type, set_date, get_calendar_string, &
                                 operator(+)

use utilities_mod,        only : E_MSG, E_ERR, error_handler, to_upper
use netcdf_utilities_mod, only : nc_check, nc_open_file_readonly, nc_close_file
use typeSizes
use netcdf

implicit none
private

public :: read_model_time, write_model_time

character(len=*), parameter :: source = 'dart_time_io_mod.f90'

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

character(len=*), parameter :: routine = 'read_model_time'
real(digits12) :: model_time, time_array(1)

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=256) :: file_calendar, dart_calendar
character(len=256) :: unitstring

! this is used in many error messages below.  set it here, and 
! don't reuse string3 here, please.
write(string3,*)'You may need to supply a model-specific "read_model_time()" to read the time.'


ncid = nc_open_file_readonly(filename, routine)

ios = nf90_inq_varid(ncid, "time", VarID)
if (ios /= NF90_NOERR) then
   write(string1,*)'Expecting to be able to read the model time from a "time" variable'
   write(string2,*)'in file "'//trim(filename)//'"'
   call error_handler(E_ERR,'read_model_time',string1,source,text2=string2,text3=string3)
endif

ios = nf90_inquire_variable(ncid, VarID, xtype=xtype, dimids=dimIDs, ndims=numdims)
call nc_check(ios, 'read_model_time', 'inquire_variable "time" from "'//trim(filename)//'"')

if (numdims > 1) then
   write(string1,*)'Expecting the "time" variable to be a single dimension.'
   write(string2,*)'in file "'//trim(filename)//'"'
   call error_handler(E_ERR,'read_model_time', string1, &
              source, text2=string2, text3=string3)
endif

if (numdims == 0) then
   ios = nf90_get_var(ncid, VarID, model_time)
   call nc_check(ios, 'read_model_time','get_var scalar time' )
else

   ! Since the time variable is known to have only 1 dimension, we know it is the first one.
   ios = nf90_inquire_dimension(ncid, dimIDs(1), len=ntimes)
   call nc_check(ios, 'read_model_time', 'inquire_dimension for time dimension from "'//trim(filename) )

   ! read the last one
   ios = nf90_get_var(ncid, VarID, time_array, start=(/ntimes/), count=(/1/))
   call nc_check(ios, 'read_model_time','get_var time' )
   model_time = time_array(1)
endif

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
      call error_handler(E_ERR, 'read_model_time:', string1, source, &
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
                source, text2=string2, text3=string3)
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
                source, text2=string2, text3=string3)
      endif
   
      ! This is the start of their calendar
      base_time  = set_date(year, month, day, hour, minute, second)
      delta_time = set_time(seconds)

      read_model_time = base_time + delta_time
   
   else
   
      write(string1, *) 'looking for "days since" or "seconds since" in the "units" attribute'
      call error_handler(E_ERR, 'read_model_time:', &
              'unable to set base time for gregorian calendar', &
              source, text2=string1, text3=string3)
   
   endif
else
   call error_handler(E_ERR, 'read_model_time:', &
    'calendar type "'//trim(dart_calendar)//'" unsupported by default read_model_time() routine', &
                      source, text2=string3)
endif

!>@todo FIXME: do we really want this to print from any
!> task without being asked?  i vote no.
call print_time(read_model_time,'read_model_time')

!>@todo FIXME:
! make print_date() return without error if calendar is no_calendar,
! and then add a call to print_date() here.  (also vote no.)

call nc_close_file(ncid, routine)

end function read_model_time


!--------------------------------------------------------------------
!> Write time to a netcdf file.
!> This routine assumes the file comes in already opened, and doesn't close it.

subroutine write_model_time(ncid, dart_time)

integer,             intent(in) :: ncid
type(time_type),     intent(in) :: dart_time

character(len=*), parameter :: routine = 'write_model_time'
integer  :: ios
integer  :: numdims, ntimes
integer  :: VarID
integer  :: dart_days, dart_seconds
integer  :: unlimitedDimId

logical :: has_unlimited, time_is_unlimited

real(digits12) :: model_time

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character     (len=NF90_MAX_NAME)     :: dart_calendar, var_calendar

! If there is no unlimited dimension, unlimitedDimID = -1
ios = nf90_inquire(ncid, unlimitedDimId=unlimitedDimId )
call nc_check(ios,routine,'checking unlimited dimension')

has_unlimited = (unlimitedDimID /= -1) 

! this is used in many error messages below.  set it here, and
! don't reuse string3 here, please.
write(string3,*)'You may need to supply a model-specific "write_model_time()" to write the time.'

! see what kind of calendar dart is currently running with.
call get_calendar_string(dart_calendar)

ios = nf90_inq_varid(ncid, "time", VarID)

! if the file doesn't already have a "time" variable, we make one
if (ios /= NF90_NOERR) then

   call error_handler(E_MSG, routine, 'no variable "time" found in file', &
              source, text2='creating one')

   ! begin define mode
   ios = nf90_Redef(ncid)
   call nc_check(ios, routine, "redef")

   ! check to see if there is a time dimension
   ! if it does not exist create it
   ios = nf90_inq_dimid(ncid, "time", dimIDs(1))
   if (ios /= NF90_NOERR) then

      ! If there is already an unlimited dimension, just make a
      ! time dimension of 'normal' size. If there is no unlimited dim already
      ! make the time variable 'unlimited'.

      if (has_unlimited) then
         ios = nf90_def_dim(ncid, "time", 1, dimIDs(1))
         call nc_check(ios, routine, 'def_dim singleton dimension time')
      else
         ios = nf90_def_dim(ncid, "time", nf90_unlimited, dimIDs(1))
         call nc_check(ios, routine, 'def_dim unlimited dimension time')
         has_unlimited = .true.
         unlimitedDimID = dimIDs(1)
      endif

   endif

   ! make the time variable be dimensioned time(time) which is the
   ! netCDF convention for coordinate variables (variables with the
   ! same name as a dimension).
   ios = nf90_def_var(ncid, name="time", xtype=nf90_double, dimids=dimIDs(1), varid=VarID)
   call nc_check(ios, routine, "time def_var")

   ! define time attributes consistent with CF convention
   ios = nf90_put_att(ncid, VarID, "long_name", "valid time of the model state")
   call nc_check(ios, routine, "time long_name")

   if (dart_calendar == 'NO_CALENDAR') then
      ios = nf90_put_att(ncid, VarID, "calendar", "none")
      call nc_check(ios, routine, "calendar long_name")

      ! ncview (actually, probably udunits2) crashes or errors out or 
      ! displays misleading plot axes if you use 'days since ...' as the units.
      ! if you simply use 'days' it works much better.

      ios = nf90_put_att(ncid, VarID, "units", "days")
      call nc_check(ios, routine, "units long_name")

   else if (dart_calendar == 'GREGORIAN') then
      ios = nf90_put_att(ncid, VarID, "calendar", "gregorian")
      call nc_check(ios, routine, "calendar long_name")

      ios = nf90_put_att(ncid, VarID, "units", "days since 1601-01-01 00:00:00")
      call nc_check(ios, routine, "units long_name")
   else
      write(string1,*) 'calendar type "'//trim(dart_calendar)// &
                       &'" unsupported by default write_model_time() routine'
      call error_handler(E_ERR, routine, string1, source, text2=string3)
   endif

   ! end define mode
   call nc_check( nf90_Enddef(ncid),routine, "Enddef" )
endif

! See if the existing time variable has a calendar and start date to consider

ios = nf90_get_att(ncid, VarID, 'calendar', var_calendar)
if (ios /= NF90_NOERR) var_calendar = 'NO_CALENDAR'
call to_upper(var_calendar)

if (dart_calendar /= var_calendar ) then
   ! allow NO_CALENDAR, NO CALENDAR, and NONE to be synonyms.
   ! replace var_calendar with what dart uses to simplify the tests below.
   if (dart_calendar == 'NO_CALENDAR' .and. (var_calendar == 'NONE'        .or. &
                                             var_calendar == 'NO CALENDAR')) then
      var_calendar = 'NO_CALENDAR'
   else
      write(string1,*)'inconsistent calendar types between DART program and input file.'
      write(string2,*)'DART initialized with: ', trim(dart_calendar), ' File uses: ', trim(var_calendar)
      call error_handler(E_ERR, routine, string1, source, text2=string2, text3=string3)
   endif
endif

! convert time to something that netCDF can store, fractional days
call get_time(dart_time, dart_seconds, dart_days)
model_time = real(dart_days,digits12) + real(dart_seconds,digits12)/86400.0_digits12

! need to know how long the time variable is and hammer the last time
ios = nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims)
call nc_check(ios, routine, 'inquire number of dimensions of variable "time"')

if (numdims == 0) then ! variable is a scalar
   ios = nf90_put_var(ncid, VarID, model_time)
   call nc_check(ios, routine, 'put_var scalar "model_time"')
   return
endif

if (numdims > 1) then
   write(string1,*)'Expecting the "time" variable to be a single dimension.'
   call error_handler(E_ERR, routine, string1, source, text2=string3)
endif

ios = nf90_inquire_dimension(ncid, dimIDs(1), len=ntimes)
call nc_check(ios, routine, 'inquire_dimension for time dimension')

if (dimIDs(1) == unlimitedDimID) time_is_unlimited = .true.

if (ntimes == 0 .and. time_is_unlimited) then
   ntimes = ntimes + 1
elseif (ntimes == 0) then
   write(string1,*)'"time" variable has length 0 but is not the unlimited dimension.'
   call error_handler(E_ERR, routine, string1, source, text2=string3)
endif

! write dart days and seconds files to netcdf file
ios = nf90_put_var(ncid, VarID, model_time, start=(/ ntimes /))
call nc_check(ios, routine, "put_var model_time")

end subroutine write_model_time

!--------------------------------------------------------------------
!> @}
end module dart_time_io_mod

