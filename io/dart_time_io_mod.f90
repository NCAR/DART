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
!>@todo some synergy with state_space_diag_mod.f90 routines ... nc_get_tindex, etc
!> @{

use types_mod,        only : r8, digits12
use time_manager_mod, only : time_type, set_time, get_time, print_time, &
                             set_calendar_type, set_date, get_calendar_string

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
!> and setting the calendar

function read_model_time(filename)

character(len=*), intent(in) :: filename
type(time_type)              :: read_model_time

integer :: ncid, ios, numdims, xtype, VarID, TimeDimID
integer :: ntimes, seconds, days
integer :: year, month, day, hour, minute, second
type(time_type) :: base_time

real(digits12), allocatable :: model_time(:)

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname,dimname
character(len=256) :: file_calendar, dart_calendar
character(len=256) :: unitstring

call nc_check( nf90_open(filename, NF90_NOWRITE, ncid), &
               'read_model_time',  'opening : "'//trim(filename)//'"')

ios = nf90_inq_varid(ncid, "time", VarID)
if (ios /= NF90_NOERR) then
   write(string1,*)'Expecting to be able to read the model time from a "time" variable'
   write(string2,*)'in "'//trim(filename)//'"'
   write(string3,*)'You may need to supply a model-specific "read_model_time()" to read the time.'
   call error_handler(E_ERR,'read_model_time', string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

ios = nf90_inquire_variable(ncid, VarID, xtype=xtype, dimids=dimIDs, ndims=numdims)
call nc_check(ios, 'read_model_time', 'inquire_variable "time" from "'//trim(filename)//'"')

if (numdims > 1) then
   write(string1,*)'Expecting the "time" variable to be a single dimension.'
   write(string2,*)'in "'//trim(filename)//'"'
   write(string3,*)'You may need to supply a model-specific "read_model_time()" to read the time.'
   call error_handler(E_ERR,'read_model_time', string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

! Since the time variable is known to have only 1 dimension, we know TimeDimID is the first one.
TimeDimID = dimids(1)

ios = nf90_inquire_dimension(ncid, TimeDimID, len=ntimes)
call nc_check(ios, 'read_model_time', 'inquire_dimension for ntimes "'//trim(filename) )

allocate( model_time(ntimes) )

ios = nf90_get_var(ncid, VarID, model_time)
call nc_check(ios, 'read_model_time','get_var time' )

!>@todo calendar madness 

call get_calendar_string(dart_calendar)

ios = nf90_get_att(ncid, VarID, 'calendar', file_calendar)
if (ios /= NF90_NOERR) file_calendar = 'no_calendar'
if (trim(file_calendar) == 'no calendar') file_calendar = 'no_calendar'
call to_upper(file_calendar)

if (trim(dart_calendar) /= trim(file_calendar) ) then
      write(*,*)trim(dart_calendar)
      write(*,*)trim(file_calendar)
      write(string1,*)'I dunno what to do.'
      call error_handler(E_ERR, 'read_model_time:', string1, source,revision, revdate)
endif

if ( trim(file_calendar) == 'NO_CALENDAR' ) then

   days    = floor(model_time(ntimes))
   seconds = (model_time(ntimes)-real(days,r8))*86400
   read_model_time = set_time(seconds, days)

   call print_time(read_model_time,'read_model_time')

   call nc_check( nf90_close(ncid) , 'read_model_time closing : ', filename)
   return
endif

!>@todo where does DART normally set the calendar now?
call set_calendar_type(trim(file_calendar))

ios = nf90_get_att(ncid, VarID, 'units', unitstring)
call nc_check(ios, 'read_model_time', 'get_att time units "'//trim(filename)//'"')

if (unitstring(1:10) == 'days since') then

   read(unitstring,'(11x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
   if (ios /= 0) then
      write(string1,*)'Unable to read dtime units. Error status was ',ios
      write(string2,*)'expected "days since YYYY-MM-DD HH:MM:SS"'
      write(string3,*)'was      "'//trim(unitstring)//'"'
      call error_handler(E_ERR, 'read_model_time:', string1, &
             source, revision, revdate, text2=string2, text3=string3)
   endif

   ! This is the start of their calendar
   base_time = set_date(year, month, day, hour, minute, second)

elseif (unitstring(1:13) == 'seconds since') then
         
   read(unitstring,'(14x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
   if (ios /= 0) then
      write(string1,*)'Unable to read dtime units. Error status was ',ios
      write(string2,*)'expected "seconds since YYYY-MM-DD HH:MM:SS"'
      write(string3,*)'was      "'//trim(unitstring)//'"'
      call error_handler(E_ERR, 'read_model_time:', string1, &
             source, revision, revdate, text2=string2, text3=string3)
   endif

   ! This is the start of their calendar
   base_time = set_date(year, month, day, hour, minute, second)

else

 !>@todo put something here

endif

!>@ todo  fix converting from our calendar to their calendar
! some_days    = floor(dstart)
! some_seconds = int(fraction(dstart) * (24*60*60))
! time_offset  = set_time(some_seconds, some_days)  ! seconds, days
! spinup_end   = base_time + time_offset


deallocate( model_time )

end function read_model_time


!--------------------------------------------------------------------
!> Write time to a netcdf file


subroutine write_model_time(ncid, dart_time)

integer,             intent(in) :: ncid
type(time_type),     intent(in) :: dart_time

integer  :: ios
integer  :: xtype, numdims, ntimes
integer  :: VarID, TimeDimID
integer  :: dart_days, dart_seconds

real(digits12) :: model_time

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname,dimname
character(len=NF90_MAX_NAME)          :: calendar

! convert time to something that netCDF can store
call get_time(dart_time, dart_seconds, dart_days)
model_time = real(dart_days,digits12) + real(dart_seconds,digits12)/86400.0_digits12

ios = nf90_inq_varid(ncid, "time", VarID)
if (ios /= NF90_NOERR) then

   call error_handler(E_MSG,'write_model_time','no time variable to exploit', &
              source, revision, revdate, text2='creating one')

   ! begin define mode
   ios = nf90_Redef(ncid)
   call nc_check(ios, "write_model_time", "redef")

   !>@todo NF90_UNLIMITED
   ios = nf90_def_var(ncid, name="time", xtype=nf90_int, varid=VarID)
   call nc_check(ios, "write_model_time", "time def_var")

   ! define time attributes consistent with CF convention
   ios = nf90_put_att(ncid, VarID, "long_name", "valid time of the model state")
   call nc_check(ios, "write_model_time", "time long_name")

   ios = nf90_put_att(ncid, VarID, "calendar", "no calendar")
   call nc_check(ios, "write_model_time", "calendar long_name")

   ios = nf90_put_att(ncid, VarID, "axis", "T")
   call nc_check(ios, "write_model_time", "axis long_name")

   ios = nf90_put_att(ncid, VarID, "cartesian_axis", "T")
   call nc_check(ios, "write_model_time", "cartesian_axis long_name")

   ios = nf90_put_att(ncid, VarID, "month_lengths", &
                  (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /))
   call nc_check(ios, "write_model_time", "month_lengths long_name")

   ios = nf90_put_att(ncid, VarID, "units", "days since 0000-01-01 00:00:00")
   call nc_check(ios, "write_model_time", "units long_name")

   ! end define mode
   call nc_check( nf90_Enddef(ncid),"write_model_time", "Enddef" )
endif

! See if the existing time dimension has a calendar and start date to consider

ios = nf90_get_att(ncid, VarID, 'calendar', calendar)
if (ios == NF90_NOERR) then
   call set_calendar_type(trim(calendar))
else
   call set_calendar_type('no calendar')
endif


! need to know how long the time variable is and hammer the last time
ios = nf90_inquire_variable(ncid, VarID, xtype=xtype, dimids=dimIDs, ndims=numdims)
call nc_check(ios, 'write_model_time', 'inquire_variable "time"')

if (numdims > 1) then
   write(string1,*)'Expecting the "time" variable to be a single dimension.'
   write(string2,*)'You may need to supply a model-specific "write_model_time()".'
   call error_handler(E_ERR,'write_model_time', string1, &
              source, revision, revdate, text2=string2)
endif

! Since the time variable is known to have only 1 dimension, we know TimeDimID is the first one.
TimeDimID = dimids(1)

ios = nf90_inquire_dimension(ncid, TimeDimID, len=ntimes)
call nc_check(ios, 'write_model_time', 'inquire_dimension for ntimes ')

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
