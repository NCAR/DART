! Netcdf reading and writing dart model time
! Temporary module for dart time.
!> @todo should this go in state_vector_io_mod or io_filename_mod?
module dart_time_io_mod

use types_mod,        only : r8, i8
use time_manager_mod, only : time_type, set_time, get_time
use utilities_mod,    only : nc_check

use typeSizes
use netcdf

private

public :: read_model_time, write_model_time

contains

!--------------------------------------------------------------------
!> read the dart time from the input file
!--------------------------------------------------------------------
function read_model_time(filename)

character(len=1024), intent(in) :: filename
type(time_type) :: read_model_time

integer :: ret !< netcdf return code
integer :: ncid, dart_secsVarID, dart_daysVarID
integer :: seconds, days

! open netcdf file
call nc_check( nf90_open(filename, NF90_NOWRITE, ncid), &
               'read_model_time opening : ', filename )


! grab dart_days from file
call nc_check( nf90_inq_varid(ncid, "dart_days", dart_daysVarID), &
               'read_model_time', 'inq_varid dart_days' )

call nc_check( nf90_get_var(ncid, dart_daysVarID, days), &
               'read_model_time','get_var dart_days' )


! grab dart_seconds from file
call nc_check( nf90_inq_varid(ncid, "dart_seconds", dart_secsVarID), &
               'read_model_time', 'inq_varid dart_days' )

call nc_check( nf90_get_var(ncid, dart_secsVarID, seconds), &
               'read_model_time','get_var dart_seconds' )

read_model_time = set_time(seconds, days)

call nc_check( nf90_close(ncid) , 'read_model_time closing : ', filename)

end function read_model_time

!--------------------------------------------------------------------
subroutine write_model_time(ncid, dart_time)

integer,             intent(in) :: ncid
type(time_type),     intent(in) :: dart_time
integer                         :: ret !< netcdf return code
integer                         :: dart_daysVarID, dart_secondsVarID
integer                         :: dart_days, dart_seconds

! begin define mode
call nc_check(nf90_Redef(ncid),"write_model_time", "redef")

! define days if it does not exist
ret = nf90_inq_varid(ncid, "dart_days", dart_daysVarID)
if (ret /= NF90_NOERR) then
   call nc_check( nf90_def_var(ncid, name="dart_days", &
                  xtype=nf90_int, varid=dart_daysVarID) , &
                  "write_model_time", "dart_days def_var")

   ! define days attributed
   call nc_check( nf90_put_att(ncid, dart_daysVarID, "long_name", "days"), &
                  "write_model_time", "dart_days long_name")

   call nc_check( nf90_put_att(ncid, dart_daysVarID, "calendar", "no calendar"), &
                  "write_model_time", "dart_days calendar")

   call nc_check( nf90_put_att(ncid, dart_daysVarID, "units", "days since 0000-00-00 00:00:00"), &
                  "write_model_time", "dart_days units")
endif

! define seconds if it does not exist
ret = nf90_inq_varid(ncid, "dart_seconds", dart_secondsVarID)
if (ret /= NF90_NOERR) then
   call nc_check( nf90_def_var(ncid, name="dart_seconds", &
                  xtype=nf90_int, varid=dart_secondsVarID) , &
                  "write_model_time", "dart_secondsdef_var")

   ! define seconds attributed
   call nc_check( nf90_put_att(ncid, dart_secondsVarID, "long_name", "seconds"), &
                  "write_model_time", "dart_seconds long_name")

   call nc_check( nf90_put_att(ncid, dart_secondsVarID, "calendar", "no calendar"), &
                  "write_model_time", "dart_seconds calendar")

   call nc_check( nf90_put_att(ncid, dart_secondsVarID, "units", "seconds since midnight"), &
                  "write_model_time", "dart_seconds units")
endif
  
! end define mode
call nc_check( nf90_Enddef(ncid),"write_model_time", "Enddef" )

! get the dart time
call get_time(dart_time, dart_seconds, dart_days)

! write dart days and seconds files to netcdf file
call nc_check( nf90_put_var(ncid, dart_daysVarID, dart_days ), &
               "write_model_time", "dart_days put_var")

call nc_check( nf90_put_var(ncid, dart_secondsVarID, dart_seconds ), &
               "write_model_time", "dart_seconds put_var")

end subroutine write_model_time

!--------------------------------------------------------------------

end module dart_time_io_mod

