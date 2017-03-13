! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module netcdf_utilities_mod

!>@todo FIXME: NEEDS MUCH MORE WORK.
!>
!> start of a module to help streamline the code that
!> calls netcdf routines.  usually adds error checking.
!> 
!> also isolates as much access to the netcdf libs to routines
!> inside this module.
!>
!> routines in this file are prefixed with nc_

use types_mod, only : r4, r8, digits12, i4, i8, PI, MISSING_R8, MISSING_I
use utilities_mod, only : error_handler, E_DBG, E_MSG, E_ALLMSG, E_WARN, E_ERR
           
use netcdf

implicit none
private

public :: nc_check, nc_add_global_attribute, nc_add_global_creation_time

interface nc_add_global_attribute
   module procedure nc_add_global_string_attr
   module procedure nc_add_global_real_attr
end interface

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: msgstring1, msgstring2, msgstring3

! do we need one of these?
!namelist /netcdf_utilities_nml/ 

contains

!------------------------------------------------------------------
!> check return code from previous call. on error, print and stop.
!> if you want to continue after an error don't use this call. 

subroutine nc_check(istatus, subr_name, context)

integer, intent (in)                   :: istatus
character(len=*), intent(in)           :: subr_name
character(len=*), intent(in), optional :: context
  

if (istatus == nf90_noerr) return

! something wrong.  construct an error string and call the handler.

! context is optional, but is very useful if specified.
! if context + error code > 129, the assignment will truncate.
if (present(context)) then
  error_msg = trim(context) // ': ' // trim(nf90_strerror(istatus))
else
  error_msg = nf90_strerror(istatus)
endif

! this does not return 
call error_handler(E_ERR, subr_name, error_msg, source, revision, revdate)
  

end subroutine nc_check

!------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine nc_add_global_char_att(ncid, name, val)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: name
character(len=*), intent(in) :: val

integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, name, val)
call nc_check(ret, 'put_global_char_att', 'adding the global attribute: '//trim(name))

end subroutine nc_add_global_char_att

!--------------------------------------------------------------------

subroutine nc_add_global_real_att(ncid, name, val)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: name
real(r8),         intent(in) :: val

integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, name, val)
call nc_check(ret, 'put_global_real_att', 'adding the global attribute: '//trim(name))

end subroutine nc_add_global_real_att

!--------------------------------------------------------------------

subroutine nc_add_global_creation_time(ncid)
integer, intent(in) :: ncid

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

character(len=128) :: str1

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_add_global_char_att(ncid, "creation_date", str1)

end subroutine nc_put_global_creation_time

!--------------------------------------------------------------------
!--------------------------------------------------------------------

end module netcdf_utilities_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
