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

public :: nc_check, nc_add_global_attribute, nc_add_global_creation_time, &
          nc_redef, nc_enddef, nc_sync

interface nc_add_global_attribute
   module procedure nc_add_global_char_att
   module procedure nc_add_global_int_att
   module procedure nc_add_global_real_att
   module procedure nc_add_global_real_array_att
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

subroutine nc_check(istatus, subr_name, context, filename)

integer, intent (in)                   :: istatus
character(len=*), intent(in)           :: subr_name
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
  

if (istatus == nf90_noerr) return

! something wrong.  construct an error string and call the handler.
msgstring1 = nf90_strerror(istatus)

! context is optional, but is very useful if specified.
if (present(context)) then
  msgstring1 = trim(context) // ': ' // trim(msgstring1)
endif

! filename is optional, but is very useful if specified.
if (present(filename)) then
  msgstring1 = trim(msgstring1) // ', file "' // trim(filename) // '"'
endif

! this does not return 
call error_handler(E_ERR, subr_name, msgstring1, source, revision, revdate)
  

end subroutine nc_check

!------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine nc_add_global_char_att(ncid, name, val, context)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: name
character(len=*), intent(in) :: val
character(len=*), intent(in), optional :: context

integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, name, val)
call nc_check(ret, 'put_global_char_att', 'adding the global attribute: '//trim(name))

end subroutine nc_add_global_char_att

!--------------------------------------------------------------------

subroutine nc_add_global_int_att(ncid, name, val, context)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: name
integer,          intent(in) :: val
character(len=*), intent(in), optional :: context

integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, name, val)
call nc_check(ret, 'put_global_char_att', 'adding the global attribute: '//trim(name))

end subroutine nc_add_global_int_att

!--------------------------------------------------------------------

subroutine nc_add_global_real_att(ncid, name, val, context)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: name
real(r8),         intent(in) :: val
character(len=*), intent(in), optional :: context

integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, name, val)
call nc_check(ret, 'put_global_real_att', 'adding the global attribute: '//trim(name))

end subroutine nc_add_global_real_att

!--------------------------------------------------------------------

subroutine nc_add_global_real_array_att(ncid, name, val, context)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: name
real(r8),         intent(in) :: val(:)
character(len=*), intent(in), optional :: context

integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, name, val)
call nc_check(ret, 'put_global_real_array_att', 'adding the global attribute: '//trim(name))

end subroutine nc_add_global_real_array_att

!--------------------------------------------------------------------

subroutine nc_add_global_creation_time(ncid, context)
integer, intent(in) :: ncid
character(len=*), intent(in), optional :: context

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

character(len=128) :: str1

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_add_global_char_att(ncid, "creation_date", str1)

end subroutine nc_add_global_creation_time

!--------------------------------------------------------------------

subroutine nc_redef(ncid, context)
integer, intent(in) :: ncid
character(len=*), intent(in), optional :: context

integer :: ret

ret = nf90_Redef(ncid)
call nc_check(ret, 'nc_redef', 'nf90_Redef')

end subroutine nc_redef

!--------------------------------------------------------------------

subroutine nc_enddef(ncid, context)
integer, intent(in) :: ncid
character(len=*), intent(in), optional :: context

integer :: ret

ret = nf90_EndDef(ncid)
call nc_check(ret, 'nc_enddef', 'nf90_EndDef')

end subroutine nc_enddef

!--------------------------------------------------------------------

subroutine nc_sync(ncid, context)
integer, intent(in) :: ncid
character(len=*), intent(in), optional :: context

integer :: ret

ret = nf90_Sync(ncid)
call nc_check(ret, 'nc_sync', 'nf90_Sync')

end subroutine nc_sync

!--------------------------------------------------------------------
!--------------------------------------------------------------------

end module netcdf_utilities_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
