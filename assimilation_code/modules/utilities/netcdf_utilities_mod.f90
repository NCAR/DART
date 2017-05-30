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
!>
!> the intent is that these are simple routines that aren't trying to
!> handle all possible options.  if you have complicated needs, write your
!> own specialized routine, either here or in the calling code.
!> for example, the 'put_var' routine here doesn't have start or count, 
!> intentionally, for array entries.


use types_mod, only : r4, r8, digits12, i4, i8, PI, MISSING_R8, MISSING_I
use utilities_mod, only : error_handler, E_DBG, E_MSG, E_ALLMSG, E_WARN, E_ERR
           
use netcdf

implicit none
private

public :: nc_check,                     &
          nc_add_global_attribute,      &
          nc_add_attribute_to_variable, &
          nc_define_dimension,          &
          nc_define_integer_variable,   &
          nc_define_real_variable,      &
          nc_define_double_variable,    &
          nc_put_variable,              &
          nc_get_variable,              &
          nc_add_global_creation_time,  &
          nc_redef,                     &
          nc_enddef,                    &
          nc_sync


interface nc_add_global_attribute
   module procedure nc_add_global_char_att
   module procedure nc_add_global_int_att
   module procedure nc_add_global_real_att
   module procedure nc_add_global_real_array_att
end interface

interface nc_add_attribute_to_variable
   module procedure nc_add_char_att_to_var
   module procedure nc_add_int_array_att_to_var
   module procedure nc_add_int_att_to_var
   module procedure nc_add_real_att_to_var
   module procedure nc_add_real_array_att_to_var
end interface

interface nc_define_integer_variable
   module procedure nc_define_var_int_1d
   module procedure nc_define_var_int_Nd
end interface

interface nc_define_real_variable
   module procedure nc_define_var_real_1d
   module procedure nc_define_var_real_Nd
end interface

interface nc_define_double_variable
   module procedure nc_define_var_double_1d
   module procedure nc_define_var_double_Nd
end interface

interface nc_put_variable
   module procedure nc_put_single_int_1d
   module procedure nc_put_int_1d
   module procedure nc_put_single_real_1d
   module procedure nc_put_real_1d
   module procedure nc_put_int_2d
   module procedure nc_put_real_2d
   module procedure nc_put_int_3d
   module procedure nc_put_real_3d
end interface

!>@todo FIXME do we need single value versions
!> of these as well?
interface nc_get_variable
   module procedure nc_get_int_1d
   module procedure nc_get_real_1d
   module procedure nc_get_int_2d
   module procedure nc_get_real_2d
   module procedure nc_get_int_3d
   module procedure nc_get_real_3d
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

subroutine nc_check(istatus, subr_name, context, context2, filename)

integer,          intent(in)           :: istatus
character(len=*), intent(in)           :: subr_name
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: context2
character(len=*), intent(in), optional :: filename
  

if (istatus == nf90_noerr) return

! something wrong.  construct an error string and call the handler.
msgstring1 = nf90_strerror(istatus)

! context is optional, but is very useful if specified.
if (present(context)) then
  msgstring1 = trim(context) // ': ' // trim(msgstring1)
endif

! this does not return 
call error_handler(E_ERR, subr_name, msgstring1, source, revision, revdate, &
                   text2=context2, text3=filename)
  

end subroutine nc_check

!------------------------------------------------------------------
!--------------------------------------------------------------------
! global attribute section

subroutine nc_add_global_char_att(ncid, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: attname
character(len=*), intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_global_char_att'
integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'adding the global attribute: '//trim(attname), context, filename)

end subroutine nc_add_global_char_att

!--------------------------------------------------------------------

subroutine nc_add_global_int_att(ncid, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: attname
integer,          intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_global_int_att'
integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'adding the global attribute: '//trim(attname), context, filename)

end subroutine nc_add_global_int_att

!--------------------------------------------------------------------

subroutine nc_add_global_real_att(ncid, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: attname
real(r8),         intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_global_real_att'
integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'adding the global attribute: '//trim(attname), context, filename)

end subroutine nc_add_global_real_att

!--------------------------------------------------------------------

subroutine nc_add_global_real_array_att(ncid, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: attname
real(r8),         intent(in) :: val(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_global_real_array_att'
integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'adding the global attribute: '//trim(attname), context, filename)

end subroutine nc_add_global_real_array_att

!------------------------------------------------------------------
!--------------------------------------------------------------------
! adding attribute to a variable section

subroutine nc_add_char_att_to_var(ncid, varname, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
character(len=*), intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_char_att_to_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename)

end subroutine nc_add_char_att_to_var

!--------------------------------------------------------------------

subroutine nc_add_int_att_to_var(ncid, varname, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
integer,          intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_int_att_to_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename)

end subroutine nc_add_int_att_to_var

!--------------------------------------------------------------------

subroutine nc_add_int_array_att_to_var(ncid, varname, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
integer,          intent(in) :: val(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_int_array_att_to_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename)

end subroutine nc_add_int_array_att_to_var

!--------------------------------------------------------------------

subroutine nc_add_real_att_to_var(ncid, varname, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
real(r8),         intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_real_att_to_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename)

end subroutine nc_add_real_att_to_var

!--------------------------------------------------------------------

subroutine nc_add_real_array_att_to_var(ncid, varname, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
real(r8),         intent(in) :: val(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_real_array_att_to_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename)

end subroutine nc_add_real_array_att_to_var

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! defining dimensions section

subroutine nc_define_dimension(ncid, dimname, dimlen, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: dimname
integer,          intent(in) :: dimlen
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_dimension'
integer :: ret, dimid

ret = nf90_def_dim(ncid, dimname, dimlen, dimid)
call nc_check(ret, routine, 'define dimension '//trim(dimname), context, filename)

end subroutine nc_define_dimension

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! defining variables section

subroutine nc_define_var_int_1d(ncid, varname, dimname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_int_1d'
integer :: ret, dimid, varid

ret = nf90_inq_dimid(ncid, dimname, dimid)
call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimname), context, filename)

ret = nf90_def_var(ncid, varname, nf90_int, dimid, varid)
call nc_check(ret, routine, 'define integer variable '//trim(varname), context, filename)

end subroutine nc_define_var_int_1d

!--------------------------------------------------------------------

subroutine nc_define_var_int_Nd(ncid, varname, dimnames, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimnames(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_int_Nd'
integer :: ret, dimid1, dimid2, dimid3, varid

if (size(dimnames) >= 1) then
   ret = nf90_inq_dimid(ncid, dimnames(1), dimid1)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(1)), context, filename)
endif

if (size(dimnames) >= 2) then
   ret = nf90_inq_dimid(ncid, dimnames(2), dimid2)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(2)), context, filename)
endif

if (size(dimnames) >= 3) then
   ret = nf90_inq_dimid(ncid, dimnames(3), dimid3)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(3)), context, filename)
endif

if (size(dimnames) >= 4) then
   call error_handler(E_ERR, routine, 'only 1d, 2d and 3d integer variables supported', &
                      source, revision, revdate, text2='variable '//trim(varname))
endif

if (size(dimnames) == 1) then
   ret = nf90_def_var(ncid, varname, nf90_int, dimid1, varid=varid)
else if (size(dimnames) == 2) then
   ret = nf90_def_var(ncid, varname, nf90_int, dimids=(/ dimid1, dimid2 /), varid=varid)
else if (size(dimnames) == 3) then
   ret = nf90_def_var(ncid, varname, nf90_int, dimids=(/ dimid1, dimid2, dimid3 /), varid=varid)
endif

call nc_check(ret, routine, 'define integer variable '//trim(varname), context, filename)

end subroutine nc_define_var_int_Nd

!--------------------------------------------------------------------

subroutine nc_define_var_real_1d(ncid, varname, dimname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_real_1d'
integer :: ret, dimid, varid

ret = nf90_inq_dimid(ncid, dimname, dimid)
call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimname), context, filename)

ret = nf90_def_var(ncid, varname, nf90_real, dimid, varid)
call nc_check(ret, routine, 'define real variable '//trim(varname), context, filename)

end subroutine nc_define_var_real_1d

!--------------------------------------------------------------------

subroutine nc_define_var_real_Nd(ncid, varname, dimnames, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimnames(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_real_Nd'
integer :: ret, dimid1, dimid2, dimid3, varid

if (size(dimnames) >= 1) then
   ret = nf90_inq_dimid(ncid, dimnames(1), dimid1)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(1)), context, filename)
endif

if (size(dimnames) >= 2) then
   ret = nf90_inq_dimid(ncid, dimnames(2), dimid2)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(2)), context, filename)
endif

if (size(dimnames) >= 3) then
   ret = nf90_inq_dimid(ncid, dimnames(3), dimid3)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(3)), context, filename)
endif

if (size(dimnames) >= 4) then
   call error_handler(E_ERR, routine, 'only 1d, 2d and 3d real variables supported', &
                      source, revision, revdate, text2='variable '//trim(varname))
endif

if (size(dimnames) == 1) then
   ret = nf90_def_var(ncid, varname, nf90_real, dimid1, varid=varid)
else if (size(dimnames) == 2) then
   ret = nf90_def_var(ncid, varname, nf90_real, dimids=(/ dimid1, dimid2 /), varid=varid)
else if (size(dimnames) == 3) then
   ret = nf90_def_var(ncid, varname, nf90_real, dimids=(/ dimid1, dimid2, dimid3 /), varid=varid)
endif

call nc_check(ret, routine, 'define real variable '//trim(varname), context, filename)

end subroutine nc_define_var_real_Nd

!--------------------------------------------------------------------

subroutine nc_define_var_double_1d(ncid, varname, dimname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_double_1d'
integer :: ret, dimid, varid

ret = nf90_inq_dimid(ncid, dimname, dimid)
call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimname), context, filename)

ret = nf90_def_var(ncid, varname, nf90_double, dimid, varid)
call nc_check(ret, routine, 'define double variable '//trim(varname), context, filename)

end subroutine nc_define_var_double_1d

!--------------------------------------------------------------------

subroutine nc_define_var_double_Nd(ncid, varname, dimnames, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimnames(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_double_Nd'
integer :: ret, dimid1, dimid2, dimid3, varid

if (size(dimnames) >= 1) then
   ret = nf90_inq_dimid(ncid, dimnames(1), dimid1)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(1)), context, filename)
endif

if (size(dimnames) >= 2) then
   ret = nf90_inq_dimid(ncid, dimnames(2), dimid2)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(2)), context, filename)
endif

if (size(dimnames) >= 3) then
   ret = nf90_inq_dimid(ncid, dimnames(3), dimid3)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(3)), context, filename)
endif

if (size(dimnames) >= 4) then
   call error_handler(E_ERR, routine, 'only 1d, 2d and 3d double variables supported', &
                      source, revision, revdate, text2='variable '//trim(varname))
endif

if (size(dimnames) == 1) then
   ret = nf90_def_var(ncid, varname, nf90_double, dimid1, varid=varid)
else if (size(dimnames) == 2) then
   ret = nf90_def_var(ncid, varname, nf90_double, dimids=(/ dimid1, dimid2 /), varid=varid)
else if (size(dimnames) == 3) then
   ret = nf90_def_var(ncid, varname, nf90_double, dimids=(/ dimid1, dimid2, dimid3 /), varid=varid)
endif

call nc_check(ret, routine, 'define double variable '//trim(varname), context, filename)

end subroutine nc_define_var_double_Nd

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! put values into variables

subroutine nc_put_single_int_1d(ncid, varname, varindex, varval, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: varindex
integer,          intent(in) :: varval
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_int_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_put_var(ncid, varid, varval, start=(/varindex/))
call nc_check(ret, routine, 'put value for '//trim(varname), context, filename)

end subroutine nc_put_single_int_1d

!--------------------------------------------------------------------

subroutine nc_put_int_1d(ncid, varname, varvals, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_int_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_put_var(ncid, varid, varvals)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename)

end subroutine nc_put_int_1d

!--------------------------------------------------------------------

subroutine nc_put_single_real_1d(ncid, varname, varindex, varval, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: varindex
real(r8),         intent(in) :: varval
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_real_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_put_var(ncid, varid, varval, start=(/varindex/))
call nc_check(ret, routine, 'put value for '//trim(varname), context, filename)

end subroutine nc_put_single_real_1d

!--------------------------------------------------------------------

subroutine nc_put_real_1d(ncid, varname, varvals, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(r8),         intent(in) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_real_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_put_var(ncid, varid, varvals)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename)

end subroutine nc_put_real_1d

!--------------------------------------------------------------------

subroutine nc_put_int_2d(ncid, varname, varvals, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_int_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_put_var(ncid, varid, varvals)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename)

end subroutine nc_put_int_2d

!--------------------------------------------------------------------

subroutine nc_put_real_2d(ncid, varname, varvals, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(r8),         intent(in) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_real_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_put_var(ncid, varid, varvals)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename)

end subroutine nc_put_real_2d

!--------------------------------------------------------------------

subroutine nc_put_int_3d(ncid, varname, varvals, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_int_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_put_var(ncid, varid, varvals)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename)

end subroutine nc_put_int_3d

!--------------------------------------------------------------------

subroutine nc_put_real_3d(ncid, varname, varvals, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(r8),         intent(in) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_real_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_put_var(ncid, varid, varvals)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename)

end subroutine nc_put_real_3d

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! get values from variables

subroutine nc_get_int_1d(ncid, varname, varvals, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_int_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_get_var(ncid, varid, varvals)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename)

end subroutine nc_get_int_1d

!--------------------------------------------------------------------

subroutine nc_get_real_1d(ncid, varname, varvals, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(r8),         intent(out) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_real_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_get_var(ncid, varid, varvals)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename)

end subroutine nc_get_real_1d

!--------------------------------------------------------------------

subroutine nc_get_int_2d(ncid, varname, varvals, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_int_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_get_var(ncid, varid, varvals)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename)

end subroutine nc_get_int_2d

!--------------------------------------------------------------------

subroutine nc_get_real_2d(ncid, varname, varvals, context, filename)
integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(r8),         intent(out) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_real_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_get_var(ncid, varid, varvals)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename)

end subroutine nc_get_real_2d

!--------------------------------------------------------------------

subroutine nc_get_int_3d(ncid, varname, varvals, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_int_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_get_var(ncid, varid, varvals)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename)

end subroutine nc_get_int_3d

!--------------------------------------------------------------------

subroutine nc_get_real_3d(ncid, varname, varvals, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(r8),         intent(out) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_real_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename)

ret = nf90_get_var(ncid, varid, varvals)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename)

end subroutine nc_get_real_3d

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! misc section: file operations, standard timestamp routine

subroutine nc_add_global_creation_time(ncid, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

character(len=128) :: str1

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_add_global_char_att(ncid, "creation_date", str1, context, filename)

end subroutine nc_add_global_creation_time

!--------------------------------------------------------------------

subroutine nc_redef(ncid, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_redef'
integer :: ret

ret = nf90_Redef(ncid)
call nc_check(ret, routine, 'begin file define mode', context, filename)

end subroutine nc_redef

!--------------------------------------------------------------------

subroutine nc_enddef(ncid, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_enddef'
integer :: ret

ret = nf90_EndDef(ncid)
call nc_check(ret, routine, 'end file define mode', context, filename)

end subroutine nc_enddef

!--------------------------------------------------------------------

subroutine nc_sync(ncid, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_sync'
integer :: ret

ret = nf90_Sync(ncid)
call nc_check(ret, routine, 'file sync', context, filename)

end subroutine nc_sync

!--------------------------------------------------------------------
!--------------------------------------------------------------------

end module netcdf_utilities_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
