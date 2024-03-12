! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module netcdf_utilities_mod

!> a module to streamline the code that calls netcdf routines. 
!>
!> calling code should NOT have to use the netcdf module, nor see
!> any netcdf ID values, types, lengths, etc.  this code does error 
!> checking and reports the offending call, filename, and user-supplied 
!> context (often the 'routine' variable with the calling subroutine name).
!> 
!> isolates access to the netcdf libs to routines inside this module.
!>
!> routines in this file are prefixed with nc_ and attempt to say
!> in english what they do.  the intent is someone who does not know
!> anything about netcdf can still use these.
!>
!> The 'put_var' and 'get_var' routines take optional start, count, 
!> stride, and map variables.
!>


use types_mod, only : r4, r8, digits12, i2, i4, i8, PI, MISSING_R8, MISSING_I
use utilities_mod, only : error_handler, E_DBG, E_MSG, E_ALLMSG, E_WARN, E_ERR
           
use netcdf

implicit none
private

public :: nc_check,                       &
          nc_add_global_attribute,        &
          nc_get_global_attribute,        &
          nc_add_attribute_to_variable,   &
          nc_get_attribute_from_variable, &
          nc_define_dimension,            &
          nc_define_unlimited_dimension,  &
          nc_get_dimension_size,          &
          nc_define_character_variable,   &
          nc_define_integer_variable,     &
          nc_define_real_variable,        &
          nc_define_double_variable,      &
          nc_define_integer_scalar,       &
          nc_define_real_scalar,          &
          nc_define_double_scalar,        &
          nc_global_attribute_exists,     &
          nc_variable_attribute_exists,   &
          nc_dimension_exists,            &
          nc_variable_exists,             &
          nc_put_variable,                &
          nc_get_variable,                &
          nc_get_variable_info,           &
          nc_add_global_creation_time,    &
          nc_get_variable_num_dimensions, &
          nc_get_variable_dimension_names, &
          nc_get_variable_size,           &
          nc_open_file_readonly,          &
          nc_open_file_readwrite,         &
          nc_create_file,                 &
          nc_close_file,                  &
          nc_begin_define_mode,           &
          nc_end_define_mode,             &
          nc_synchronize_file,            &
          NF90_MAX_NAME, NF90_MAX_VAR_DIMS, NF90_FILL_REAL


! note here that you only need to distinguish between
! r4 (float) and r8 (double) when defining or adding
! a new variable or attribute.  the get and query routines 
! will coerce the values to the destination precision correctly.

interface nc_add_global_attribute
   module procedure nc_add_global_char_att
   module procedure nc_add_global_int_att
   module procedure nc_add_global_float_att
   module procedure nc_add_global_double_att
   module procedure nc_add_global_float_array_att
   module procedure nc_add_global_double_array_att
end interface

interface nc_get_global_attribute
   module procedure nc_get_global_char_att
   module procedure nc_get_global_int_att
   module procedure nc_get_global_float_att
   module procedure nc_get_global_double_att
   module procedure nc_get_global_float_array_att
   module procedure nc_get_global_double_array_att
end interface

interface nc_add_attribute_to_variable
   module procedure nc_add_char_att_to_var
   module procedure nc_add_int_array_att_to_var
   module procedure nc_add_int_att_to_var
   module procedure nc_add_float_att_to_var
   module procedure nc_add_double_att_to_var
   module procedure nc_add_float_array_att_to_var
   module procedure nc_add_double_array_att_to_var
end interface

interface nc_get_attribute_from_variable
   module procedure nc_get_char_att_from_var
   module procedure nc_get_int_array_att_from_var
   module procedure nc_get_int_att_from_var
   module procedure nc_get_float_att_from_var
   module procedure nc_get_double_att_from_var
   module procedure nc_get_float_array_att_from_var
   module procedure nc_get_double_array_att_from_var
end interface

interface nc_define_character_variable
   module procedure nc_define_var_char_1d
   module procedure nc_define_var_char_Nd
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
   module procedure nc_put_char_1d
   module procedure nc_put_char_2d
   module procedure nc_put_int_0d
   module procedure nc_put_int_1d
   module procedure nc_put_int_2d
   module procedure nc_put_int_3d
   module procedure nc_put_int_4d
   module procedure nc_put_real_0d
   module procedure nc_put_real_1d
   module procedure nc_put_real_2d
   module procedure nc_put_real_3d
   module procedure nc_put_real_4d
   module procedure nc_put_double_0d
   module procedure nc_put_double_1d
   module procedure nc_put_double_2d
   module procedure nc_put_double_3d
   module procedure nc_put_double_4d
end interface

interface nc_get_variable
   module procedure nc_get_short_0d
   module procedure nc_get_short_1d
   module procedure nc_get_short_2d
   module procedure nc_get_short_3d
   module procedure nc_get_int_0d
   module procedure nc_get_int_1d
   module procedure nc_get_int_2d
   module procedure nc_get_int_3d
   module procedure nc_get_int_4d
   module procedure nc_get_real_0d
   module procedure nc_get_real_1d
   module procedure nc_get_real_2d
   module procedure nc_get_real_3d
   module procedure nc_get_real_4d
   module procedure nc_get_double_0d
   module procedure nc_get_double_1d
   module procedure nc_get_double_2d
   module procedure nc_get_double_3d
   module procedure nc_get_double_4d
end interface

interface nc_get_variable_size
   module procedure nc_get_variable_size_1d
   module procedure nc_get_variable_size_Nd
end interface

character(len=*), parameter :: source = 'netcdf_utilities_mod.f90'

character(len=512) :: msgstring1

!> make a derived type that is (ncid, filename) 
!> store filename on file open, delete it on file close. cache the
!> last N filenames - look them up on error and stop
!> having to keep the filename around.

! NOTE this is the max number of concurrently 
! open netcdf files on a single task.
integer, parameter :: MAX_NCFILES = 50
integer, parameter :: FH_EMPTY = -1

type ncinfo_type
  integer :: file_handle = FH_EMPTY
  character(len=256) :: file_name = ''
end type

! for now hardcode max size.  could make it allocatable
! and extend it if you need more open slots

type(ncinfo_type) :: ncinfo(MAX_NCFILES)

! do we need one of these?
!namelist /netcdf_utilities_nml/ 

contains

!------------------------------------------------------------------
!> check return code from previous call. on error, print and stop.
!> if you want to continue after an error don't use this call. 

subroutine nc_check(istatus, subr_name, context, context2, filename, ncid)

integer,          intent(in)           :: istatus
character(len=*), intent(in)           :: subr_name
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: context2
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: ncid
  
character(len=256) :: saved_filename

if (istatus == nf90_noerr) return

! something wrong.  construct an error string and call the handler.
msgstring1 = nf90_strerror(istatus)

! context is optional, but is very useful if specified.
if (present(context)) then
  msgstring1 = trim(context) // ': ' // trim(msgstring1)
endif

! you can give this routine a file handle and it will try
! to extract the filename from it.
if (.not. present(filename) .and. present(ncid)) then
   call find_name_from_fh(ncid, saved_filename)
else if (present(filename)) then
   saved_filename = filename
else
   saved_filename = ''
endif

! this does not return 
call error_handler(E_ERR, subr_name, msgstring1, source, &
                   text2=context2, text3=saved_filename)


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
call nc_check(ret, routine, 'adding the global attribute: '//trim(attname), context, filename, ncid)

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
call nc_check(ret, routine, 'adding the global attribute: '//trim(attname), context, filename, ncid)

end subroutine nc_add_global_int_att

!--------------------------------------------------------------------

subroutine nc_add_global_float_att(ncid, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: attname
real(r4),         intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_global_float_att'
integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'adding the global attribute: '//trim(attname), context, filename, ncid)

end subroutine nc_add_global_float_att

!--------------------------------------------------------------------

subroutine nc_add_global_double_att(ncid, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: attname
real(digits12),   intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_global_double_att'
integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'adding the global attribute: '//trim(attname), context, filename, ncid)

end subroutine nc_add_global_double_att

!--------------------------------------------------------------------

subroutine nc_add_global_float_array_att(ncid, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: attname
real(r4),         intent(in) :: val(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_global_float_array_att'
integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'adding the global attribute: '//trim(attname), context, filename, ncid)

end subroutine nc_add_global_float_array_att

!--------------------------------------------------------------------

subroutine nc_add_global_double_array_att(ncid, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: attname
real(digits12),   intent(in) :: val(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_global_double_array_att'
integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'adding the global attribute: '//trim(attname), context, filename, ncid)

end subroutine nc_add_global_double_array_att

!------------------------------------------------------------------

subroutine nc_get_global_char_att(ncid, attname, val, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: attname
character(len=*), intent(out) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_global_char_att'
integer :: ret

ret = nf90_get_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'getting the global attribute: '//trim(attname), context, filename, ncid)

end subroutine nc_get_global_char_att

!--------------------------------------------------------------------

subroutine nc_get_global_int_att(ncid, attname, val, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: attname
integer,          intent(out) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_global_int_att'
integer :: ret

ret = nf90_get_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'getting the global attribute: '//trim(attname), context, filename, ncid)

end subroutine nc_get_global_int_att

!--------------------------------------------------------------------

subroutine nc_get_global_float_att(ncid, attname, val, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: attname
real(r4),         intent(out) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_global_float_att'
integer :: ret

ret = nf90_get_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'getting the global attribute: '//trim(attname), context, filename, ncid)

end subroutine nc_get_global_float_att

!--------------------------------------------------------------------

subroutine nc_get_global_double_att(ncid, attname, val, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: attname
real(digits12),         intent(out) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_global_double_att'
integer :: ret

ret = nf90_get_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'getting the global attribute: '//trim(attname), context, filename, ncid)

end subroutine nc_get_global_double_att

!--------------------------------------------------------------------

subroutine nc_get_global_float_array_att(ncid, attname, val, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: attname
real(r4),         intent(out) :: val(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_global_float_array_att'
integer :: ret

ret = nf90_get_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'getting the global attribute: '//trim(attname), context, filename, ncid)

end subroutine nc_get_global_float_array_att

!--------------------------------------------------------------------

subroutine nc_get_global_double_array_att(ncid, attname, val, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: attname
real(digits12),   intent(out) :: val(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_global_double_array_att'
integer :: ret

ret = nf90_get_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'getting the global attribute: '//trim(attname), context, filename, ncid)

end subroutine nc_get_global_double_array_att

!--------------------------------------------------------------------
! attributes on specific variables section

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
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

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
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

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
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_add_int_array_att_to_var

!--------------------------------------------------------------------

subroutine nc_add_float_att_to_var(ncid, varname, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
real(r4),         intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_float_att_to_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_add_float_att_to_var

!--------------------------------------------------------------------

subroutine nc_add_double_att_to_var(ncid, varname, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
real(digits12),   intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_double_att_to_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_add_double_att_to_var

!--------------------------------------------------------------------

subroutine nc_add_float_array_att_to_var(ncid, varname, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
real(r4),         intent(in) :: val(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_float_array_att_to_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_add_float_array_att_to_var

!--------------------------------------------------------------------

subroutine nc_add_double_array_att_to_var(ncid, varname, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
real(digits12),   intent(in) :: val(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_double_array_att_to_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_add_double_array_att_to_var

!--------------------------------------------------------------------

subroutine nc_get_char_att_from_var(ncid, varname, attname, val, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
character(len=*), intent(in)  :: attname
character(len=*), intent(out) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_char_att_from_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_get_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'getting the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_get_char_att_from_var

!--------------------------------------------------------------------

subroutine nc_get_int_att_from_var(ncid, varname, attname, val, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
character(len=*), intent(in)  :: attname
integer,          intent(out) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_int_att_from_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_get_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'getting the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_get_int_att_from_var

!--------------------------------------------------------------------

subroutine nc_get_int_array_att_from_var(ncid, varname, attname, val, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
character(len=*), intent(in)  :: attname
integer,          intent(out) :: val(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_int_array_att_from_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_get_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'getting the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_get_int_array_att_from_var

!--------------------------------------------------------------------

subroutine nc_get_float_att_from_var(ncid, varname, attname, val, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
character(len=*), intent(in)  :: attname
real(r4),         intent(out) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_float_att_from_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_get_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'getting the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_get_float_att_from_var

!--------------------------------------------------------------------

subroutine nc_get_double_att_from_var(ncid, varname, attname, val, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
character(len=*), intent(in)  :: attname
real(digits12),   intent(out) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_double_att_from_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_get_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'getting the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_get_double_att_from_var

!--------------------------------------------------------------------

subroutine nc_get_float_array_att_from_var(ncid, varname, attname, val, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
character(len=*), intent(in)  :: attname
real(r4),         intent(out) :: val(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_float_array_att_from_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_get_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'getting the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_get_float_array_att_from_var

!--------------------------------------------------------------------

subroutine nc_get_double_array_att_from_var(ncid, varname, attname, val, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
character(len=*), intent(in)  :: attname
real(digits12),   intent(out) :: val(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_double_array_att_from_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_get_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'getting the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_get_double_array_att_from_var

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! dimensions section

subroutine nc_define_dimension(ncid, dimname, dimlen, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: dimname
integer,          intent(in) :: dimlen
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_dimension'
integer :: ret, dimid

ret = nf90_def_dim(ncid, dimname, dimlen, dimid)
call nc_check(ret, routine, 'define dimension '//trim(dimname), context, filename, ncid)

end subroutine nc_define_dimension

!--------------------------------------------------------------------

subroutine nc_define_unlimited_dimension(ncid, dimname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: dimname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_unlimited_dimension'
integer :: ret, dimid

ret = nf90_def_dim(ncid, dimname, NF90_UNLIMITED, dimid)
call nc_check(ret, routine, 'define unlimited dimension '//trim(dimname), context, filename, ncid)

end subroutine nc_define_unlimited_dimension

!--------------------------------------------------------------------

function nc_get_dimension_size(ncid, dimname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: dimname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer :: nc_get_dimension_size

character(len=*), parameter :: routine = 'nc_get_dimension_size'
integer :: ret, dimid

ret = nf90_inq_dimid(ncid, dimname, dimid)
call nc_check(ret, routine, 'inq dimid '//trim(dimname), context, filename, ncid)

ret = nf90_inquire_dimension(ncid, dimid, len=nc_get_dimension_size)
call nc_check(ret, routine, 'inquire dimension '//trim(dimname), context, filename, ncid)

end function nc_get_dimension_size

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! defining variables section

!> unfortunately, the scalar versions of these routines cannot be
!> overloaded with the Nd versions.  the optional arguments make
!> the signatures (combinations of arguments) inseperable. 
!> it's less common to define scalars in netcdf files, so those ones 
!> get a separate entry point.

!--------------------------------------------------------------------

subroutine nc_define_var_char_1d(ncid, varname, dimname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_char_1d'
integer :: ret, dimid, varid

ret = nf90_inq_dimid(ncid, dimname, dimid)
call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimname), context, filename, ncid)

ret = nf90_def_var(ncid, varname, nf90_char, dimid, varid)
call nc_check(ret, routine, 'define character variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_char_1d

!--------------------------------------------------------------------

subroutine nc_define_var_char_Nd(ncid, varname, dimnames, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimnames(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_char_Nd'
integer :: ret, dimid1, dimid2, dimid3, varid

if (size(dimnames) >= 1) then
   ret = nf90_inq_dimid(ncid, dimnames(1), dimid1)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(1)), context, filename, ncid)
endif

if (size(dimnames) >= 2) then
   ret = nf90_inq_dimid(ncid, dimnames(2), dimid2)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(2)), context, filename, ncid)
endif

if (size(dimnames) >= 3) then
   ret = nf90_inq_dimid(ncid, dimnames(3), dimid3)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(3)), context, filename, ncid)
endif

if (size(dimnames) >= 4) then
   call error_handler(E_ERR, routine, 'only 1d, 2d and 3d character variables supported', &
                      source, text2='variable '//trim(varname))
endif

if (size(dimnames) == 1) then
   ret = nf90_def_var(ncid, varname, nf90_char, dimid1, varid=varid)
else if (size(dimnames) == 2) then
   ret = nf90_def_var(ncid, varname, nf90_char, dimids=(/ dimid1, dimid2 /), varid=varid)
else if (size(dimnames) == 3) then
   ret = nf90_def_var(ncid, varname, nf90_char, dimids=(/ dimid1, dimid2, dimid3 /), varid=varid)
endif

call nc_check(ret, routine, 'define character variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_char_Nd

!--------------------------------------------------------------------

subroutine nc_define_integer_scalar(ncid, varname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_integer_scalar'
integer :: ret, varid

ret = nf90_def_var(ncid, varname, nf90_int, varid=varid)
call nc_check(ret, routine, 'define scalar integer variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_integer_scalar

!--------------------------------------------------------------------

subroutine nc_define_var_int_1d(ncid, varname, dimname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_int_1d'
integer :: ret, dimid, varid

ret = nf90_inq_dimid(ncid, dimname, dimid)
call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimname), context, filename, ncid)

ret = nf90_def_var(ncid, varname, nf90_int, dimid, varid)
call nc_check(ret, routine, 'define integer variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_int_1d

!--------------------------------------------------------------------
! fortran supports up to 7 dimensional arrays.  but in this set of 
! routines we only go up to 4D arrays.

subroutine nc_define_var_int_Nd(ncid, varname, dimnames, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimnames(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_int_Nd'
integer :: i, ret, ndims, varid, dimids(NF90_MAX_VAR_DIMS)

ndims = size(dimnames)
if (ndims > 4) then
   call error_handler(E_ERR, routine, 'only 1d, 2d, 3d and 4d integer variables supported', &
                      source, text2='variable '//trim(varname))
endif

do i=1, ndims
   write(msgstring1,*)'"'//trim(varname)//'" inquire dimension id for dim "'//trim(dimnames(i))//'"'
   ret = nf90_inq_dimid(ncid, dimnames(i), dimids(i))
   call nc_check(ret, routine, msgstring1, context, filename, ncid)
enddo

ret = nf90_def_var(ncid, varname, nf90_int, dimids(1:ndims), varid=varid)
call nc_check(ret, routine, 'define integer variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_int_Nd

!--------------------------------------------------------------------

subroutine nc_define_real_scalar(ncid, varname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_real_scalar'
integer :: ret, varid

ret = nf90_def_var(ncid, varname, nf90_real, varid=varid)
call nc_check(ret, routine, 'define scalar real variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_real_scalar

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
call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimname), context, filename, ncid)

ret = nf90_def_var(ncid, varname, nf90_real, dimid, varid)
call nc_check(ret, routine, 'define real variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_real_1d

!--------------------------------------------------------------------

subroutine nc_define_var_real_Nd(ncid, varname, dimnames, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimnames(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_real_Nd'
integer :: i, ret, ndims, varid, dimids(NF90_MAX_VAR_DIMS)

ndims = size(dimnames)
if (ndims > 4) then
   call error_handler(E_ERR, routine, 'only 1d, 2d, 3d and 4d real variables supported', &
                      source, text2='variable '//trim(varname))
endif

do i=1, ndims
   write(msgstring1,*)'"'//trim(varname)//'" inquire dimension id for dim "'//trim(dimnames(i))//'"'
   ret = nf90_inq_dimid(ncid, dimnames(i), dimids(i))
   call nc_check(ret, routine, msgstring1, context, filename, ncid)
enddo

ret = nf90_def_var(ncid, varname, nf90_real, dimids(1:ndims), varid=varid)
call nc_check(ret, routine, 'define real variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_real_Nd

!--------------------------------------------------------------------

subroutine nc_define_double_scalar(ncid, varname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_double_scalar'
integer :: ret, varid

ret = nf90_def_var(ncid, varname, nf90_double, varid=varid)
call nc_check(ret, routine, 'define scalar double variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_double_scalar

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
call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimname), context, filename, ncid)

ret = nf90_def_var(ncid, varname, nf90_double, dimid, varid)
call nc_check(ret, routine, 'define double variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_double_1d

!--------------------------------------------------------------------

subroutine nc_define_var_double_Nd(ncid, varname, dimnames, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimnames(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_double_Nd'
integer :: i, ret, ndims, varid, dimids(NF90_MAX_VAR_DIMS)

ndims = size(dimnames)
if (ndims > 4) then
   call error_handler(E_ERR, routine, 'only 1d, 2d, 3d and 4d double variables supported', &
                      source, text2='variable '//trim(varname))
endif

do i=1, ndims
   write(msgstring1,*)'"'//trim(varname)//'" inquire dimension id for dim "'//trim(dimnames(i))//'"'
   ret = nf90_inq_dimid(ncid, dimnames(i), dimids(i))
   call nc_check(ret, routine, msgstring1, context, filename, ncid)
enddo

ret = nf90_def_var(ncid, varname, nf90_double, dimids(1:ndims), varid=varid)
call nc_check(ret, routine, 'define double variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_double_Nd

!--------------------------------------------------------------------
! check if vars, dims, or global atts exist (without error if not)
! these are functions, unlike the rest of these routines.

function nc_global_attribute_exists(ncid, attname)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: attname
logical                      :: nc_global_attribute_exists

character(len=*), parameter :: routine = 'nc_global_attribute_exists'
integer :: ret

ret = nf90_inquire_attribute(ncid, NF90_GLOBAL, attname)
nc_global_attribute_exists = (ret == NF90_NOERR) 

end function nc_global_attribute_exists

!--------------------------------------------------------------------

function nc_variable_attribute_exists(ncid, varname, attname)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
logical                      :: nc_variable_attribute_exists

character(len=*), parameter :: routine = 'nc_variable_attribute_exists'
integer :: varid, ret

ret = nf90_inq_varid(ncid, varname, varid)
nc_variable_attribute_exists = (ret == NF90_NOERR) 
if (ret /= NF90_NOERR) return

ret = nf90_inquire_attribute(ncid, varid, attname)
nc_variable_attribute_exists = (ret == NF90_NOERR) 

end function nc_variable_attribute_exists

!--------------------------------------------------------------------

function nc_dimension_exists(ncid, dimname)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: dimname
logical                      :: nc_dimension_exists

character(len=*), parameter :: routine = 'nc_dimension_exists'
integer :: ret, dimid

ret = nf90_inq_dimid(ncid, dimname, dimid)
nc_dimension_exists = (ret == NF90_NOERR)

end function nc_dimension_exists

!--------------------------------------------------------------------

function nc_variable_exists(ncid, varname)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
logical                      :: nc_variable_exists

character(len=*), parameter :: routine = 'nc_variable_exists'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
nc_variable_exists = (ret == NF90_NOERR) 

end function nc_variable_exists

!--------------------------------------------------------------------

subroutine nc_put_char_1d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: varvals
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_char_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_char_1d

!--------------------------------------------------------------------

subroutine nc_put_int_0d(ncid, varname, varindex, varval, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: varindex
integer,          intent(in) :: varval
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_int_0d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varval, start=(/varindex/))
call nc_check(ret, routine, 'put value for '//trim(varname), context, filename, ncid)

end subroutine nc_put_int_0d

!--------------------------------------------------------------------

subroutine nc_put_int_1d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_int_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_int_1d

!--------------------------------------------------------------------

subroutine nc_put_real_0d(ncid, varname, varval, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(r4),         intent(in) :: varval
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_real_0d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varval)
call nc_check(ret, routine, 'put value for '//trim(varname), context, filename, ncid)

end subroutine nc_put_real_0d

!--------------------------------------------------------------------

subroutine nc_put_real_1d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(r4),         intent(in) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_real_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_real_1d

!--------------------------------------------------------------------

subroutine nc_put_char_2d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_char_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_char_2d

!--------------------------------------------------------------------

subroutine nc_put_int_2d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_int_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_int_2d

!--------------------------------------------------------------------

subroutine nc_put_real_2d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(r4),         intent(in) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_real_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_real_2d

!--------------------------------------------------------------------

subroutine nc_put_int_3d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_int_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_int_3d

!--------------------------------------------------------------------

subroutine nc_put_real_3d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(r4),         intent(in) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_real_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_real_3d

!--------------------------------------------------------------------

subroutine nc_put_int_4d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: varvals(:,:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_int_4d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_int_4d

!--------------------------------------------------------------------

subroutine nc_put_real_4d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(r4),         intent(in) :: varvals(:,:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_real_4d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_real_4d

!--------------------------------------------------------------------

subroutine nc_put_double_0d(ncid, varname, varval, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(digits12),   intent(in) :: varval
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_double_0d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varval)
call nc_check(ret, routine, 'put value for '//trim(varname), context, filename, ncid)

end subroutine nc_put_double_0d

!--------------------------------------------------------------------

subroutine nc_put_double_1d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(digits12),   intent(in) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_double_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_double_1d

!--------------------------------------------------------------------

subroutine nc_put_double_2d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(digits12),   intent(in) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_double_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_double_2d

!--------------------------------------------------------------------

subroutine nc_put_double_3d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(digits12),   intent(in) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_double_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_double_3d

!--------------------------------------------------------------------

subroutine nc_put_double_4d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(digits12),   intent(in) :: varvals(:,:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_double_4d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_double_4d

!--------------------------------------------------------------------
! get values from variables
!
! check for scale/offset attributes and for now, error out.  eventually
! we could support in the _get_ routines to scale/offset them - but
! the return type won't be the same as the input.  e.g. might compute a
! real from input of short plus offset, scale factors.  maybe we don't
! ever want to support these - just punt and make the caller drop down
! into direct calls to the netcdf lib.

!--------------------------------------------------------------------

subroutine nc_get_short_0d(ncid, varname, varval, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer(i2),      intent(out) :: varval
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_short_0d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varval)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_short_0d

!--------------------------------------------------------------------

subroutine nc_get_short_1d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer(i2),      intent(out) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_short_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_short_1d

!--------------------------------------------------------------------

subroutine nc_get_int_0d(ncid, varname, varval, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varval
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_int_0d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varval)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_int_0d

!--------------------------------------------------------------------

subroutine nc_get_int_1d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_int_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_int_1d

!--------------------------------------------------------------------

subroutine nc_get_real_0d(ncid, varname, varval, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(r4),         intent(out) :: varval
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_real_0d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varval)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_real_0d

!--------------------------------------------------------------------

subroutine nc_get_real_1d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(r4),         intent(out) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_real_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_real_1d

!--------------------------------------------------------------------

subroutine nc_get_real_2d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)
integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(r4),         intent(out) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_real_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_real_2d

!--------------------------------------------------------------------

subroutine nc_get_real_3d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(r4),         intent(out) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_real_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_real_3d

!--------------------------------------------------------------------

subroutine nc_get_real_4d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(r4),         intent(out) :: varvals(:,:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_real_4d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_real_4d

!--------------------------------------------------------------------

subroutine nc_get_double_0d(ncid, varname, varvals, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(digits12),   intent(out) :: varvals
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_double_0d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_double_0d

!--------------------------------------------------------------------

subroutine nc_get_double_1d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(digits12),   intent(out) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_double_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_double_1d

!--------------------------------------------------------------------

subroutine nc_get_double_2d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)
integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(digits12),   intent(out) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_double_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_double_2d

!--------------------------------------------------------------------

subroutine nc_get_double_3d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(digits12),   intent(out) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_double_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_double_3d

!--------------------------------------------------------------------

subroutine nc_get_double_4d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(digits12),   intent(out) :: varvals(:,:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_double_4d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_double_4d

!--------------------------------------------------------------------

subroutine nc_get_short_2d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer(i2),      intent(out) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_short_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_short_2d

!--------------------------------------------------------------------

subroutine nc_get_int_2d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_int_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_int_2d

!--------------------------------------------------------------------

subroutine nc_get_short_3d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer(i2),      intent(out) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_short_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_short_3d

!--------------------------------------------------------------------

subroutine nc_get_int_3d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_int_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_int_3d

!--------------------------------------------------------------------

subroutine nc_get_int_4d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varvals(:,:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_int_4d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

! don't support variables which are supposed to have the values multiplied and shifted.
if (has_scale_off(ncid, varid)) call no_scale_off(ncid, routine, varname, context, filename)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_int_4d

!--------------------------------------------------------------------
!> Query and return information about a netCDF variable given the variable name.
!> Optionally returns the type of variable, the number of dimensions, 
!> the dimension names and lengths, the number of attributes (but not the attribute values (yet))

subroutine nc_get_variable_info(ncid, varname, xtype, ndims, dimlens, dimnames, natts, &
                                context, filename)

integer,          intent(in)            :: ncid
character(len=*), intent(in)            :: varname
integer,          intent(out), optional :: xtype
integer,          intent(out), optional :: ndims
integer,          intent(out), optional :: dimlens(:)
character(len=*), intent(out), optional :: dimnames(:)
integer,          intent(out), optional :: natts
character(len=*), intent(in) , optional :: context
character(len=*), intent(in) , optional :: filename

character(len=*), parameter :: routine = 'nc_get_variable_info'
integer :: ret, varid, dimid, ii

integer :: myndims
integer :: mydimids(NF90_MAX_VAR_DIMS)
integer :: mydimlens(NF90_MAX_VAR_DIMS)
character(len=NF90_MAX_NAME) :: mydimnames(NF90_MAX_VAR_DIMS)

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inq_varid for '//trim(varname), context, filename)

ret = nf90_inquire_variable(ncid, varid, xtype=xtype, ndims=myndims, &
                            dimids=mydimids, natts=natts) 
call nc_check(ret, routine, 'inquire_variable for '//trim(varname), context, filename)

if (present(dimlens) .or. present(dimnames)) then  ! more work to do 

   !>@todo do we want to make sure dimlens, dimnames are long enough
   if (present(dimlens))  dimlens  = 0
   if (present(dimnames)) dimnames = 'null'

   do ii = 1,myndims

      dimid = mydimids(ii)
      ret = nf90_inquire_dimension(ncid, dimid, name=mydimnames(ii), len=mydimlens(ii))

      write(msgstring1,*)'inquire_dimension ',ii,' for "'//trim(varname)//'"'
      call nc_check(ret, routine, msgstring1, context, filename)

   enddo
endif

if (present(   ndims)) ndims             = myndims
if (present( dimlens)) dimlens(1:ndims)  = mydimlens(1:ndims)
if (present(dimnames)) dimnames(1:ndims) = mydimnames(1:ndims)

end subroutine nc_get_variable_info

!--------------------------------------------------------------------

subroutine nc_get_variable_size_1d(ncid, varname, varsize, context, filename)      

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varsize
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_variable_size_1d'
integer :: ret, ndims, varid, dimids(NF90_MAX_VAR_DIMS)

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_inquire_variable(ncid, varid, dimids=dimids, ndims=ndims)
call nc_check(ret, routine, 'inquire dimensions for variable '//trim(varname), context, filename, ncid)

if (ndims /= 1) &
   call nc_check(NF90_EDIMSIZE, routine, 'variable '//trim(varname)//' not 1 dimensional', &
                 context, filename, ncid)

ret = nf90_inquire_dimension(ncid, dimids(1), len=varsize)
call nc_check(ret, routine, 'inquire dimension length for dimension 1', context, filename, ncid)

end subroutine nc_get_variable_size_1d

!--------------------------------------------------------------------

subroutine nc_get_variable_size_Nd(ncid, varname, varsize, context, filename)      

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varsize(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_variable_size_Nd'
integer :: ret, i, ndims, varid, dimids(NF90_MAX_VAR_DIMS)


ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_inquire_variable(ncid, varid, dimids=dimids, ndims=ndims)
call nc_check(ret, routine, 'inquire dimensions for variable '//trim(varname), context, filename, ncid)

if (ndims > size(varsize)) &
   call nc_check(NF90_EDIMSIZE, routine, 'variable '//trim(varname)//' return varsize array too small', &
                 context, filename, ncid)

! in case the var is larger than ndims, set unused dims to -1
varsize(:) = -1
do i=1, ndims
   ret = nf90_inquire_dimension(ncid, dimids(i), len=varsize(i))
   call nc_check(ret, routine, 'inquire dimension length for variable '//trim(varname), &
                 context, filename, ncid)
enddo

end subroutine nc_get_variable_size_Nd

!------------------------------------------------------------------

subroutine nc_get_variable_num_dimensions(ncid, varname, numdims, context, filename) 

integer, intent(in) :: ncid
character(len=*), intent(in):: varname
integer, intent(out) :: numdims
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_variable_num_dimensions'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_inquire_variable(ncid, varid, ndims=numdims)
call nc_check(ret, routine, 'inquire dimensions for variable '//trim(varname), context, filename, ncid)

end subroutine nc_get_variable_num_dimensions 

!------------------------------------------------------------------

subroutine nc_get_variable_dimension_names(ncid, varname, dimnames, context, filename) 

integer, intent(in) :: ncid
character(len=*), intent(in):: varname
character(len=*), intent(out) :: dimnames(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_variable_dimension_names'
integer :: ret, i, ndims, varid, dimids(NF90_MAX_VAR_DIMS)


ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_inquire_variable(ncid, varid, dimids=dimids, ndims=ndims)
call nc_check(ret, routine, 'inquire dimensions for variable '//trim(varname), context, filename, ncid)

if (ndims > size(dimnames)) &
   call nc_check(NF90_EDIMSIZE, routine, 'variable '//trim(varname)//' return dimnames array too small', &
                 context, filename, ncid)

! in case the var is larger than ndims, set unused dims to ""
dimnames(:) = ""
do i=1, ndims
   ret = nf90_inquire_dimension(ncid, dimids(i), name=dimnames(i))
   call nc_check(ret, routine, 'inquire dimension name for variable '//trim(varname), &
                 context, filename, ncid)
enddo

end subroutine nc_get_variable_dimension_names 

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

! the opens are the only routines in which filename is not the last argument.
! all other start with ncid.  this one starts with filename and it's required.
! it is also a function that returns an integer.  (caller doesn't need netcdf)

function nc_open_file_readonly(filename, context)

character(len=*), intent(in)  :: filename
character(len=*), intent(in), optional :: context
integer                       :: nc_open_file_readonly

character(len=*), parameter :: routine = 'nc_open_file_readonly'
integer :: ret, ncid

ret = nf90_open(filename, NF90_NOWRITE, ncid)
call nc_check(ret, routine, 'open '//trim(filename)//' read only', context)

call add_fh_to_list(ncid, filename)
nc_open_file_readonly = ncid

end function nc_open_file_readonly

!--------------------------------------------------------------------

! the opens are the only routines in which filename is not the last argument.
! all other start with ncid.  this one starts with filename and it's required.

function nc_open_file_readwrite(filename, context)

character(len=*), intent(in)  :: filename
character(len=*), intent(in), optional :: context
integer                       :: nc_open_file_readwrite

character(len=*), parameter :: routine = 'nc_open_file_readwrite'
integer :: ret, ncid, oldmode

ret = nf90_open(filename, NF90_WRITE, ncid)
call nc_check(ret, routine, 'open '//trim(filename)//' read/write', context)

call add_fh_to_list(ncid, filename)
nc_open_file_readwrite = ncid

! faster if we don't fill the vars first with 'fill' value.
! this works if we are planning to write all vars.  (which we are.)

ret = nf90_set_fill(ncid, NF90_NOFILL, oldmode)
call nc_check(ret, routine, 'set nofill mode', context, ncid=ncid)

end function nc_open_file_readwrite

!--------------------------------------------------------------------

! the opens are the only routines in which filename is not the last argument.
! all other start with ncid.  this one starts with filename and it's required.

function nc_create_file(filename, context)

character(len=*), intent(in)  :: filename
character(len=*), intent(in), optional :: context
integer                       :: nc_create_file

character(len=*), parameter :: routine = 'nc_create_file'
integer :: ret, ncid, oldmode

ret = nf90_create(filename, NF90_CLOBBER, ncid)
call nc_check(ret, routine, 'create '//trim(filename)//' read/write', context)

call add_fh_to_list(ncid, filename)
nc_create_file = ncid

! faster if we don't fill the vars first with 'fill' value.
! this works if we are planning to write all vars.  (which we are.)

ret = nf90_set_fill(ncid, NF90_NOFILL, oldmode)
call nc_check(ret, routine, 'set nofill mode', context, ncid=ncid)

end function nc_create_file

!--------------------------------------------------------------------

subroutine nc_close_file(ncid, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_close_file'
integer :: ret

ret = nf90_close(ncid)
call nc_check(ret, routine, 'close file', context, filename, ncid)

call del_fh_from_list(ncid)

end subroutine nc_close_file

!--------------------------------------------------------------------

subroutine nc_begin_define_mode(ncid, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_begin_define_mode'
integer :: ret

ret = nf90_redef(ncid)
call nc_check(ret, routine, 'begin file define mode', context, filename, ncid)

end subroutine nc_begin_define_mode

!--------------------------------------------------------------------

subroutine nc_end_define_mode(ncid, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_end_define_mode'
integer :: ret

ret = nf90_enddef(ncid)
call nc_check(ret, routine, 'end file define mode', context, filename, ncid)

end subroutine nc_end_define_mode

!--------------------------------------------------------------------

subroutine nc_synchronize_file(ncid, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_synchronize_file'
integer :: ret

ret = nf90_sync(ncid)
call nc_check(ret, routine, 'synchronize file contents', context, filename, ncid)

end subroutine nc_synchronize_file

!--------------------------------------------------------------------
!> check for the existence of either (or both) scale/offset attributes

function has_scale_off(ncid, varid)
integer, intent(in) :: ncid
integer, intent(in) :: varid
logical :: has_scale_off

integer :: ret

ret = nf90_inquire_attribute(ncid, varid, "scale_factor")
has_scale_off = (ret == NF90_NOERR) 
if (ret == NF90_NOERR) return

ret = nf90_inquire_attribute(ncid, varid, "add_offset")
has_scale_off = (ret == NF90_NOERR) 

end function has_scale_off

!------------------------------------------------------------------
!> don't support reading in that variable with this code
!> if either or both are found.  caller needs to go straight
!> to the netcdf interfaces to get the right types

subroutine no_scale_off(ncid, routine, varname, context, filename)
integer, intent(in) :: ncid
character(len=*), intent(in) :: routine
character(len=*), intent(in) :: varname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

call nc_check(NF90_ERANGE, routine, &
   'variable '//trim(varname)//' attrs "scale_factor","add_offset" unsupported', &
   context, filename, ncid)

end subroutine no_scale_off

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! internal bookkeeping of mapping between ncids and filenames

!------------------------------------------------------------------
!> add file handle and filename to an available slot

subroutine add_fh_to_list(ncid, filename)
integer, intent(in) :: ncid
character(len=*), intent(in) :: filename

integer :: i

do i=1, MAX_NCFILES
   if (ncinfo(i)%file_handle /= FH_EMPTY) cycle

   ncinfo(i)%file_handle = ncid
   ncinfo(i)%file_name = filename
   return 
enddo

end subroutine add_fh_to_list

!------------------------------------------------------------------
!> remove an entry when file is closed

subroutine del_fh_from_list(ncid)
integer, intent(in) :: ncid

integer :: i

do i=1, MAX_NCFILES
   if (ncinfo(i)%file_handle /= ncid) cycle

   ncinfo(i)%file_handle = FH_EMPTY
   ncinfo(i)%file_name = ''
   return 
enddo

end subroutine del_fh_from_list

!------------------------------------------------------------------
!> look up and return the filename for this handle

subroutine find_name_from_fh(ncid, filename)
integer, intent(in) :: ncid
character(len=*), intent(out) :: filename

integer :: i

do i=1, MAX_NCFILES
   if (ncinfo(i)%file_handle /= ncid) cycle

   filename = ncinfo(i)%file_name 
   return 
enddo

filename = ''

end subroutine find_name_from_fh

!------------------------------------------------------------------

end module netcdf_utilities_mod

