!-------------------------------------------------------------------------------
! module for common routines shared by both
! direct_netcdf_mpi_mod and direct_netcdf_no_mpi_mod
!-------------------------------------------------------------------------------

module direct_netcdf_common_mod

use types_mod,           only : r8, missing_r8, digits12
use utilities_mod,       only : E_MSG, error_handler, nc_check
use time_manager_mod,    only : time_type
use state_structure_mod, only : get_variable_name, get_clamping_maxval,   &
                                get_clamping_minval, do_clamping,         &
                                get_num_dims, get_dim_lengths,            &
                                get_variable_size, get_num_unique_dims,   &
                                get_unique_dim_name, get_dim_name,        &
                                get_unique_dim_length, get_num_variables, &
                                set_var_id
use model_mod,           only : write_model_time

use netcdf

private

public :: clamp_variable,        &
          read_variables,        &
          write_variables,       &
          write_variables_clamp, &
          create_state_output

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_trunk/io/direct_netcdf_common_mod.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 8813 $"
character(len=128), parameter :: revdate  = "$Date: 2015-10-13 17:26:07 -0600 (Tue, 13 Oct 2015) $"

character(len=256) :: string1

contains

!-------------------------------------------------------------------------------
!> Check a variable for out of bounds and clamp or fail if needed
!-------------------------------------------------------------------------------
subroutine clamp_variable(dom_id, var_index, variable)

integer,     intent(in) :: dom_id ! domain id
integer,     intent(in) :: var_index ! variable index
real(r8), intent(inout) :: variable(:) ! variable

real(r8) :: minclamp, maxclamp
character(len=NF90_MAX_NAME) :: varname ! for debugging only

! is lower bound set
minclamp = get_clamping_minval(dom_id, var_index)
if ( minclamp /= missing_r8 ) then
   !>@todo Calculating both minval, and max might be expensive for
   !> large variables.  Should we handle this in a different way?
   !> do we want to have a message if the variable is being clamped?
   if ( minval(variable) < minclamp ) then
      
      !>@todo info about clamped variables is only printed on task 0
      !> E_MSG, but clamping may be done by other tasks
      varname = get_variable_name(dom_id, var_index) 
      write(string1, *) 'min val = ', minval(variable), &
                        'min bounds = ', minclamp
      call error_handler(E_MSG, 'clamp_variable', &
                  'Clamping '//trim(varname)//', values out of bounds.', &
                   source,revision,revdate, text2=string1)

      variable = max(minclamp, variable)
   endif
endif ! min range set

! is upper bound set
maxclamp = get_clamping_maxval(dom_id, var_index)
if ( maxclamp /= missing_r8 ) then
   if ( maxval(variable) > maxclamp ) then
      varname = get_variable_name(dom_id, var_index) 
      write(string1, *) 'max val = ', maxval(variable), &
                        'max bounds = ', maxclamp
      call error_handler(E_MSG, 'clamp_variable', &
                  'Clamping '//trim(varname)//', values out of bounds.', &
                   source,revision,revdate, text2=string1)

      variable = min(maxclamp, variable)
   endif
endif ! max range set

end subroutine clamp_variable

!-------------------------------------------------------------------------------
!> Read in variables from start_var to end_var
!> FIXME: At the moment, this code is assuming that the variables in the state start
!> at (1,1,1) and that the whole variable is read. This is not the case for 
!> Tiegcm and CLM.  
!-------------------------------------------------------------------------------
subroutine read_variables(ncfile_in, var_block, start_var, end_var, domain)

integer,  intent(in)    :: ncfile_in
real(r8), intent(inout) :: var_block(:)
integer,  intent(in)    :: start_var
integer,  intent(in)    :: end_var
integer,  intent(in)    :: domain

integer :: i
integer :: start_in_var_block, end_in_var_block
integer :: var_size
integer, allocatable :: dims(:)
integer :: ret, var_id

start_in_var_block = 1

do i = start_var, end_var

   var_size = get_variable_size(domain, i)
   end_in_var_block = start_in_var_block + var_size - 1

   ! number of dimensions and length of each
   allocate(dims(get_num_dims(domain, i)))

   dims = get_dim_lengths(domain, i)

   ret = nf90_inq_varid(ncfile_in, get_variable_name(domain, i), var_id)
   call nc_check(ret, 'read_variables: nf90_inq_varid',trim(get_variable_name(domain,i)) )

   ret = nf90_get_var(ncfile_in, var_id, var_block(start_in_var_block:end_in_var_block), count=dims)
   call nc_check(ret, 'read_variables: nf90_get_var',trim(get_variable_name(domain,i)) )

   start_in_var_block = start_in_var_block + var_size

   deallocate(dims)

enddo

end subroutine read_variables

!-------------------------------------------------------------------------------
!> Write variables from start_var to end_var no clamping
!-------------------------------------------------------------------------------
subroutine write_variables(ncfile_out, var_block, start_var, end_var, domain)

integer,  intent(in)    :: ncfile_out
real(r8), intent(inout) :: var_block(:)
integer,  intent(in)    :: start_var
integer,  intent(in)    :: end_var
integer,  intent(in)    :: domain 

integer :: start_in_var_block, end_in_var_block
integer, allocatable :: dims(:)
integer :: i, ret, var_id, var_size

start_in_var_block = 1
do i = start_var, end_var

   var_size = get_variable_size(domain, i)
   end_in_var_block = start_in_var_block + var_size - 1

   ! number of dimensions and length of each
   allocate(dims(get_num_dims(domain, i)))

   dims = get_dim_lengths(domain, i)

   ret = nf90_inq_varid(ncfile_out, get_variable_name(domain, i), var_id)
   call nc_check(ret, 'write_variables', 'getting variable id')

   ret = nf90_put_var(ncfile_out, var_id, var_block(start_in_var_block:end_in_var_block), count=dims)
   call nc_check(ret, 'write_variables', 'writing')
   start_in_var_block = start_in_var_block + var_size

   deallocate(dims)

enddo

end subroutine write_variables

!-------------------------------------------------------------------------------
!> Write variables from start_var to end_var for actual ensemble members
!-------------------------------------------------------------------------------
subroutine write_variables_clamp(ncfile_out, var_block, start_var, end_var, domain)

integer,  intent(in) :: ncfile_out
real(r8), intent(inout) :: var_block(:)
integer,  intent(in) :: start_var
integer,  intent(in) :: end_var
integer,  intent(in) :: domain 

integer :: start_in_var_block, end_in_var_block
integer, allocatable :: dims(:)
integer :: i, ret, var_id, var_size

start_in_var_block = 1
do i = start_var, end_var

   var_size = get_variable_size(domain, i)
   end_in_var_block = start_in_var_block + var_size - 1

   ! check whether you have to do anything to the variable, clamp or fail
   if ( do_clamping(domain, i) ) then
      call clamp_variable(domain, i, var_block(start_in_var_block:end_in_var_block))
   endif

   ! number of dimensions and length of each
   allocate(dims(get_num_dims(domain, i)))

   dims = get_dim_lengths(domain, i)

   ret = nf90_inq_varid(ncfile_out, get_variable_name(domain, i), var_id)
   call nc_check(ret, 'write_variables_clamp', 'getting variable id')

   ret = nf90_put_var(ncfile_out, var_id, var_block(start_in_var_block:end_in_var_block), count=dims)
   call nc_check(ret, 'write_variables_clamp', 'writing')
   start_in_var_block = start_in_var_block + var_size

   deallocate(dims)

enddo

end subroutine write_variables_clamp

!-------------------------------------------------------------------------------
!> Create the output files
!> I have removed fresh_netcdf_file, since filter_write_restart_direct can add a 
!> blank domain.
!>
!> A 'blank' domain is one variable called state, with dimension = model size.
!> It is used when the model has not suppled any netdcf info but 
!>     direct_netcdf_write = .true.
!>@todo The file is not closed here.
!-------------------------------------------------------------------------------
subroutine create_state_output(ncfile_out, dom_id, dart_time, single_precision_output)

integer,         intent(in) :: ncfile_out
integer,         intent(in) :: dom_id !< domain, not sure whether you need this?
type(time_type), intent(in) :: dart_time
logical,         intent(in) :: single_precision_output

integer :: ret !> netcdf return code
integer :: create_mode
integer :: i, j ! loop variables
integer :: new_dimid
integer :: new_varid
integer :: ndims
integer :: xtype ! precision for netcdf file
integer :: dimids(NF90_MAX_VAR_DIMS)

! define dimensions, loop around unique dimensions
do i = 1, get_num_unique_dims(dom_id)
   ret = nf90_def_dim(ncfile_out, get_unique_dim_name(dom_id, i), get_unique_dim_length(dom_id, i), new_dimid)
   !> @todo if we already have a unique names we can take this test out
   if(ret /= NF90_NOERR .and. ret /= NF90_ENAMEINUSE) then
      call nc_check(ret, 'create_state_output', 'defining dimensions')
   endif
enddo

! define variables
do i = 1, get_num_variables(dom_id) ! loop around state variables
   ! double or single precision?
   ndims = get_num_dims(dom_id, i)

   if (single_precision_output) then
      xtype = nf90_real
   else ! write output that is the precision of filter
      if (r8 == digits12) then ! datasize = MPI_REAL8  ! What should we be writing?
         xtype = nf90_double
      else
         xtype = nf90_real
      endif
   endif

   ! query the dimension ids
   do j = 1, get_num_dims(dom_id, i)
      ret = nf90_inq_dimid(ncfile_out, get_dim_name(dom_id, i, j), dimids(j))
      call nc_check(ret, 'create_state_output', 'querying dimensions')
   enddo

   ret = nf90_def_var(ncfile_out, trim(get_variable_name(dom_id, i)), &
                      xtype=xtype, dimids=dimids(1:get_num_dims(dom_id, i)), &
                      varid=new_varid)

   call nc_check(ret, 'create_state_output', 'defining variable')
      !variable_ids(i, dom_id) = new_varid
   call set_var_id(dom_id, i, new_varid)
enddo

ret = nf90_enddef(ncfile_out)
call nc_check(ret, 'create_state_output', 'end define mode')

call write_model_time(ncfile_out, dart_time)

end subroutine create_state_output

!-------------------------------------------------------------------------------

end module direct_netcdf_common_mod

