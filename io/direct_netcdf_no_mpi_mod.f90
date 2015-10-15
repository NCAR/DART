!> At the moment this is to have a null version of read_transpose
!> and transpose_write for programs like closest_member_tool that
!> read in an state vector (old school) but don't use mpi.  
!> These programs need to be able to use state_vector_io_mod but 
!> without mpi.
!> Not sure if this is the best way to organize the code.
!> Should this module reaad from netcdf and put into a copies array on a single task?
module direct_netcdf_mod

use types_mod,            only : r8, digits12, missing_r8
use ensemble_manager_mod, only : ensemble_type

use state_structure_mod,  only : get_domain_size, get_variable_name, &
                                 get_dim_lengths, get_num_variables, &
                                 get_variable_size, get_num_dims, &
                                 get_dim_name, get_unique_dim_length, &
                                 get_unique_dim_name, get_num_unique_dims, &
                                 set_var_id, get_clamping_minval, &
                                 get_clamping_maxval,     &
                                 do_clamping

use model_mod,            only : write_model_time

use copies_on_off_mod,    only : query_read_copy, query_write_copy

use io_filenames_mod,     only : get_input_file, get_output_file

use utilities_mod,        only : error_handler, nc_check, &
                                 E_MSG, E_ERR, file_exist

use time_manager_mod,     only : time_type, print_time

use netcdf

implicit none

private

public :: read_transpose, transpose_write

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=256) :: string1
character(len=256) :: string2


integer :: ret !< netcdf return code
integer :: ncfile !< netcdf input file identifier
integer :: ncfile_out !< netcdf output file handle
character(len=256) :: netcdf_filename
character(len=256) :: netcdf_filename_out
logical :: single_precision_output = .false. ! Allows you to write r4 netcdf files even if filter is double precision


contains

!-------------------------------------------------
!> Single processor version of read_transpose.  Reads ens_size whole vectors from
!> netcdf files and fills up a row of %copies for each file.
subroutine read_transpose(state_ens_handle, domain, dart_index, read_limit_mem, read_limit_procs)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: domain
integer,             intent(inout) :: dart_index !< This is for mulitple domains
! The following two arguments are not used in the non-mpi version
integer,             intent(in)    :: read_limit_mem !< How many state elements you can read at once
integer,             intent(in)    :: read_limit_procs !< How many processors are involved in the transpose


real(r8), allocatable :: vector(:)

integer :: block_size
integer :: starting_point, ending_point
integer :: copy
integer :: start_var

starting_point = dart_index ! position in state_ens_handle%vars
block_size = 0

! need to read into a tempory array, then fill up copies
allocate(vector(get_domain_size(domain)))

COPIES: do copy = 1, state_ens_handle%my_num_copies

   start_var = 1 ! read first variable first

   ! open netcdf file
   if (query_read_copy(copy)) then
      netcdf_filename = get_input_file(copy, domain)
      ret = nf90_open(netcdf_filename, NF90_NOWRITE, ncfile)
      call nc_check(ret, 'read_restart_netcdf opening', netcdf_filename)
   endif

   block_size = get_domain_size(domain)

   ending_point = starting_point + block_size -1

   if (query_read_copy(copy)) then
      call read_variables(vector, 1, get_num_variables(domain), domain)
      ! close netcdf file
      ret = nf90_close(ncfile)
      call nc_check(ret, 'read_restart_netcdf closing', netcdf_filename)
      state_ens_handle%copies(copy, starting_point:ending_point) = vector

   endif

enddo COPIES

! update starting point
starting_point = starting_point + block_size

dart_index = starting_point

deallocate(vector)

end subroutine read_transpose
!-------------------------------------------------

!> Single processor version of transpose write.  Takes copies array one row
!> at a time and writes copy to a netcdf file.
subroutine transpose_write(state_ens_handle, num_extras, domain, dart_index, isprior, write_limit_mem, write_limit_procs, write_single_precision)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: num_extras ! non restart copies
integer,             intent(in)    :: domain
integer,             intent(inout) :: dart_index
logical,             intent(in)    :: isprior
integer,             intent(in)    :: write_limit_mem !< How many state elements you can read at once
integer,             intent(in)    :: write_limit_procs !< How many processors are involved in the transpose
logical,             intent(in)    :: write_single_precision


real(r8), allocatable :: vector(:)

integer :: block_size
integer :: starting_point, ending_point
integer :: copy
integer :: start_var

character(len=256)      :: msgstring
type(time_type) :: dart_time

! need to read into a tempory array to fill with one copies
allocate(vector(get_domain_size(domain)))

single_precision_output = write_single_precision

starting_point = dart_index ! position in state_ens_handle%vars
block_size = 0

! need to read into a tempory array, then fill up copies

COPIES: do copy = 1, state_ens_handle%my_num_copies

   start_var = 1 ! read first variable first

   ! open netcdf file
   if (query_write_copy(copy)) then
      netcdf_filename_out = get_output_file(copy, domain, isprior)

      if(file_exist(netcdf_filename_out)) then
         ret = nf90_open(netcdf_filename_out, NF90_WRITE, ncfile_out)
         call nc_check(ret, 'transpose_write: opening', trim(netcdf_filename_out))
         else ! create output file if it does not exist
            write(msgstring, *) 'Creating output file ', trim(netcdf_filename_out)
         call error_handler(E_MSG,'state_vector_io_mod:', msgstring)
         dart_time = state_ens_handle%time(1)
         call create_state_output(netcdf_filename_out, domain, dart_time)
      endif

   endif

   block_size = get_domain_size(domain)

   ending_point = starting_point + block_size -1

   if (query_write_copy(copy)) then

      vector = state_ens_handle%copies(copy, starting_point:ending_point)

      if (copy <= state_ens_handle%num_copies - state_ens_handle%num_extras) then
         ! actual copy, may need clamping
         call write_variables_clamp(vector, 1, get_num_variables(domain), domain)
      else ! extra copy, don't clamp
         call write_variables(vector, 1, get_num_variables(domain), domain)
      endif
      ! close netcdf file
      ret = nf90_close(ncfile_out)
      call nc_check(ret, 'transpose_write closing', netcdf_filename_out)
   endif

enddo COPIES

! update starting point
starting_point = starting_point + block_size

dart_index = starting_point

deallocate(vector)

end subroutine transpose_write
!-------------------------------------------------

!-------------------------------------------------
! Read_variables and write_variables is duplicate code
! from direct_netcdf_mpi_mod.f90
!-------------------------------------------------
!> Read in variables from start_var to end_var
!> FIXME: At the moment, this code is assuming that the variables in the state start
!> at (1,1,1) and that the whole variable is read. This is not the case for 
!> Tiegcm and CLM.  
subroutine read_variables(var_block, start_var, end_var, domain)

real(r8),           intent(inout) :: var_block(:)
integer,            intent(in)    :: start_var
integer,            intent(in)    :: end_var
integer,            intent(in)    :: domain

integer :: i
integer :: start_in_var_block
integer :: var_size
integer, allocatable :: dims(:)
integer :: var_id

start_in_var_block = 1

do i = start_var, end_var

   var_size = get_variable_size(domain, i)

   ! number of dimensions and length of each
   allocate(dims(get_num_dims(domain, i)))

   dims = get_dim_lengths(domain, i)

   ret = nf90_inq_varid(ncfile, get_variable_name(domain, i), var_id)
   call nc_check(ret, 'read_variables','inquire variable id')

   ret = nf90_get_var(ncfile, var_id, var_block(start_in_var_block:start_in_var_block+var_size-1), count=dims)
   call nc_check(ret, 'read_variables','reading')

   start_in_var_block = start_in_var_block + var_size

   deallocate(dims)

enddo

end subroutine read_variables
!-------------------------------------------------
!> Write variables from start_var to end_var
!> no clamping
subroutine write_variables(var_block, start_var, end_var, domain)

real(r8), intent(inout) :: var_block(:)
integer,  intent(in) :: start_var
integer,  intent(in) :: end_var
integer,  intent(in) :: domain 

integer :: i
integer :: count_displacement
integer :: start_in_var_block
integer :: var_size
integer, allocatable :: dims(:)
integer :: var_id 

start_in_var_block = 1
do i = start_var, end_var

   var_size = get_variable_size(domain, i)

   ! number of dimensions and length of each
   allocate(dims(get_num_dims(domain, i)))

   dims = get_dim_lengths(domain, i)

   ret = nf90_inq_varid(ncfile_out, get_variable_name(domain, i), var_id)
   call nc_check(ret, 'write_variables', 'getting variable id')

   ret = nf90_put_var(ncfile_out, var_id, var_block(start_in_var_block:start_in_var_block+var_size-1), count=dims)
   call nc_check(ret, 'write_variables', 'writing')
   start_in_var_block = start_in_var_block + var_size

   deallocate(dims)

enddo

end subroutine write_variables

!-------------------------------------------------
!> Write variables from start_var to end_var
!> For actual ensemble members
subroutine write_variables_clamp(var_block, start_var, end_var, domain)

real(r8), intent(inout) :: var_block(:)
integer,  intent(in) :: start_var
integer,  intent(in) :: end_var
integer,  intent(in) :: domain 

integer :: i
integer :: count_displacement
integer :: start_in_var_block
integer :: var_size
integer, allocatable :: dims(:)
integer :: var_id 

start_in_var_block = 1
do i = start_var, end_var

   var_size = get_variable_size(domain, i)

   ! check whether you have to do anything to the variable, clamp or fail
   if (do_clamping(domain, i)) then
      call clamp_variable(domain, i, var_block(start_in_var_block:start_in_var_block+var_size-1))
   endif

   ! number of dimensions and length of each
   allocate(dims(get_num_dims(domain, i)))

   dims = get_dim_lengths(domain, i)

   ret = nf90_inq_varid(ncfile_out, get_variable_name(domain, i), var_id)
   call nc_check(ret, 'write_variables_clamp', 'getting variable id')

   ret = nf90_put_var(ncfile_out, var_id, var_block(start_in_var_block:start_in_var_block+var_size-1), count=dims)
   call nc_check(ret, 'write_variables_clamp', 'writing')
   start_in_var_block = start_in_var_block + var_size

   deallocate(dims)

enddo

end subroutine write_variables_clamp

!-------------------------------------------------
!> Create the output files
!> ncfile_out is global - is this ok?
!> It is used when the model has not suppled any netdcf info but direct_netcdf_write = .true.
!>@todo The file is not closed here.
subroutine create_state_output(filename, dom, dart_time)

character(len=256), intent(in) :: filename
integer,            intent(in) :: dom !< domain, not sure whether you need this?
type(time_type),    intent(in) :: dart_time

integer :: ret !> netcdf return code
integer :: create_mode
integer :: i, j ! loop variables
integer :: new_dimid
integer :: new_varid
integer :: ndims
integer :: xtype ! precision for netcdf file
logical :: time_dimension_exists
integer :: dimids(NF90_MAX_VAR_DIMS)

time_dimension_exists = .false.

! What file options do you want
create_mode = ior(NF90_CLOBBER, NF90_64BIT_OFFSET)
ret = nf90_create(filename, create_mode, ncfile_out)
call nc_check(ret, 'create_state_output', 'creating')

! define dimensions, loop around unique dimensions
do i = 1, get_num_unique_dims(dom)
   ret = nf90_def_dim(ncfile_out, get_unique_dim_name(dom, i), get_unique_dim_length(dom, i), new_dimid)
   !> @todo if we already have a unique names we can take this test out
   if(ret /= NF90_NOERR .and. ret /= NF90_ENAMEINUSE) then
      call nc_check(ret, 'create_state_output', 'defining dimensions')
   endif
enddo

! define variables
do i = 1, get_num_variables(dom) ! loop around state variables
   ! double or single precision?
   ndims = get_num_dims(dom, i)

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
   do j = 1, get_num_dims(dom, i)
      ret = nf90_inq_dimid(ncfile_out, get_dim_name(dom, i, j), dimids(j))
      call nc_check(ret, 'create_state_output', 'querying dimensions')
   enddo

   ret = nf90_def_var(ncfile_out, trim(get_variable_name(dom, i)), xtype=xtype, dimids=dimids(1:get_num_dims(dom, i)), varid=new_varid)
   call nc_check(ret, 'create_state_output', 'defining variable')
      !variable_ids(i, dom) = new_varid
   call set_var_id(dom, i, new_varid)
enddo

ret = nf90_enddef(ncfile_out)
call nc_check(ret, 'create_state_output', 'end define mode')

call write_model_time(ncfile_out, dart_time)

end subroutine create_state_output
!-------------------------------------------------

!-------------------------------------------------------
!> Check a variable for out of bounds and clamp or fail if
!> needed
subroutine clamp_variable(dom, var_index, variable)

integer,     intent(in) :: dom ! domain index
integer,     intent(in) :: var_index ! variable index
real(r8), intent(inout) :: variable(:) ! variable

real(r8) :: minclamp, maxclamp
character(len=NF90_MAX_NAME) :: varname ! for debugging only

! is lower bound set
minclamp = get_clamping_minval(dom, var_index)
if ( minclamp /= missing_r8 ) then
   if ( minval(variable) < minclamp ) then
      varname = get_variable_name(dom, var_index) 
      write(string1, *) 'min data val = ', minval(variable), &
                        'min data bounds = ', minclamp
      call error_handler(E_MSG, 'clamp_variable', &
                  'Clamping '//trim(varname)//', values out of bounds.', &
                   source,revision,revdate, text2=string1)

      variable = max(minclamp, variable)
   endif
endif ! min range set

! is upper bound set
maxclamp = get_clamping_maxval(dom, var_index)
if ( maxclamp /= missing_r8 ) then
   if ( maxval(variable) > maxclamp ) then
      varname = get_variable_name(dom, var_index) 
      write(string1, *) 'max data val = ', maxval(variable), &
                        'max data bounds = ', maxclamp
      call error_handler(E_MSG, 'clamp_variable', &
                  'Clamping '//trim(varname)//', values out of bounds.', &
                   source,revision,revdate, text2=string1)

      variable = min(maxclamp, variable)
   endif
endif ! max range set

end subroutine clamp_variable

!-------------------------------------------------

end module direct_netcdf_mod