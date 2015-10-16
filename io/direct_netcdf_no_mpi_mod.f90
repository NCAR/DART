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

use state_structure_mod,  only : get_domain_size, get_num_variables 
                                 

use copies_on_off_mod,    only : query_read_copy, query_write_copy

use io_filenames_mod,     only : get_input_file, get_output_file

use direct_netcdf_common_mod,  only : write_variables, write_variables_clamp, &
                                      read_variables, create_state_output

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
      call read_variables(ncfile, vector, 1, get_num_variables(domain), domain)
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

integer :: create_mode

! need to read into a tempory array to fill with one copies
allocate(vector(get_domain_size(domain)))

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

         ! What file options do you want
         create_mode = ior(NF90_CLOBBER, NF90_64BIT_OFFSET)
         ret = nf90_create(netcdf_filename_out, create_mode, ncfile_out)
         call nc_check(ret, 'transpose_write: creating', trim(netcdf_filename_out))

         dart_time = state_ens_handle%time(1)
         call create_state_output(ncfile_out, domain, dart_time, write_single_precision)
      endif

   endif

   block_size = get_domain_size(domain)

   ending_point = starting_point + block_size -1

   if (query_write_copy(copy)) then

      vector = state_ens_handle%copies(copy, starting_point:ending_point)

      if (copy <= state_ens_handle%num_copies - state_ens_handle%num_extras) then
         ! actual copy, may need clamping
         call write_variables_clamp(ncfile_out, vector, 1, get_num_variables(domain), domain)
      else ! extra copy, don't clamp
         call write_variables(ncfile_out, vector, 1, get_num_variables(domain), domain)
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

end module direct_netcdf_mod
