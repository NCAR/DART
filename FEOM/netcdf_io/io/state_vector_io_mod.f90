! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> IO for the state vector. \n
!> Idea is to be generic. \n
!> Don't want to go to the filesystem twice wrf => wrf_to_dart => dart => dart_to_wrf \n
!> Get a list of variables in the state \n
!> Get info about their dimenions \n

module state_vector_io_mod

use types_mod,            only : r8, MISSING_R8, digits12

use mpi_utilities_mod,    only : task_count, send_to, receive_from, my_task_id

use ensemble_manager_mod, only : ensemble_type, map_pe_to_task

use utilities_mod,        only : error_handler, E_ERR, nc_check, check_namelist_read, &
                                 find_namelist_in_file, nmlfileunit, do_nml_file, do_nml_term, file_exist, &
                                 E_MSG

use assim_model_mod,      only : get_model_size, aread_state_restart, awrite_state_restart, &
                                 open_restart_read, open_restart_write, close_restart, &
                                 clamp_or_fail_it, do_clamp_or_fail

use time_manager_mod,     only : time_type

use netcdf

use io_filenames_mod,     only : restart_files_in, restart_files_out

use copies_on_off_mod


implicit none

private

public :: state_vector_io_init, &
          initialize_arrays_for_read, &
          get_state_variable_info, &
          read_restart_netcdf, &
          write_restart_netcdf, &
          netcdf_filename_out, &
          setup_read_write, &
          turn_read_copy_on, turn_write_copy_on, &
          turn_read_copies_off, turn_write_copy_off

integer :: ret !< netcdf return code
integer :: ncfile !< netcdf input file identifier
integer :: ncfile_out !< netcdf output file handle
character(len=256) :: netcdf_filename !< needs to be different for each task
character(len=256) :: netcdf_filename_out !< needs to be different for each task
integer, parameter :: MAXDIMS = NF90_MAX_VAR_DIMS 

integer,   allocatable :: variable_ids(:, :)
integer,   allocatable :: variable_sizes(:, :)
integer,   allocatable :: dimensions_and_lengths(:, :, :) !< number of dimensions and length of each dimension

integer :: num_state_variables
integer :: num_domains

! record dimensions for output netcdf file,
! (state_variable, dimension, domain)
integer,            allocatable :: dimIds(:, :, :) !< dimension ids
integer,            allocatable :: copy_dimIds(:, :, :) !< dimension ids copy
integer,            allocatable :: length(:, :, :) !< dimension length
character(len=256), allocatable :: dim_names(:, :, :)

! list of variables names in the state
character(len=256), allocatable :: global_variable_names(:)

! namelist variables with default values
logical :: single_precision_output = .false. ! Allows you to write r4 netcdf files even if filter is double precision
logical :: create_restarts = .false. ! what if the restart files exist?
logical :: time_unlimited = .true. ! You need to keep track of the time.

namelist /  state_vector_io_nml / single_precision_output, create_restarts, time_unlimited

contains

!-------------------------------------------------
!> Initialize model 
!> so you can read the namelist
subroutine state_vector_io_init()

integer :: iunit, io

!call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "state_vector_io_nml", iunit)
read(iunit, nml = state_vector_io_nml, iostat = io)
call check_namelist_read(iunit, io, "state_vector_io_nml")

! Write the namelist values to the log file
if (do_nml_file()) write(nmlfileunit, nml=state_vector_io_nml)
if (do_nml_term()) write(     *     , nml=state_vector_io_nml)

end subroutine state_vector_io_init

!-------------------------------------------------
!> Initialize arrays.  Need to know the number of domains
subroutine initialize_arrays_for_read(n, n_domains)

integer, intent(in) :: n !< number of state variables
integer, intent(in) :: n_domains !< number of domains (and therfore netcdf files) to read

if(allocated(variable_ids))   call error_handler(E_ERR, 'initialize_arrays_for_read', 'already called this routine')
if(allocated(variable_sizes)) call error_handler(E_ERR, 'initialize_arrays_for_read', 'already called this routine')

num_state_variables = n
num_domains = n_domains

allocate(variable_ids(n, num_domains), variable_sizes(n, num_domains))
allocate(dimIds(n, MAXDIMS, num_domains), length(n, MAXDIMS, num_domains), dim_names(n, MAXDIMS, num_domains))
allocate(copy_dimIds(n, MAXDIMS, num_domains))
allocate(dimensions_and_lengths(n, MAXDIMS +1, num_domains))
allocate(global_variable_names(n))

dimensions_and_lengths = -1  ! initialize to a nonsense value
dimIds = -1

end subroutine initialize_arrays_for_read

!-------------------------------------------------
!> Need list of variables in the state
!> Each task grabs the variable dimensions from the netcdf file
subroutine get_state_variable_info(n, variable_names, domain, domain_size)

integer,            intent(in)  :: n !< number of state variables
character(len=256), intent(in)  :: variable_names(n)
integer,            intent(in)  :: domain !< which domain info you are grabbing
integer,            intent(out) :: domain_size

integer :: i

! load up module storage with variable names
global_variable_names = variable_names

! open netcdf file - all restart files have the same info?
netcdf_filename = restart_files_in(1,domain)
ret = nf90_open(netcdf_filename, NF90_NOWRITE, ncfile)
call nc_check(ret, 'get_state_variable_info', 'opening')

! get all variable ids
call get_variable_ids(variable_names, domain, variable_ids(:, domain))

! get all variable sizes, only readers store dimensions?
variable_sizes(:, domain) = total_size(n, variable_ids(:, domain), domain)
domain_size = sum(variable_sizes(:, domain))

! close netcdf file
ret = nf90_close(ncfile)
call nc_check(ret, 'get_state_variable_info', 'closing')

end subroutine get_state_variable_info

!-------------------------------------------------
!> Get netcdf variable ids
subroutine get_variable_ids(variable_names, domain, variable_ids)

character(len=*), intent(in)    :: variable_names(:) !< netcdf variable names
integer,          intent(in)    :: domain ! where is this used.
integer,          intent(inout) :: variable_ids(:) !< netcdf variable ids

integer :: n !< number of variables in the state
integer :: i !< loop variable

n = size(variable_ids)

do i = 1, n

   ret = nf90_inq_varid(ncfile, variable_names(i), variable_ids(i))
   call nc_check(ret, 'get_variable_ids', 'inq_var_id')

enddo

end subroutine get_variable_ids

!-------------------------------------------------
!> Get state variable size, i.e. how many elements are in the variable
!> readers need to store the dimension information
!> Every task is storing the dimension information at the moment.
!> FIXME: this should be a subroutine since it has side effects.
function total_size(n, varId, domain)

integer, intent(in)  :: n !< number of variables in the state vector
integer, intent(in)  :: varId(n) !< variable ids
integer, intent(in)  :: domain
integer              :: total_size(n) ! NOTE RETURN VALUE LAST

integer              :: ndims !< number of dimensions
integer              :: i, j !< loop variable
integer              :: xtype !< variable type: NF90_FLOAT, NF90_DOUBLE, etc. I don't think we need to store this

do i = 1, n

   ret = nf90_inquire_variable(ncfile, varId(i), ndims=ndims, dimids=dimIds(i, :, domain), xtype=xtype)
   call nc_check(ret, 'totalsize', 'inq_var')

   do j = 1, ndims
      ret = nf90_inquire_dimension(ncfile, dimIds(i, j, domain), name=dim_names(i, j, domain), len=length(i, j, domain))
      call nc_check(ret, 'totalsize', 'inq_dimlen')
   enddo

   total_size(i) = product(length(i, 1:ndims, domain))

   dimensions_and_lengths(i, 1, domain) = ndims
   dimensions_and_lengths(i, 2:ndims+1, domain) = length(i, 1:ndims, domain)

enddo

end function total_size


!-------------------------------------------------
!> Read in variables from model restart file into state_ens_handle%vars
subroutine read_restart_netcdf(state_ens_handle, restart_in_file_name, domain, dart_index)

type(ensemble_type), intent(inout) :: state_ens_handle
character(len=129),  intent(in)    :: restart_in_file_name
integer,             intent(in)    :: domain
integer,             intent(inout) :: dart_index !< This is for mulitple domains

integer :: i
integer :: start_var, end_var !< start/end variables in a read block
integer :: my_pe !< task or pe?
integer :: recv_pe, sending_pe
real(r8), allocatable :: var_block(:) !< for reading in variables
integer :: block_size !< number of state elements in a block
integer :: count !< number of elements to send
integer :: starting_point!< position in state_ens_handle%copies
integer :: ending_point
integer :: ens_size !< ensemble size
integer :: remainder
integer :: start_rank
integer :: my_copy !< which copy a pe is reading
integer :: c !< copies_read loop index
integer :: copies_read

! single file
integer :: iunit
type(time_type) :: ens_time

ens_size = state_ens_handle%num_copies ! have the extras, incase you need to read inflation restarts

my_pe = state_ens_handle%my_pe

copies_read = 0

starting_point = dart_index ! position in state_ens_handle%vars
block_size = 0

COPIES: do c = 1, state_ens_handle%my_num_copies

   start_var = 1 ! read first variable first
   my_copy = state_ens_handle%my_copies(c)

   ! open netcdf file
   ! You have already opened this once to read the variable info. Should you just leave it open
   ! on the readers?
   if (query_read_copy(my_copy)) then
      netcdf_filename = restart_files_in(my_copy, domain)
      !print*, 'opening netcdf_filename ', trim(netcdf_filename)
      ret = nf90_open(netcdf_filename, NF90_NOWRITE, ncfile)
      call nc_check(ret, 'read_restart_netcdf opening', netcdf_filename)
   endif

   block_size = sum(variable_sizes(:, domain))

   ending_point = starting_point + block_size -1

   if (query_read_copy(my_copy)) then
      call read_variables(state_ens_handle%vars(starting_point:ending_point, c), 1, num_state_variables, domain)
      ! close netcdf file
      ret = nf90_close(ncfile)
      call nc_check(ret, 'read_restart_netcdf closing', netcdf_filename)
   endif

enddo COPIES

! update starting point
starting_point = starting_point + block_size

dart_index = starting_point

end subroutine read_restart_netcdf

!-------------------------------------------------
!> Write directly to netcdf file from %vars
!> The suffix argument can be used when writing state members for the prior
subroutine write_restart_netcdf(state_ens_handle, restart_out_file_name, domain, dart_index, isprior)

type(ensemble_type), intent(inout) :: state_ens_handle
character(len=129),  intent(in)    :: restart_out_file_name
integer,             intent(in)    :: domain
integer,             intent(inout) :: dart_index
logical,             intent(in)    :: isprior

integer :: i
integer :: start_var, end_var !< start/end variables in a read block
integer :: my_pe !< task or pe?
integer :: num_vars !< number of variables in a block
integer :: starting_point!< position in state_ens_handle%copies
integer :: ending_point
integer :: ens_size !< ensemble size
integer :: dummy_loop, j
integer :: my_copy !< which copy a pe is reading, starting from 0 to match pe
integer :: c !< copies_read loop index

character(len=256)      :: msgstring

! single file
integer :: iunit
type(time_type) :: ens_time

ens_size = state_ens_handle%num_copies ! have the extras incase you want to read inflation restarts
my_pe = state_ens_handle%my_pe

starting_point = dart_index ! position in state_ens_handle%copies
num_vars = 0

COPIES : do c = 1, state_ens_handle%my_num_copies

   my_copy = state_ens_handle%my_copies(c)

   ! writers open netcdf output file. This is a copy of the input file
   if ( query_write_copy(my_copy)) then
         if(isprior) then
            netcdf_filename_out = trim(restart_files_out((my_copy), domain, 1))
         else
            netcdf_filename_out = trim(restart_files_out((my_copy), domain, 2))
         endif
      if (create_restarts) then ! How do you want to do create restarts
         call create_state_output(netcdf_filename_out, domain)
      else
         if(file_exist(netcdf_filename_out)) then
            ret = nf90_open(netcdf_filename_out, NF90_WRITE, ncfile_out)
            call nc_check(ret, 'write_restart_netcdf opening', trim(netcdf_filename_out))
         else ! create output file if it does not exist
            write(msgstring, *) 'Creating output file ', trim(netcdf_filename_out)
            call error_handler(E_MSG,'state_vector_io_mod:', msgstring)
            call create_state_output(netcdf_filename_out, domain)
         endif
      endif
   endif

   num_vars = sum(variable_sizes(1:num_state_variables, domain))

   ending_point = starting_point + num_vars -1

   if (query_write_copy(my_copy)) then
      if (my_copy <= state_ens_handle%num_copies -6) then ! actual copy, may need clamping
         call write_variables_clamp(state_ens_handle%vars(starting_point:ending_point, c), 1, num_state_variables, domain)
      else ! extra copy, don't clamp
         call write_variables(state_ens_handle%vars(starting_point:ending_point, c), 1, num_state_variables, domain)
      endif
      ret = nf90_close(ncfile_out)
      call nc_check(ret, 'write_restart_netcdf', 'closing')
   endif

enddo COPIES

! update starting point
starting_point = starting_point + num_vars

dart_index = starting_point

end subroutine write_restart_netcdf

!-------------------------------------------------
!> Read in variables from start_var to end_var
!> FIXME: At the moment, this code is assuming that the variables in the state start
!> at (1,1,1) and that the whole variable is read. This is not the case for 
!> Tiegcm and CLM.  
subroutine read_variables(var_block, start_var, end_var, domain)

real(r8), intent(inout) :: var_block(:)
integer,  intent(in)    :: start_var
integer,  intent(in)    :: end_var
integer,  intent(in)    :: domain

integer :: i
integer :: start_in_var_block
integer :: var_size
integer, allocatable :: dims(:)
integer :: var_id

start_in_var_block = 1

do i = start_var, end_var

   var_size = variable_sizes(i, domain)

   ! number of dimensions and length of each
   allocate(dims(dimensions_and_lengths(i,1, domain)))

   dims = dimensions_and_lengths(i, 2:dimensions_and_lengths(i,1, domain) + 1, domain)
   ret = nf90_inq_varid(ncfile, global_variable_names(i), var_id)
   call nc_check(ret, 'read_variables','inquire variable id')

   ret = nf90_get_var(ncfile, var_id, var_block(start_in_var_block:start_in_var_block+var_size-1), count=dims)
   call nc_check(ret, 'read_variables','reading')

   start_in_var_block = start_in_var_block + var_size

   deallocate(dims)

enddo

end subroutine read_variables


!-------------------------------------------------
!> Create the output files
!> Assume same variables in each domain
!> Overwriting reading arrays with writing info.  Is this what you want to do?
!> ncfile_out is global - is this ok?
subroutine create_state_output(filename, dom)

character(len=256), intent(in) :: filename
integer,            intent(in) :: dom !< domain, not sure whether you need this?

integer :: ret !> netcdf return code
integer :: create_mode
integer :: i, j ! loop variables
integer :: new_dimid
integer :: new_varid
integer :: ndims
integer :: xtype ! precision for netcdf file
logical :: time_dimension_exists

time_dimension_exists = .false.

! What file options do you want
create_mode = NF90_64BIT_OFFSET
ret = nf90_create(filename, create_mode, ncfile_out)
call nc_check(ret, 'create_state_output', 'creating')

! make a copy of the dimIds array
copy_dimIds = dimIds

! define dimensions
do i = 1, num_state_variables ! loop around state variables
   ! check if their dimensions exist
   do j = 1, dimensions_and_lengths(i, 1, dom) ! ndims
      if (time_unlimited .and. (dim_names(i, j, dom) == 'Time')) then ! case sensitive
         ret = nf90_def_dim(ncfile_out, dim_names(i, j, dom), NF90_UNLIMITED, new_dimid) ! does this do nothing if the dimension already exists?
         time_dimension_exists = .true.
      else
         ret = nf90_def_dim(ncfile_out, dim_names(i, j, dom), dimensions_and_lengths(i, j+1, dom), new_dimid) ! does this do nothing if the dimension already exists?
      endif
      if(ret == NF90_NOERR) then ! successfully created, store this dimenion id
         where (dimIds == dimIds(i, j, dom))
            copy_dimIds = new_dimid
         end where
      endif

   enddo
enddo

if ((.not. time_dimension_exists) .and. (time_unlimited)) then ! create unlimlited dimension time
   ret = nf90_def_dim(ncfile_out, 'Time', NF90_UNLIMITED, new_dimid) !> @todo Case sensitive?
   call nc_check(ret, 'create_state_output', 'creating time as the unlimited dimension')

   call add_time_unlimited(new_dimid)

endif

! overwrite dimIds
dimIds = copy_dimIds

! define variables
do i = 1, num_state_variables ! loop around state variables
   ! double or single precision?
   ndims = dimensions_and_lengths(i, 1, dom)

   if (single_precision_output) then
      xtype = nf90_real
   else ! write output that is the precision of filter
      if (r8 == digits12) then ! datasize = MPI_REAL8  ! What should we be writing?
         xtype = nf90_double
      else
         xtype = nf90_real
      endif
   endif

   ret = nf90_def_var(ncfile_out, trim(global_variable_names(i)), xtype=xtype, dimids=dimIds(i, 1:ndims, dom), varid=new_varid)
   call nc_check(ret, 'create_state_output', 'defining variable')
      variable_ids(i, dom) = new_varid

enddo

ret = nf90_enddef(ncfile_out)
call nc_check(ret, 'create_state_output', 'end define mode')

end subroutine create_state_output

!-------------------------------------------------
!> Write variables from start_var to end_var
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

   var_size = variable_sizes(i, domain)

   ! number of dimensions and length of each
   allocate(dims(dimensions_and_lengths(i, 1, domain)))
   dims = dimensions_and_lengths(i, 2:dimensions_and_lengths(i,1, domain) + 1, domain)

   ret = nf90_inq_varid(ncfile_out, global_variable_names(i), var_id)
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

   var_size = variable_sizes(i, domain)

   ! check whether you have to do anything to the variable, clamp or fail
   if (do_clamp_or_fail(i, domain)) then
      call clamp_or_fail_it(i, domain, var_block(start_in_var_block:start_in_var_block+var_size-1))
   endif

   ! number of dimensions and length of each
   allocate(dims(dimensions_and_lengths(i, 1, domain)))
   dims = dimensions_and_lengths(i, 2:dimensions_and_lengths(i,1, domain) + 1, domain)

   ret = nf90_inq_varid(ncfile_out, global_variable_names(i), var_id)
   call nc_check(ret, 'write_variables_clamp', 'getting variable id')

   ret = nf90_put_var(ncfile_out, var_id, var_block(start_in_var_block:start_in_var_block+var_size-1), count=dims)
   call nc_check(ret, 'write_variables_clamp', 'writing')
   start_in_var_block = start_in_var_block + var_size

   deallocate(dims)

enddo

end subroutine write_variables_clamp


!-------------------------------------------------------
!> Adding space for an unlimited dimension in the dimesion arrays
!> The unlimited dimension needs to be last in the list for def_var
subroutine add_time_unlimited(unlimited_dimId)

integer, intent(in)  :: unlimited_dimId

integer :: i !> loop variable

! add a dimension
dimensions_and_lengths(:, 1, :) = dimensions_and_lengths(:, 1, :) + 1

do i = 1, num_state_variables
   dimensions_and_lengths(i, dimensions_and_lengths(i, 1, 1) +1, :) = 1  ! unlimited dimension is length 1?
   copy_dimIds(i, dimensions_and_lengths(i, 1, 1), :) = unlimited_dimId
enddo

end subroutine add_time_unlimited

!-------------------------------------------------------
end module state_vector_io_mod
