module transform_state_mod

use netcdf
use types_mod,            only : r4, r8, varnamelength
use netcdf_utilities_mod, only : nc_open_file_readonly, nc_open_file_readwrite, nc_close_file, &
                                 nc_create_file, nc_end_define_mode
use utilities_mod,        only : open_file, close_file, find_namelist_in_file, &
                                 check_namelist_read, error_handler, E_ERR

implicit none
private

public :: initialize_transform_state_mod, &
          finalize_transform_state_mod, &
          model_to_dart, &
          dart_to_model, &
          integer_to_string, &
          file_type, &
          zero_fill

integer  :: iunit, io

character(len=4) :: ensemble_member

type :: file_type
character(len=256) :: file_path
integer :: ncid, ncstatus, unlimitedDimId, nDimensions, nVariables, nAttributes, formatNum
end type file_type

type(file_type), allocatable, dimension(:) :: block_files
type(file_type) :: filter_file

integer :: nblocks, nhalos
character(len=256) :: restart_file_prefix, restart_file_middle, restart_file_suffix, &
                      filter_input_prefix, filter_input_suffix
namelist /transform_state_nml/ nblocks, nhalos, restart_file_prefix, restart_file_middle, &
                               restart_file_suffix, filter_input_prefix, filter_input_suffix

character(len=256) :: restart_directory, grid_directory, filter_directory
namelist /directory_nml/ restart_directory, grid_directory, filter_directory

contains

subroutine initialize_transform_state_mod()

   ensemble_member = get_ensemble_member_from_command_line()

   call find_namelist_in_file('input.nml', 'transform_state_nml', iunit)
   read(iunit, nml = transform_state_nml, iostat = io)
   call check_namelist_read(iunit, io, 'transform_state_nml')

   block_files = assign_block_files_array(nblocks, ensemble_member, restart_directory, &
                                          restart_file_prefix, restart_file_middle, &
                                          restart_file_suffix)

   filter_file = assign_filter_file(ensemble_member, filter_directory, filter_input_prefix, filter_input_suffix)

end subroutine initialize_transform_state_mod

subroutine finalize_transform_state_mod()
   
   integer :: iblock

   ! Close all of the files

   do iblock = 1, nblocks
      call nc_close_file(block_files(iblock)%ncid)
   end do

   call nc_close_file(filter_file%ncid)

end subroutine finalize_transform_state_mod

subroutine model_to_dart()

   integer :: iblock
   integer :: dimid, dart_dimid
   integer :: ix, iy, iz, icol
   integer :: varid
   character(len=NF90_MAX_NAME) :: name
   character(len=NF90_MAX_NAME) :: attribute
   integer :: length
   integer :: xtype, nDimensions, nAtts
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimids
   integer :: ntimes
   integer :: nxs_per_block, nys_per_block, truncated_nxs_per_block, truncated_nys_per_block, total_truncated_ncols
   integer :: nzs

   integer, dimension(3) :: time_lev_col_dims
   integer, allocatable, dimension (:) :: dart_varids

   ! The time variable in the block files is a double
   real(r8), allocatable, dimension(:) :: time_array
   ! The other variables are floats
   real(r4), allocatable, dimension(:, :, :) :: block_array
   real(r4), allocatable, dimension(:) :: spatial_array
   real(r4), allocatable, dimension(:, :, :) :: variable_array

   ! The block files are read only
   do iblock = 1, nblocks
      block_files(iblock)%ncid = nc_open_file_readonly(block_files(iblock)%file_path)
   end do

   ! The dart file is create
   filter_file%ncid = nc_create_file(filter_file%file_path)

   ! The first set of nested loops iterates through all of the block files and all of the dimensions
   ! of each block file and counts the lengths of each dimension.

   do iblock = 1, nblocks
      ! There doesn't seem to be a helper procedure corresponding to nf90_inquire in
      ! netcdf_utilities_mod so this uses the external function directly from the netcdf library
      block_files(iblock)%ncstatus = nf90_inquire(block_files(iblock)%ncid, &
                                                  block_files(iblock)%nDimensions, &
                                                  block_files(iblock)%nVariables, &
                                                  block_files(iblock)%nAttributes, &
                                                  block_files(iblock)%unlimitedDimId, &
                                                  block_files(iblock)%formatNum)
      
      if (iblock == 1) then
         do dimid = 1, block_files(iblock)%nDimensions
            ! There doesn't seem to be a helper procedure corresponding to nf90_inquire_dimension that
            ! assigns name and length in netcdf_utilities_mod so this uses the external function
            ! directly from the netcdf library
            block_files(iblock)%ncstatus = nf90_inquire_dimension(block_files(iblock)%ncid, dimid, name, length)

            if (trim(name) == 'time') then
               ntimes = length
            else if (trim(name) == 'x') then
               truncated_nxs_per_block = length-2*nhalos
               nxs_per_block = length
            else if (trim(name) == 'y') then
               truncated_nys_per_block = length-2*nhalos
               nys_per_block = length
            else if (trim(name) == 'z') then
               nzs = length
            end if
         end do
      end if
   end do

   total_truncated_ncols = truncated_nxs_per_block*truncated_nys_per_block*nblocks

   ! All of the lengths have been counted properly, create each dimension in the filter_file and save
   ! the dimensions to the time_x_y_z and x_y_z arrays used during variable definition
   filter_file%ncstatus = nf90_def_dim(filter_file%ncid, 'time', ntimes, dart_dimid)
   time_lev_col_dims(3) = dart_dimid

   filter_file%ncstatus = nf90_def_dim(filter_file%ncid, 'z', nzs, dart_dimid)
   time_lev_col_dims(2) = dart_dimid

   filter_file%ncstatus = nf90_def_dim(filter_file%ncid, 'col', total_truncated_ncols, dart_dimid)
   time_lev_col_dims(1) = dart_dimid

   ! Allocate all of the storage arrays
   allocate(time_array(ntimes))
   allocate(block_array(nzs, nys_per_block, nxs_per_block))
   allocate(spatial_array(total_truncated_ncols))
   allocate(variable_array(total_truncated_ncols, nzs, ntimes))
   allocate(dart_varids(block_files(1)%nVariables))

   block_array(:, :, :) = 0
   spatial_array(:) = 0
   variable_array(:, :, :) = 0

   ! The filter_file is still in define mode. Create all of the variables before entering data mode.
   do varid = 1, block_files(1)%nVariables
      block_files(1)%ncstatus = nf90_inquire_variable(block_files(1)%ncid, varid, name, xtype, nDimensions, dimids, nAtts)
      if (trim(name) == 'time') then
         filter_file%ncstatus = nf90_def_var(filter_file%ncid, name, xtype, time_lev_col_dims(3), dart_varids(varid))
      else if (trim(name) == 'z') then
         ! Rename the 'z' variable as 'alt' so there isn't a dimension and a variable with the same name
         filter_file%ncstatus = nf90_def_var(filter_file%ncid, 'alt', xtype, time_lev_col_dims(2), dart_varids(varid))
      else if ((trim(name) == 'lon') .or. (trim(name) == 'lat')) then
         filter_file%ncstatus = nf90_def_var(filter_file%ncid, name, xtype, time_lev_col_dims(1), dart_varids(varid))
      else
         filter_file%ncstatus = nf90_def_var(filter_file%ncid, name, xtype, time_lev_col_dims, dart_varids(varid))
      end if

      ! In the block files, time does not have units
      if (trim(name) /= 'time') then
         block_files(iblock)%ncstatus = nf90_get_att(block_files(1)%ncid, varid, 'units', attribute)
         filter_file%ncstatus = nf90_put_att(filter_file%ncid, dart_varids(varid), 'units', attribute)
      end if

      ! In the block files, only lon, lat and z have long_name
      if ((trim(name) == 'lon') .or. (trim(name) == 'lat') .or. (trim(name) == 'z')) then
         block_files(iblock)%ncstatus = nf90_get_att(block_files(1)%ncid, varid, 'long_name', attribute)
         filter_file%ncstatus = nf90_put_att(filter_file%ncid, dart_varids(varid), 'long_name', attribute)
      end if

      ! print *, 'name: ' // name
      ! print *, 'dart_varids(varid): ' // integer_to_string(dart_varids(varid))

   end do

   call nc_end_define_mode(filter_file%ncid)

   ! The second set of nested loops has a different loop order. The outer loop is all of the
   ! variables while the inner loop is all of the blocks. The order is switched because all of the
   ! ncid pointers to each of the block files have already been assigned and it is more
   ! straightforward to assign all of the elements in the variable arrays if the blocks are the
   ! inner loop.

   do varid = 1, block_files(1)%nVariables
      icol = 0
      do iblock = 1, nblocks
         block_files(iblock)%ncstatus = nf90_inquire_variable(block_files(iblock)%ncid, varid, name, xtype, nDimensions, dimids, nAtts)
         
         if (trim(name) == 'time') then
            ! This is a 1-D time array
            if (iblock == 1) then
               block_files(iblock)%ncstatus = nf90_get_var(block_files(iblock)%ncid, varid, time_array)
               filter_file%ncstatus = nf90_put_var(filter_file%ncid, dart_varids(varid), time_array)
            end if
         else if (trim(name) == 'z') then
            if (iblock == 1) then
               block_files(iblock)%ncstatus = nf90_get_var(block_files(iblock)%ncid, varid, block_array)
               filter_file%ncstatus = nf90_put_var(filter_file%ncid, dart_varids(varid), block_array(:,1,1))
            end if
         else
            ! All of the variables besides time can be read into the block array
            block_files(iblock)%ncstatus = nf90_get_var(block_files(iblock)%ncid, varid, block_array)
            
            if ((trim(name) == 'lon') .or. (trim(name) == 'lat')) then
               do iy = 1, truncated_nys_per_block
                  do ix = 1, truncated_nxs_per_block
                     icol = icol + 1
                     spatial_array(icol) = block_array(1, nhalos+iy, nhalos+ix)
                  end do
               end do
               
               if (iblock == nblocks) then
                  filter_file%ncstatus = nf90_put_var(filter_file%ncid, dart_varids(varid), spatial_array)
               end if
            else
               ! This is one of the other non-spatial variables
               do iz = 1, nzs
                  do iy = 1, truncated_nys_per_block
                     do ix = 1, truncated_nxs_per_block
                        icol = icol + 1
                        variable_array(icol, iz, 1) = block_array(iz, nhalos+iy, nhalos+ix)
                     end do
                  end do
               end do

               if (iblock == nblocks) then
                  filter_file%ncstatus = nf90_put_var(filter_file%ncid, dart_varids(varid), variable_array)
               end if
            end if
         end if
      end do
   end do

end subroutine model_to_dart

subroutine dart_to_model()

   integer :: iblock, idim

   ! The block files are read/write
   do iblock = 1, nblocks
      block_files(iblock)%ncid = nc_open_file_readwrite(block_files(iblock)%file_path)
   end do
   ! The dart file is read only
   filter_file%ncid = nc_open_file_readonly(filter_file%file_path)

end subroutine dart_to_model

function get_ensemble_member_from_command_line() result(ensemble_member)
   ! Calls Fortran intrinsic subroutine get_command_argument and returns
   ! a string with four characters

   character(len=4) :: ensemble_member
   integer :: nargs

   nargs = command_argument_count()

   if (nargs /= 1) then
      call error_handler(E_ERR, 'get_ensemble_member_from_command_line', &
                         'ensemble member must be passed as a command line argument')
   end if

   call get_command_argument(1, ensemble_member)

end function get_ensemble_member_from_command_line

function assign_block_files_array(nblocks, ensemble_member, restart_directory, &
                                  restart_file_prefix, restart_file_middle, restart_file_suffix) &
                                  result(block_files)
   
   integer, intent(in)            :: nblocks
   character(len=4), intent(in)   :: ensemble_member
   character(len=*), intent(in) :: restart_directory
   character(len=*), intent(in) :: restart_file_prefix
   character(len=*), intent(in) :: restart_file_middle
   character(len=*), intent(in) :: restart_file_suffix
   type(file_type), allocatable, dimension(:) :: block_files
   character(len=4) :: block_name
   integer :: iblock

   allocate(block_files(nblocks))

   do iblock = 1, nblocks
      block_name = zero_fill(integer_to_string(iblock-1), 4)
      block_files(iblock)%file_path = trim(restart_directory) // trim(restart_file_prefix) // &
                                      ensemble_member // trim(restart_file_middle) // &
                                      block_name // trim(restart_file_suffix)
   end do

end function assign_block_files_array

function assign_filter_file(ensemble_member, filter_directory, filter_input_prefix, filter_input_suffix) &
                          result(filter_file)
   
   character(len=4), intent(in) :: ensemble_member
   character(len=*), intent(in) :: filter_directory
   character(len=*), intent(in) :: filter_input_prefix
   character(len=*), intent(in) :: filter_input_suffix
   type(file_type)              :: filter_file
   
   filter_file%file_path = trim(filter_directory) // trim(filter_input_prefix) // ensemble_member // trim(filter_input_suffix)

end function assign_filter_file

function integer_to_string(int) result(string)

   integer, intent(in) :: int
   character(len=varnamelength) :: string

   write(string,'(I0)') int
   string = trim(string)

end function integer_to_string

function zero_fill(string, desired_length) result(filled_string)

   character(len=*), intent(in) :: string
   integer, intent(in) :: desired_length

   character(len=varnamelength) :: filled_string
   integer :: length_of_string
   integer :: string_index, difference_of_string_lengths

   filled_string = ''
   length_of_string = len_trim(string)
   difference_of_string_lengths = desired_length - length_of_string

   if (difference_of_string_lengths < 0) then
      print *, 'Error: input string is longer than the desired output string.'
      stop
   else if (difference_of_string_lengths > 0) then
      do string_index = 1, difference_of_string_lengths
         filled_string(string_index:string_index) = '0'
      end do
   end if

   filled_string(difference_of_string_lengths+1:desired_length) = trim(string)

end function zero_fill

end module transform_state_mod
