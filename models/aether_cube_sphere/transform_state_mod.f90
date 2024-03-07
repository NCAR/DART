module transform_state_mod

use netcdf
use types_mod,            only : r4, r8, varnamelength
use netcdf_utilities_mod, only : nc_open_file_readonly, nc_open_file_readwrite, nc_close_file, &
                                 nc_create_file, nc_end_define_mode
use utilities_mod,        only : open_file, close_file, find_namelist_in_file, &
                                 check_namelist_read, error_handler, E_ERR

implicit none

character(len=4) :: ensemble_member
integer :: nblocks, nhalos
character(len=256) :: restart_directory, restart_file_prefix, restart_file_middle, & 
                      restart_file_suffix, dart_directory, dart_file_prefix, dart_file_suffix

type :: file_type
character(len=256) :: file_path
integer :: ncid, ncstatus, unlimitedDimId, nDimensions, nVariables, nAttributes, formatNum
end type file_type

type(file_type), allocatable, dimension(:) :: block_files
type(file_type) :: dart_file

namelist /transform_state_nml/ &
   nblocks, &
   nhalos, &
   restart_directory, &
   restart_file_prefix, &
   restart_file_middle, &
   restart_file_suffix, &
   dart_directory, &
   dart_file_prefix, &
   dart_file_suffix

contains

subroutine initialize_transform_state_mod()

   ensemble_member = get_ensemble_member_from_command_line()

   call read_namelist()

   block_files = assign_block_files_array(nblocks, ensemble_member, restart_directory, &
                                          restart_file_prefix, restart_file_middle, &
                                          restart_file_suffix)

   dart_file = assign_dart_file(ensemble_member, dart_directory, dart_file_prefix, dart_file_suffix)

end subroutine initialize_transform_state_mod

subroutine finalize_transform_state_mod()
   
   integer :: iblock

   ! Close all of the files

   do iblock = 1, nblocks
      call nc_close_file(block_files(iblock)%ncid)
   end do

   call nc_close_file(dart_file%ncid)

end subroutine finalize_transform_state_mod

subroutine model_to_dart()

   integer :: iblock
   integer :: dimid, dart_dimid, dart_varid
   integer :: varid
   character(len=NF90_MAX_NAME) :: name
   character(len=NF90_MAX_NAME) :: attribute
   integer :: length
   integer :: xtype, nDimensions, nAtts
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimids
   integer :: ntimes = 0
   integer :: nxs = 0
   integer :: nys = 0
   integer :: nzs = 0
   integer :: nxs_per_block, nys_per_block
   integer, dimension(nblocks+1) :: cumulative_nxs, cumulative_nys

   integer, dimension(4) :: time_x_y_z_dims
   integer, dimension(3) :: x_y_z_dims

   ! The time variable in the block files is a double
   real(r8), allocatable, dimension(:) :: time_array
   ! The other variables are floats
   real(r4), allocatable, dimension(:, :, :) :: block_array
   real(r4), allocatable, dimension(:, :, :) :: spatial_array
   real(r4), allocatable, dimension(:, :, :, :) :: variable_array

   cumulative_nxs(1) = 0
   cumulative_nys(1) = 0

   ! The block files are read only
   do iblock = 1, nblocks
      block_files(iblock)%ncid = nc_open_file_readonly(block_files(iblock)%file_path)
   end do

   ! The dart file is create
   dart_file%ncid = nc_create_file(dart_file%file_path)

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
      
      do dimid = 1, block_files(iblock)%nDimensions
         ! There doesn't seem to be a helper procedure corresponding to nf90_inquire_dimension that
         ! assigns name and length in netcdf_utilities_mod so this uses the external function
         ! directly from the netcdf library
         block_files(iblock)%ncstatus = nf90_inquire_dimension(block_files(iblock)%ncid, dimid, name, length)

         ! Incrementing the length of dimensions requires detailed logic detailed in 
         ! in the comment in each statement
         if ((trim(name) == 'time') .and. (iblock == 1)) then
            ! The time dimension should only be incremented by the length of the time 
            ! dimension in the first file because the block_files are all output for the same times.
            ntimes = length
         else if (trim(name) == 'x') then
            ! The lon dimension should only be incremented by a length after
            ! the halos are removed
            nxs = nxs + (length-2*nhalos)
            if (iblock == 1) then
               nxs_per_block = length
            end if
            cumulative_nxs(iblock+1) = nxs
         else if (trim(name) == 'y') then
            ! The lat dimension should only be incremented by a length after
            ! the halos are removed
            nys = nys + (length-2*nhalos)
            if (iblock == 1) then
               nys_per_block = length
            end if
            cumulative_nys(iblock+1) = nys
         else if ((trim(name) == 'z') .and. (iblock == 1)) then
            ! The z dimension should only be incremented by the length of the z
            ! dimension in the first file because the block_files are all output at the same height.
            nzs = length
         end if
      end do
   end do

   ! All of the lengths have been counted properly, create each dimension in the dart_file and save
   ! the dimensions to the time_x_y_z and x_y_z arrays used during variable definition
   dart_file%ncstatus = nf90_def_dim(dart_file%ncid, 'time', ntimes, dart_dimid)
   time_x_y_z_dims(4) = dart_dimid

   dart_file%ncstatus = nf90_def_dim(dart_file%ncid, 'x', cumulative_nxs(nblocks+1), dart_dimid)
   time_x_y_z_dims(3) = dart_dimid
   x_y_z_dims(3) = dart_dimid

   dart_file%ncstatus = nf90_def_dim(dart_file%ncid, 'y', cumulative_nys(nblocks+1), dart_dimid)
   time_x_y_z_dims(2) = dart_dimid
   x_y_z_dims(2) = dart_dimid

   dart_file%ncstatus = nf90_def_dim(dart_file%ncid, 'z', nzs, dart_dimid)
   time_x_y_z_dims(1) = dart_dimid
   x_y_z_dims(1) = dart_dimid
   
   ! Allocate all of the storage arrays
   allocate(time_array(ntimes))
   allocate(block_array(nzs, nys_per_block, nxs_per_block))
   allocate(spatial_array(nzs, cumulative_nys(nblocks+1), cumulative_nxs(nblocks+1)))
   allocate(variable_array(nzs, cumulative_nys(nblocks+1), cumulative_nxs(nblocks+1), ntimes))

   ! The dart_file is still in define mode. Create all of the variables before entering data mode.
   do varid = 1, block_files(1)%nVariables
      block_files(1)%ncstatus = nf90_inquire_variable(block_files(1)%ncid, varid, name, xtype, nDimensions, dimids, nAtts)
      if (trim(name) == 'time') then
         dart_file%ncstatus = nf90_def_var(dart_file%ncid, name, xtype, time_x_y_z_dims(4), dart_varid)
      else if (trim(name) == 'z') then
         ! Rename the 'z' variable as 'alt' so there isn't a dimension and a variable with the same name
         dart_file%ncstatus = nf90_def_var(dart_file%ncid, 'alt', xtype, x_y_z_dims, dart_varid)
      else if ((trim(name) == 'lon') .or. (trim(name) == 'lat')) then
         dart_file%ncstatus = nf90_def_var(dart_file%ncid, name, xtype, x_y_z_dims, dart_varid)
      else
         dart_file%ncstatus = nf90_def_var(dart_file%ncid, name, xtype, time_x_y_z_dims, dart_varid)
      end if

      ! In the block files, time does not have units
      if (trim(name) /= 'time') then
         block_files(iblock)%ncstatus = nf90_get_att(block_files(1)%ncid, varid, 'units', attribute)
         dart_file%ncstatus = nf90_put_att(dart_file%ncid, dart_varid, 'units', attribute)
      end if

      ! In the block files, only lon, lat and z have long_name
      if ((trim(name) == 'lon') .or. (trim(name) == 'lat') .or. (trim(name) == 'z')) then
         block_files(iblock)%ncstatus = nf90_get_att(block_files(1)%ncid, varid, 'long_name', attribute)
         dart_file%ncstatus = nf90_put_att(dart_file%ncid, dart_varid, 'long_name', attribute)
      end if
   end do

   call nc_end_define_mode(dart_file%ncid)

   ! The second set of nested loops has a different loop order. The outer loop is all of the
   ! variables while the inner loop is all of the blocks. The order is switched because all of the
   ! ncid pointers to each of the block files have already been assigned and it is more
   ! straightforward to assign all of the elements in the variable arrays if the blocks are the
   ! inner loop.

   do varid = 1, block_files(1)%nVariables
      do iblock = 1, nblocks

         block_files(iblock)%ncstatus = nf90_inquire_variable(block_files(iblock)%ncid, varid, name, xtype, nDimensions, dimids, nAtts)

         if (trim(name) == 'time') then
            ! This is a 1-D time array
            if (iblock == 1) then
               block_files(iblock)%ncstatus = nf90_get_var(block_files(iblock)%ncid, varid, time_array)
            end if
   
            if (iblock == nblocks) then
               dart_file%ncstatus = nf90_put_var(dart_file%ncid, varid, time_array)
            end if
         else
            ! All of the variables besides time can be read into the block array
            block_files(iblock)%ncstatus = nf90_get_var(block_files(iblock)%ncid, varid, block_array)


            if ((trim(name) == 'lon') .or. (trim(name) == 'lat') .or. (trim(name) == 'z')) then
               spatial_array(:, cumulative_nys(iblock)+1:cumulative_nys(iblock+1), cumulative_nxs(iblock)+1:cumulative_nxs(iblock+1)) = block_array(:, nhalos+1:nys_per_block-nhalos, nhalos+1:nxs_per_block-nhalos)
            
               if (iblock == nblocks) then
                  dart_file%ncstatus = nf90_put_var(dart_file%ncid, varid, spatial_array)
               end if
            else
               variable_array(:, cumulative_nys(iblock)+1:cumulative_nys(iblock+1), cumulative_nxs(iblock)+1:cumulative_nxs(iblock+1), 1) = block_array(:, nhalos+1:nys_per_block-nhalos, nhalos+1:nxs_per_block-nhalos)
               if (iblock == nblocks) then
                  dart_file%ncstatus = nf90_put_var(dart_file%ncid, varid, variable_array)
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
   dart_file%ncid = nc_open_file_readonly(dart_file%file_path)

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

subroutine read_namelist()
   integer :: io, iunit

   call find_namelist_in_file('input.nml', 'transform_state_nml', iunit)
   read(iunit, nml = transform_state_nml, iostat = io)
   call check_namelist_read(iunit, io, 'transform_state_nml')

end subroutine read_namelist

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
   character(len=4) :: cube_face
   integer :: iface

   allocate(block_files(nblocks))

   do iface = 1, nblocks
      cube_face = zero_fill(integer_to_string(iface-1), 4)
      block_files(iface)%file_path = trim(restart_directory) // trim(restart_file_prefix) // &
                               ensemble_member // trim(restart_file_middle) // cube_face // &
                               trim(restart_file_suffix)
   end do

end function assign_block_files_array

function assign_dart_file(ensemble_member, dart_directory, dart_file_prefix, dart_file_suffix) &
                          result(dart_file)
   
   character(len=4), intent(in) :: ensemble_member
   character(len=*), intent(in) :: dart_directory
   character(len=*), intent(in) :: dart_file_prefix
   character(len=*), intent(in) :: dart_file_suffix
   type(file_type)              :: dart_file
   
   dart_file%file_path = trim(dart_directory) // trim(dart_file_prefix) // ensemble_member // trim(dart_file_suffix)

end function assign_dart_file

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
