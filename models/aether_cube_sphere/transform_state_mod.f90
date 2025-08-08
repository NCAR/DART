module transform_state_mod

use netcdf
use types_mod,            only : r4, r8, varnamelength, DEG2RAD
use netcdf_utilities_mod, only : nc_open_file_readonly, nc_open_file_readwrite, nc_close_file, &
                                 nc_create_file, nc_end_define_mode
use utilities_mod,        only : open_file, close_file, find_namelist_in_file, &
                                 check_namelist_read, error_handler, E_ERR, string_to_integer

! TEMPORARY USE OF MODEL MODE UNTIL GEOMETRY IS MOVED TO ITS OWN MODULE
use model_mod,           only : static_init_model, lat_lon_to_col_index

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

character(len=4) :: restart_ensemble_member, dart_ensemble_member

type :: file_type
   character(len=256) :: file_path
   integer            :: ncid, ncstatus, unlimitedDimId, nDimensions, nVariables, nAttributes, formatNum
end type file_type

type(file_type), allocatable :: block_files(:)
type(file_type)              :: filter_input_file, filter_output_file

integer            :: nblocks, nhalos
character(len=256) :: restart_file_prefix, restart_file_middle, restart_file_suffix, &
                      filter_input_prefix, filter_input_suffix, filter_output_prefix, &
                      filter_output_suffix
namelist /transform_state_nml/ nblocks, nhalos, restart_file_prefix, restart_file_middle, &
                               restart_file_suffix, filter_input_prefix, filter_input_suffix, &
                               filter_output_prefix, filter_output_suffix

character(len=256) :: restart_directory, grid_directory, filter_directory
namelist /directory_nml/ restart_directory, grid_directory, filter_directory

contains

!---------------------------------------------------------------

subroutine initialize_transform_state_mod()

restart_ensemble_member = get_ensemble_member_from_command_line()
dart_ensemble_member = zero_fill(integer_to_string(string_to_integer(restart_ensemble_member)+1), 4)

call find_namelist_in_file('input.nml', 'transform_state_nml', iunit)
read(iunit, nml = transform_state_nml, iostat = io)
call check_namelist_read(iunit, io, 'transform_state_nml')

call find_namelist_in_file('input.nml', 'directory_nml', iunit)
read(iunit, nml = directory_nml, iostat = io)
call check_namelist_read(iunit, io, 'directory_nml')

block_files = assign_block_files_array(nblocks, restart_ensemble_member, restart_directory, &
   restart_file_prefix, restart_file_middle, restart_file_suffix)

end subroutine initialize_transform_state_mod

!---------------------------------------------------------------

subroutine finalize_transform_state_mod()
   
integer :: iblock

! Close all of the files

do iblock = 1, nblocks
   call nc_close_file(block_files(iblock)%ncid)
end do

end subroutine finalize_transform_state_mod

!---------------------------------------------------------------

subroutine model_to_dart()

! JLA
integer :: iblock, dimid, length, ncols, dart_dimid(3), varid, xtype, nDimensions, nAtts
integer :: lat_varid, lon_varid, i, j, k, txs, tys
integer :: ntimes(nblocks), nzs(nblocks), nxs(nblocks), nys(nblocks)
integer :: haloed_nxs(nblocks), haloed_nys(nblocks)
integer :: dimids(NF90_MAX_VAR_DIMS)
real(r8) :: blat, blon
character(len=NF90_MAX_NAME) :: name, attribute
integer, allocatable :: dart_varids(:), col_index(:, :)
real(r8), allocatable :: block_lats(:, :, :), block_lons(:, :, :)

do iblock = 1, nblocks
   ! Open the block files, read only
   block_files(iblock)%ncid = nc_open_file_readonly(block_files(iblock)%file_path)
   ! There doesn't seem to be a helper procedure corresponding to nf90_inquire in
   ! netcdf_utilities_mod so this uses the external function directly from the netcdf library
   block_files(iblock)%ncstatus = nf90_inquire(block_files(iblock)%ncid, &
      block_files(iblock)%nDimensions, block_files(iblock)%nVariables, &
      block_files(iblock)%nAttributes,  block_files(iblock)%unlimitedDimId, &
      block_files(iblock)%formatNum)

   ! Allow the files to be of different size, but all must have the same number of halos from namelist

   ! Loop through each of the dimensions to find the metadata values
   do dimid = 1, block_files(iblock)%nDimensions
      ! There doesn't seem to be a helper procedure corresponding to nf90_inquire_dimension that
      ! assigns name and length in netcdf_utilities_mod so this uses the external function
      ! directly from the netcdf library
      block_files(iblock)%ncstatus = nf90_inquire_dimension(block_files(iblock)%ncid, dimid, name, length)

      if (trim(name) == 'time') then
         !JLA Error if more than one time???
         ntimes(iblock) = length
      else if (trim(name) == 'x') then
         nxs(iblock)         = length-2*nhalos
         haloed_nxs(iblock)  = length
      else if (trim(name) == 'y') then
         nys(iblock)        = length-2*nhalos
         haloed_nys(iblock) = length
      else if (trim(name) == 'z') then
         nzs(iblock) = length
      end if
   end do
   write(*, *) 'nxs, nys, nzs , nVars', iblock, nxs(iblock), nys(iblock), nzs(iblock), block_files(iblock)%nVariables

   ! Create temporary storage for this blocks variables (storage inefficient) 
end do

! Do some consistency checks to make sure the files have the same number of variables, levels and times
if(any(ntimes - ntimes(1) .ne. 0)) then
   write(*, *) 'inconsistent ntimes'
   stop
endif
if(any(nzs - nzs(1) .ne. 0)) then
   write(*, *) 'inconsistent ntimes'
   stop
endif
if(any(block_files(:)%nVariables - block_files(1)%nVariables .ne. 0)) then
   write(*, *) 'inconsistent ntimes'
   stop
endif

! Space for identifiers for variables
allocate(dart_varids(block_files(1)%nVariables))

! Compute the number of columns (without haloes)
ncols = sum(nxs(1:nblocks) * nys(1:nblocks))

! Initialize the filter input file that will be created
filter_input_file = assign_filter_file(dart_ensemble_member, filter_directory, &
   filter_input_prefix, '.nc')
! Create the filter file
filter_input_file%ncid = nc_create_file(filter_input_file%file_path)

! Create dimensions in filter_input_file; save for use during variable definition
filter_input_file%ncstatus = nf90_def_dim(filter_input_file%ncid, 'time', NF90_UNLIMITED, dart_dimid(3))
filter_input_file%ncstatus = nf90_def_dim(filter_input_file%ncid, 'z',    nzs(1),         dart_dimid(2))
filter_input_file%ncstatus = nf90_def_dim(filter_input_file%ncid, 'col',  ncols,          dart_dimid(1))

! The filter_input_file is still in define mode. Create all of the variables before entering data mode.
do varid = 1, block_files(1)%nVariables
   block_files(1)%ncstatus = nf90_inquire_variable(block_files(1)%ncid, varid, name, xtype, nDimensions, dimids, nAtts)
   if (trim(name) == 'time') then
      filter_input_file%ncstatus = nf90_def_var(filter_input_file%ncid, name, xtype, dart_dimid(3), dart_varids(varid))
   else if (trim(name) == 'z') then
      ! Rename the 'z' variable as 'alt' so there isn't a dimension and a variable with the same name
      filter_input_file%ncstatus = nf90_def_var(filter_input_file%ncid, 'alt', xtype, dart_dimid(2), dart_varids(varid))
   else if (trim(name) == 'lat') then
      filter_input_file%ncstatus = nf90_def_var(filter_input_file%ncid, name, xtype, dart_dimid(1), dart_varids(varid))
      lat_varid = varid
   else if (trim(name) == 'lon') then
      filter_input_file%ncstatus = nf90_def_var(filter_input_file%ncid, name, xtype, dart_dimid(1), dart_varids(varid))
      lon_varid = varid
   else
      filter_input_file%ncstatus = nf90_def_var(filter_input_file%ncid, name, xtype, dart_dimid, dart_varids(varid))
   end if

   ! In the block files, time does not have units
   if (trim(name) /= 'time') then
      block_files(1)%ncstatus = nf90_get_att(block_files(1)%ncid, varid, 'units', attribute)
      filter_input_file%ncstatus = nf90_put_att(filter_input_file%ncid, dart_varids(varid), 'units', attribute)
   end if

   ! In the block files, only lon, lat and z have long_name
   if ((trim(name) == 'lon') .or. (trim(name) == 'lat') .or. (trim(name) == 'z')) then
      block_files(1)%ncstatus = nf90_get_att(block_files(1)%ncid, varid, 'long_name', attribute)
      filter_input_file%ncstatus = nf90_put_att(filter_input_file%ncid, dart_varids(varid), 'long_name', attribute)
   end if
end do

! End of define mode, ready to add data
call nc_end_define_mode(filter_input_file%ncid)


! Temporary until grid geometry routines have their own module
call static_init_model

! Loop through the blocks to get latitude and longitude arrays
k = 0
do iblock = 1, nblocks
   ! Get the latitude and longtidue values 
   ! JLA BE SURE I UNNDERSTAND THE X Y ORDER ON THESE
   txs = nxs(iblock) + 2*nhalos
   tys = nys(iblock) + 2*nhalos
   ! Careful of the storage order, z, y, x is read by get_var
   ! The array of column_index given y and x; storage order consistent with blocks
   allocate(block_lats(nzs(1), tys, txs), block_lons(nzs(1), tys, txs))
   allocate(col_index(nys(iblock), nxs(iblock)))

   block_files(iblock)%ncstatus = nf90_get_var(block_files(iblock)%ncid, lat_varid, block_lats)
   block_files(iblock)%ncstatus = nf90_get_var(block_files(iblock)%ncid, lon_varid, block_lons)

   ! Need to avoid initializing model_mod, should move grid routines to a separate module
   do i = 1, nxs(iblock)
      do j = 1, nys(iblock)
         blat = block_lats(1, nhalos + j, nhalos + i);
         blon = block_lons(1, nhalos + j, nhalos + i);
         col_index(j, i) = lat_lon_to_col_index(DEG2RAD*blat, DEG2RAD*blon)
         !write(*, *) i, j, blat, blon, col_index(j, i)
      end do
   end do

   ! Print them out in the lon lat order for the non halo as a check
   do j = 1, nys(iblock)
      do i = 1, nxs(iblock)
k = k+1
         write(*, *) i, j, col_index(j, i)
if(k .ne. col_index(j, i)) then
   write(*, *) 'k col ', j, i, k, col_index(j, i)
   write(*, *) 'a mess'
   stop
endif
      end do
   end do

   deallocate(block_lats, block_lons)
   deallocate(col_index)
end do

stop

end subroutine model_to_dart

!---------------------------------------------------------------

subroutine model_to_dart_orig()

integer :: iblock, dimid, dart_dimid, ix, iy, iz, icol, varid
integer :: length, xtype, nDimensions, nAtts, ntimes, nzs, nxs_per_block, nys_per_block
integer :: truncated_nxs_per_block, truncated_nys_per_block, total_truncated_ncols
integer :: time_lev_col_dims(3)
integer :: dimids(NF90_MAX_VAR_DIMS)
integer, allocatable :: dart_varids(:)

character(len=NF90_MAX_NAME) :: name
character(len=NF90_MAX_NAME) :: attribute


! The time variable in the block files is a double
real(r8), allocatable, dimension(:) :: time_array
! The other variables are floats JLA Can't we use R8?
real(r4), allocatable :: block_array(:, :, :)
real(r4), allocatable :: spatial_array(:)
real(r4), allocatable :: variable_array(:, :, :)

filter_input_file = assign_filter_file(dart_ensemble_member, filter_directory, &
   filter_input_prefix, filter_input_suffix)

! The block files are read only
do iblock = 1, nblocks
   block_files(iblock)%ncid = nc_open_file_readonly(block_files(iblock)%file_path)
end do

! The dart file is create
filter_input_file%ncid = nc_create_file(filter_input_file%file_path)

! The first set of nested loops iterates through all of the block files and all of the dimensions
! of each block file and counts the lengths of each dimension.

do iblock = 1, nblocks
   ! There doesn't seem to be a helper procedure corresponding to nf90_inquire in
   ! netcdf_utilities_mod so this uses the external function directly from the netcdf library
!JLA DO I REALLY NEED A LOOP OVER IBLOCK HERE OR CAN I JUST GET THIS FROM THE FIRST BLOCK
   block_files(iblock)%ncstatus = nf90_inquire(block_files(iblock)%ncid, &
      block_files(iblock)%nDimensions, block_files(iblock)%nVariables, &
      block_files(iblock)%nAttributes,  block_files(iblock)%unlimitedDimId, &
      block_files(iblock)%formatNum)
      
   if (iblock == 1) then
      do dimid = 1, block_files(iblock)%nDimensions
         ! There doesn't seem to be a helper procedure corresponding to nf90_inquire_dimension that
         ! assigns name and length in netcdf_utilities_mod so this uses the external function
         ! directly from the netcdf library
         block_files(iblock)%ncstatus = nf90_inquire_dimension(block_files(iblock)%ncid, dimid, name, length)

         if (trim(name) == 'time') then
            !JLA Error if more than one time???
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

! Could do some gross error checks on things we are taking for granted in these block files at this point JLA

! JLA: Need to check with Aaron to make sure blocks are homogenous, maybe just do each block separately
total_truncated_ncols = truncated_nxs_per_block*truncated_nys_per_block*nblocks

! All of the lengths have been counted properly, create each dimension in the filter_input_file and save
! the dimensions to the time_x_y_z and x_y_z arrays used during variable definition
filter_input_file%ncstatus = nf90_def_dim(filter_input_file%ncid, 'time', NF90_UNLIMITED, dart_dimid)
time_lev_col_dims(3) = dart_dimid

filter_input_file%ncstatus = nf90_def_dim(filter_input_file%ncid, 'z', nzs, dart_dimid)
time_lev_col_dims(2) = dart_dimid

filter_input_file%ncstatus = nf90_def_dim(filter_input_file%ncid, 'col', total_truncated_ncols, dart_dimid)
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

! The filter_input_file is still in define mode. Create all of the variables before entering data mode.
do varid = 1, block_files(1)%nVariables
   block_files(1)%ncstatus = nf90_inquire_variable(block_files(1)%ncid, varid, name, xtype, nDimensions, dimids, nAtts)
   if (trim(name) == 'time') then
      filter_input_file%ncstatus = nf90_def_var(filter_input_file%ncid, name, xtype, time_lev_col_dims(3), dart_varids(varid))
   else if (trim(name) == 'z') then
      ! Rename the 'z' variable as 'alt' so there isn't a dimension and a variable with the same name
      filter_input_file%ncstatus = nf90_def_var(filter_input_file%ncid, 'alt', xtype, time_lev_col_dims(2), dart_varids(varid))
   else if ((trim(name) == 'lon') .or. (trim(name) == 'lat')) then
      filter_input_file%ncstatus = nf90_def_var(filter_input_file%ncid, name, xtype, time_lev_col_dims(1), dart_varids(varid))
   else
      filter_input_file%ncstatus = nf90_def_var(filter_input_file%ncid, name, xtype, time_lev_col_dims, dart_varids(varid))
   end if

   ! Add attribute from block file except for time which has none
   if (trim(name) /= 'time') then
      block_files(1)%ncstatus = nf90_get_att(block_files(1)%ncid, varid, 'units', attribute)
      filter_input_file%ncstatus = nf90_put_att(filter_input_file%ncid, dart_varids(varid), 'units', attribute)
   end if

   ! In the block files, only lon, lat and z have long_name
   if ((trim(name) == 'lon') .or. (trim(name) == 'lat') .or. (trim(name) == 'z')) then
      block_files(1)%ncstatus = nf90_get_att(block_files(1)%ncid, varid, 'long_name', attribute)
      filter_input_file%ncstatus = nf90_put_att(filter_input_file%ncid, dart_varids(varid), 'long_name', attribute)
   end if

   ! print *, 'name: ' // name
   ! print *, 'dart_varids(varid): ' // integer_to_string(dart_varids(varid))

end do

call nc_end_define_mode(filter_input_file%ncid)

! The second set of nested loops has a different loop order. The outer loop is all of the
! variables while the inner loop is all of the blocks. The order is switched because all of the
! ncid pointers to each of the block files have already been assigned and it is more
! straightforward to assign all of the elements in the variable arrays if the blocks are the
! inner loop.

do varid = 1, block_files(1)%nVariables
   icol = 0
   do iblock = 1, nblocks
      block_files(iblock)%ncstatus = nf90_inquire_variable(block_files(iblock)%ncid, &
         varid, name, xtype, nDimensions, dimids, nAtts)
      
      if (trim(name) == 'time') then
         ! This is a 1-D time array
         if (iblock == 1) then
            block_files(iblock)%ncstatus = nf90_get_var(block_files(iblock)%ncid, varid, time_array)
            filter_input_file%ncstatus = nf90_put_var(filter_input_file%ncid, &
               dart_varids(varid), time_array)
         end if
      else if (trim(name) == 'z') then
         if (iblock == 1) then
            block_files(iblock)%ncstatus = nf90_get_var(block_files(iblock)%ncid, varid, block_array)
            filter_input_file%ncstatus = nf90_put_var(filter_input_file%ncid, &
               dart_varids(varid), block_array(:,1,1))
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
               filter_input_file%ncstatus = nf90_put_var(filter_input_file%ncid, &
                  dart_varids(varid), spatial_array)
            end if
         else
            ! This is one of the other non-spatial variables
            
            do iy = 1, truncated_nys_per_block
               do ix = 1, truncated_nxs_per_block
                  icol = icol + 1
                  do iz = 1, nzs
                     variable_array(icol, iz, 1) = block_array(iz, nhalos+iy, nhalos+ix)
                  end do
               end do
            end do

            if (iblock == nblocks) then
               filter_input_file%ncstatus = nf90_put_var(filter_input_file%ncid, &
                  dart_varids(varid), variable_array)
            end if
         end if
      end if
   end do
end do

call nc_close_file(filter_input_file%ncid)

end subroutine model_to_dart_orig

!---------------------------------------------------------------

subroutine dart_to_model()

integer :: iblock, icol, ix, iy, iz, ntimes, nxs, nys, nzs, ncols
integer :: filter_varid, block_varid, dimid, filter_xtype, filter_nDimensions, filter_nAtts
integer :: filter_dimids(NF90_MAX_VAR_DIMS)
character(len=NF90_MAX_NAME) :: filter_name, block_name
real(r4), allocatable :: filter_array(:, :, :), block_array(:, :, :)

filter_output_file = assign_filter_file(dart_ensemble_member, filter_directory, &
   filter_output_prefix, filter_output_suffix)

! The block files are read/write
do iblock = 1, nblocks
   block_files(iblock)%ncid = nc_open_file_readwrite(block_files(iblock)%file_path)
end do
! The dart file is read only
filter_output_file%ncid = nc_open_file_readonly(filter_output_file%file_path)

! Get the variables list from the filter_output_file

filter_output_file%ncstatus = nf90_inquire(filter_output_file%ncid, &
   filter_output_file%nDimensions, filter_output_file%nVariables, &
   filter_output_file%nAttributes, filter_output_file%unlimitedDimId, &
   filter_output_file%formatNum)

filter_output_file%ncstatus = nf90_inq_dimid(filter_output_file%ncid, 'time', dimid)
filter_output_file%ncstatus = nf90_inquire_dimension(filter_output_file%ncid, dimid, &
    filter_name, ntimes)

filter_output_file%ncstatus = nf90_inq_dimid(filter_output_file%ncid, 'z', dimid)
filter_output_file%ncstatus = nf90_inquire_dimension(filter_output_file%ncid, dimid, &
    filter_name, nzs)

filter_output_file%ncstatus = nf90_inq_dimid(filter_output_file%ncid, 'col', dimid)
filter_output_file%ncstatus = nf90_inquire_dimension(filter_output_file%ncid, dimid, &
    filter_name, ncols)

allocate(filter_array(ncols, nzs, ntimes))
filter_array(:, :, :) = 0
   
! We need full blocks from the block files
do filter_varid = 1, filter_output_file%nVariables
   icol = 0
   filter_output_file%ncstatus = nf90_inquire_variable(filter_output_file%ncid, &
      filter_varid, filter_name, filter_xtype, filter_nDimensions, filter_dimids, &
      filter_nAtts)
      
   if (filter_name /= 'time') then
      filter_output_file%ncstatus = nf90_get_var(filter_output_file%ncid, &
         filter_varid, filter_array)
         
      do iblock = 1, nblocks
         if (filter_varid == 1 .and. iblock == 1) then
            block_files(iblock)%ncstatus = nf90_inq_dimid(block_files(iblock)%ncid, &
               'x', dimid)
            block_files(iblock)%ncstatus = nf90_inquire_dimension(block_files(iblock)%ncid, &
               dimid, block_name, nxs)

            block_files(iblock)%ncstatus = nf90_inq_dimid(block_files(iblock)%ncid, &
               'y', dimid)
            block_files(iblock)%ncstatus = nf90_inquire_dimension(block_files(iblock)%ncid, &
                dimid, block_name, nys)

            allocate(block_array(nzs, nys, nxs))
            block_array(:, :, :) = 0
         end if

         block_files(iblock)%ncstatus = nf90_inq_varid(block_files(iblock)%ncid, &
            filter_name, block_varid)
         block_files(iblock)%ncstatus = nf90_get_var(block_files(iblock)%ncid, &
            block_varid, block_array)

         do iy = 1, nys-2*nhalos
            do ix = 1, nxs-2*nhalos
               icol = icol + 1
               do iz = 1, nzs
                  block_array(iz, nhalos+iy, nhalos+ix) = filter_array(icol, iz, 1)
               end do
            end do
         end do

         block_files(iblock)%ncstatus = nf90_put_var(block_files(iblock)%ncid, &
            block_varid, block_array)
         print *, block_files(iblock)%ncstatus
      
      end do
   end if
end do
   
call nc_close_file(filter_output_file%ncid)

end subroutine dart_to_model

!---------------------------------------------------------------

function get_ensemble_member_from_command_line() result(ensemble_member)
! Calls Fortran intrinsic subroutine get_command_argument and returns
! a string with four characters

character(len=4) :: ensemble_member
integer          :: nargs

nargs = command_argument_count()

if (nargs /= 1) then
   call error_handler(E_ERR, 'get_ensemble_member_from_command_line', &
      'ensemble member must be passed as a command line argument')
end if

call get_command_argument(1, ensemble_member)

end function get_ensemble_member_from_command_line

!---------------------------------------------------------------

function assign_block_files_array(nblocks, ensemble_member, restart_directory, &
   restart_file_prefix, restart_file_middle, restart_file_suffix) result(block_files)
   
integer,          intent(in) :: nblocks
character(len=4), intent(in)  :: ensemble_member
character(len=*), intent(in) :: restart_directory
character(len=*), intent(in) :: restart_file_prefix
character(len=*), intent(in) :: restart_file_middle
character(len=*), intent(in) :: restart_file_suffix
type(file_type), allocatable :: block_files(:)

character(len=4) :: block_name
integer          :: iblock

allocate(block_files(nblocks))

do iblock = 1, nblocks
   block_name = zero_fill(integer_to_string(iblock-1), 4)
   block_files(iblock)%file_path = trim(restart_directory) // &
      trim(restart_file_prefix) // ensemble_member // trim(restart_file_middle) // &
      block_name // trim(restart_file_suffix)
end do

end function assign_block_files_array

!---------------------------------------------------------------

function assign_filter_file(ensemble_member, filter_directory, filter_input_prefix, &
   filter_input_suffix) result(filter_file)
   
character(len=4), intent(in) :: ensemble_member
character(len=*), intent(in) :: filter_directory
character(len=*), intent(in) :: filter_input_prefix
character(len=*), intent(in) :: filter_input_suffix
type(file_type)              :: filter_file
   
filter_file%file_path = trim(filter_directory) // trim(filter_input_prefix) // &
   ensemble_member // trim(filter_input_suffix)

end function assign_filter_file

!---------------------------------------------------------------

function integer_to_string(int) result(string)

integer, intent(in)          :: int
character(len=varnamelength) :: string

write(string,'(I0)') int
string = trim(string)

end function integer_to_string

!---------------------------------------------------------------

function zero_fill(string, desired_length) result(filled_string)

character(len=*), intent(in) :: string
integer,          intent(in) :: desired_length

integer :: length_of_string
integer :: string_index, difference_of_string_lengths
character(len=varnamelength) :: filled_string

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

!---------------------------------------------------------------

end module transform_state_mod
