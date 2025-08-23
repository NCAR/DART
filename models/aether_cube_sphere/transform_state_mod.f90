module transform_state_mod

use netcdf
use types_mod,            only : r4, r8, varnamelength, DEG2RAD, RAD2DEG
use netcdf_utilities_mod, only : nc_open_file_readonly, nc_open_file_readwrite, nc_close_file, &
                                 nc_create_file, nc_end_define_mode
use utilities_mod,        only : open_file, close_file, find_namelist_in_file, &
                                 check_namelist_read, error_handler, E_ERR, string_to_integer

use cube_sphere_grid_tools_mod, only : lat_lon_to_col_index, get_grid_delta

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

type(file_type), allocatable :: block_files(:), grid_files(:)
type(file_type)              :: filter_input_file, filter_output_file

! It would be nice to get this information from the Aether input files, but that may not be possible
integer            :: np, nblocks, nhalos
character(len=256) :: restart_file_prefix, restart_file_middle, restart_file_suffix, &
                      filter_input_prefix, filter_input_suffix, filter_output_prefix, &
                      filter_output_suffix
namelist /transform_state_nml/ np, nblocks, nhalos, restart_file_prefix, restart_file_middle, &
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

grid_files = assign_grid_files_array(nblocks)

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

integer :: iblock, dimid, length, ncols, dart_dimid(3), varid, xtype, nDimensions, nAtts
integer :: lat_varid, lon_varid, ix, iy, iz, txs, tys, icol
integer :: ntimes(nblocks), nzs(nblocks), nxs(nblocks), nys(nblocks)
integer :: final_ntimes, final_nxs, final_nys, final_nzs
integer :: haloed_nxs(nblocks), haloed_nys(nblocks)
integer :: dimids(NF90_MAX_VAR_DIMS)
real(r8) :: blat, blon, cube_side, del, half_del
character(len=NF90_MAX_NAME) :: name, attribute
integer,  allocatable         :: dart_varids(:), col_index(:, :, :)
! The time variable in the block files is a double
real(r8), allocatable, dimension(:) :: time_array
! File for reading in variables from block file; These can probably be R4
real(r4), allocatable :: block_array(:, :, :), spatial_array(:), variable_array(:, :, :)
real(r4), allocatable :: block_lats(:, :, :), block_lons(:, :, :)

! Get grid spacing from number of points across each face
call get_grid_delta(np, del, half_del)

!==================================== get info from grid file block ===================

do iblock = 1, nblocks
   ! Open the block files, read only
   grid_files(iblock)%ncid = nc_open_file_readonly(grid_files(iblock)%file_path)
   ! There doesn't seem to be a helper procedure corresponding to nf90_inquire in
   ! netcdf_utilities_mod so this uses the external function directly from the netcdf library
   grid_files(iblock)%ncstatus = nf90_inquire(grid_files(iblock)%ncid, &
      grid_files(iblock)%nDimensions, grid_files(iblock)%nVariables, &
      grid_files(iblock)%nAttributes,  grid_files(iblock)%unlimitedDimId, &
      grid_files(iblock)%formatNum)

   ! The number of variables should be 4: longitude, latitude, altitude, time
   if(grid_files(iblock)%nVariables .ne. 4) then
      write(*, *) 'nunmber of vars in grid files should be 4', grid_files(iblock)%nVariables
      stop
   endif

   ! Allow the files to be of different size, but all must have the same number of halos from namelist

   ! Loop through each of the dimensions to find the metadata values
   do dimid = 1, grid_files(iblock)%nDimensions
      ! There doesn't seem to be a helper procedure corresponding to nf90_inquire_dimension that
      ! assigns name and length in netcdf_utilities_mod so this uses the external function
      ! directly from the netcdf library
      grid_files(iblock)%ncstatus = nf90_inquire_dimension(grid_files(iblock)%ncid, dimid, name, length)

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

   ! Create temporary storage for this blocks variables (storage inefficient) 
end do

! Do some consistency checks to make sure the files have the same number of variables, levels and times
if(any(ntimes - ntimes(1) .ne. 0)) then
   write(*, *) 'inconsistent ntimes'
   stop
endif
if(any(nzs - nzs(1) .ne. 0)) then
   write(*, *) 'inconsistent number of vertical levels'
   stop
endif

! Final consistent times, x, y and levels
final_ntimes = ntimes(1)
final_nxs = nxs(1)
final_nys = nys(1)
final_nzs = nzs(1)
write(*, *) 'times, nx, ny, nz', final_ntimes, final_nxs, final_nys, final_nzs

! Compute the number of columns (without haloes)
ncols = sum(nxs(1:nblocks) * nys(1:nblocks))
write(*, *) 'ncols is ', ncols

! Allocate ncols size temporary storage
allocate(spatial_array(ncols), variable_array(ncols, nzs(1), ntimes(1)))
spatial_array = 0.0_r8

!==================================== end of get info from grid file block ===================

! Start with block files, then switch to ions and neutrals

do iblock = 1, nblocks
   ! Open the block files, read only
   block_files(iblock)%ncid = nc_open_file_readonly(block_files(iblock)%file_path)
write(*, *) 'block files ', iblock, block_files(iblock)%file_path
   ! There doesn't seem to be a helper procedure corresponding to nf90_inquire in
   ! netcdf_utilities_mod so this uses the external function directly from the netcdf library
   block_files(iblock)%ncstatus = nf90_inquire(block_files(iblock)%ncid, &
      block_files(iblock)%nDimensions, block_files(iblock)%nVariables, &
      block_files(iblock)%nAttributes,  block_files(iblock)%unlimitedDimId, &
      block_files(iblock)%formatNum)

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
! Comparison is to grid files for times and levels, but just among block files for number of variables
if(any(ntimes - final_ntimes .ne. 0)) then
   write(*, *) 'inconsistent ntimes'
   stop
endif
if(any(nzs - final_nzs .ne. 0)) then
   write(*, *) 'inconsistent of vertical levels'
   stop
endif
if(any(block_files(:)%nVariables - block_files(1)%nVariables .ne. 0)) then
   write(*, *) 'inconsistent number of variables'
   stop
endif

! Space for identifiers for variables
allocate(dart_varids(block_files(1)%nVariables), time_array(ntimes(1)))

! Initialize the filter input file that will be created
filter_input_file = assign_filter_file(dart_ensemble_member, filter_directory, &
   filter_input_prefix, '.nc')
! Create the filter file
filter_input_file%ncid = nc_create_file(filter_input_file%file_path)

! Create dimensions in filter_input_file; save for use during variable definition
filter_input_file%ncstatus = nf90_def_dim(filter_input_file%ncid, 'time', NF90_UNLIMITED, dart_dimid(3))
filter_input_file%ncstatus = nf90_def_dim(filter_input_file%ncid, 'z',    nzs(1),         dart_dimid(2))
filter_input_file%ncstatus = nf90_def_dim(filter_input_file%ncid, 'col',  ncols,          dart_dimid(1))

! The ions and neutrals files have time and all their physical variables, but not latitude, longitude, or altitude

!=========================================================
! Write out the fields from the grid files first

! The filter_input_file is still in define mode. Create all of the variables before entering data mode.
do varid = 1, 4
   grid_files(1)%ncstatus = nf90_inquire_variable(grid_files(1)%ncid, varid, name, xtype, nDimensions, dimids, nAtts)
   if (trim(name) == 'time') then
      filter_input_file%ncstatus = nf90_def_var(filter_input_file%ncid, name, xtype, dart_dimid(3), dart_varids(varid))
   else if (trim(name) == 'Altitude') then
      ! Rename the 'z' variable as 'alt' so there isn't a dimension and a variable with the same name
      filter_input_file%ncstatus = nf90_def_var(filter_input_file%ncid, 'alt', xtype, dart_dimid(2), dart_varids(varid))
      filter_input_file%ncstatus = nf90_put_att(filter_input_file%ncid, dart_varids(varid), 'long_name', &
         'height above mean sea level')
      filter_input_file%ncstatus = nf90_put_att(filter_input_file%ncid, dart_varids(varid), 'units', 'crazym')
   else if (trim(name) == 'Latitude') then
      filter_input_file%ncstatus = nf90_def_var(filter_input_file%ncid, 'lat', xtype, dart_dimid(1), dart_varids(varid))
      filter_input_file%ncstatus = nf90_put_att(filter_input_file%ncid, dart_varids(varid), 'long_name', 'latitude')
      filter_input_file%ncstatus = nf90_put_att(filter_input_file%ncid, dart_varids(varid), 'units', 'degrees north')
      lat_varid = varid
   else if (trim(name) == 'Longitude') then
      filter_input_file%ncstatus = nf90_def_var(filter_input_file%ncid, 'lon', xtype, dart_dimid(1), dart_varids(varid))
      filter_input_file%ncstatus = nf90_put_att(filter_input_file%ncid, dart_varids(varid), 'long_name', 'longitude')
      filter_input_file%ncstatus = nf90_put_att(filter_input_file%ncid, dart_varids(varid), 'units', 'degrees east')
      lon_varid = varid
   else
      write(*, *) 'Unexpected variable name in grid file', trim(name)
      stop
   end if

   ! In the block files, time does not have units
   if (trim(name) /= 'time') then
      grid_files(1)%ncstatus = nf90_get_att(grid_files(1)%ncid, varid, 'units', attribute)
      filter_input_file%ncstatus = nf90_put_att(filter_input_file%ncid, dart_varids(varid), 'units', attribute)
   end if

end do

!=========================================================

! Now get the other data fields from the block files (soon to be ions and neutrals)

! The filter_input_file is still in define mode. Create all of the variables before entering data mode.
do varid = 1, block_files(1)%nVariables
   block_files(1)%ncstatus = nf90_inquire_variable(block_files(1)%ncid, varid, name, xtype, nDimensions, dimids, nAtts)
   ! Only time should occur once we switch to ions and neutrals files
   if (trim(name) == 'time' .or. trim(name) == 'z' .or. trim(name) == 'lat' .or. trim(name) == 'lon') then
      ! Time has already been incorporated
   else
      filter_input_file%ncstatus = nf90_def_var(filter_input_file%ncid, name, xtype, dart_dimid, dart_varids(varid))
   end if
end do

! End of define mode, ready to add data
call nc_end_define_mode(filter_input_file%ncid)

!=========== Block to get lat lon and alt from grid files =================

! The col_index array will keep track of mapping from x and y for each block to final columns
allocate(col_index(nblocks, maxval(nys), maxval(nxs)))

! Choose to minimize storage rather than redundant computation
! Will get full spatial field for one variable at a time
do varid = 1, 4

   ! Loop through all the blocks for this variable
   do iblock = 1, nblocks
      txs = nxs(iblock) + 2*nhalos
      tys = nys(iblock) + 2*nhalos
      ! Careful of the storage order, z, y, x is read by get_var
      ! Allocate storage for the latitude and longitude from the blocks
      allocate(block_lats(final_nzs, tys, txs), block_lons(final_nzs, tys, txs), &
         block_array(final_nzs, tys, txs))
     
      ! Get the latitude and longitude full arrays 
      grid_files(iblock)%ncstatus = nf90_get_var(grid_files(iblock)%ncid, lat_varid, block_lats)
      grid_files(iblock)%ncstatus = nf90_get_var(grid_files(iblock)%ncid, lon_varid, block_lons)
  
      ! Compute the col_index for each of the horizontal locations in this block 
      do ix = 1, nxs(iblock)
         do iy = 1, nys(iblock)
            blat = block_lats(1, nhalos + iy, nhalos + ix);
            blon = block_lons(1, nhalos + iy, nhalos + ix);
            col_index(iblock, iy, ix) = lat_lon_to_col_index(blat, blon, del, half_del, np)
         end do
      end do
 
      ! Get metadata for this variable from block file 
      grid_files(iblock)%ncstatus = nf90_inquire_variable(grid_files(iblock)%ncid, &
         varid, name, xtype, nDimensions, dimids, nAtts)

      ! This is a 1-D time array
      if (trim(name) == 'time') then
         ! Time must be the same in all files, so just deal with it from the first one
         if (iblock == 1) then
            grid_files(iblock)%ncstatus = nf90_get_var(grid_files(iblock)%ncid, varid, time_array)
            filter_input_file%ncstatus = nf90_put_var(filter_input_file%ncid, &
               dart_varids(varid), time_array)
         end if
      else if (trim(name) == 'Altitude') then
         ! Vertical coords also must be the same in all files, so just set from the first
         if (iblock == 1) then
            grid_files(iblock)%ncstatus = nf90_get_var(grid_files(iblock)%ncid, varid, block_array)
            filter_input_file%ncstatus = nf90_put_var(filter_input_file%ncid, &
        ! JLA REVIEW THE DART_VARIDS STUFF
               dart_varids(varid), block_array(:,1,1))
         endif
      else

         ! The lat, and lon variables can be read into the full 3Dblock array
         grid_files(iblock)%ncstatus = nf90_get_var(grid_files(iblock)%ncid, varid, block_array)
         
         if ((trim(name) == 'Longitude') .or. (trim(name) == 'Latitude')) then
            do iy = 1, nys(iblock)
               do ix = 1, nxs(iblock)
                  icol = col_index(iblock, iy, ix)
                  spatial_array(icol) = block_array(1, nhalos+iy, nhalos+ix) * RAD2DEG
write(*, *) 'iy, ix, icol, value ', iy, ix, icol, spatial_array(icol)
               end do
            end do

            if (iblock == nblocks) then
write(*, *) 'spatial array ', trim(name), spatial_array
               filter_input_file%ncstatus = nf90_put_var(filter_input_file%ncid, &
                  dart_varids(varid), spatial_array)
            end if
         else
            ! Error if it isn't one of the 4 fields
            write(*, *) 'Unknown variable in grid file ', trim(name)
            stop
         end if
      end if
 
      deallocate(block_lats, block_lons, block_array)
   end do
end do

!=========== End of Block to get lat lon and alt from grid files =================

! Choose to minimize storage rather than redundant computation
! Will get full spatial field for one variable at a time
do varid = 1, block_files(1)%nVariables

   ! Loop through all the blocks for this variable
   do iblock = 1, nblocks
      txs = nxs(iblock) + 2*nhalos
      tys = nys(iblock) + 2*nhalos
      ! Careful of the storage order, z, y, x is read by get_var
      ! Allocate storage for the fields
      allocate(block_array(final_nzs, tys, txs))
     
      ! Get metadata for this variable from block file 
      block_files(iblock)%ncstatus = nf90_inquire_variable(block_files(iblock)%ncid, &
         varid, name, xtype, nDimensions, dimids, nAtts)

      ! Time was taken care of already fron grid file
      if (trim(name) .ne. 'time') then

         ! All of the other variables can be read into the full 3Dblock array
         block_files(iblock)%ncstatus = nf90_get_var(block_files(iblock)%ncid, varid, block_array)
         
         ! This is one of the other non-spatial variables
!!!write(*, *) 'field name ', trim(name), varid, block_array(1, 1, 1)
         do iy = 1, nys(iblock)
            do ix = 1, nxs(iblock)
               icol = col_index(iblock, iy, ix)
               do iz = 1, nzs(1)
                  variable_array(icol, iz, 1) = block_array(iz, nhalos+iy, nhalos+ix)
!write(*, *) icol, iz, iy, ix, variable_array(icol, iz, 1)
               end do
            end do
         end do

         if (iblock == nblocks) then
            filter_input_file%ncstatus = nf90_put_var(filter_input_file%ncid, &
               dart_varids(varid), variable_array)
         end if
      end if
 
      deallocate(block_array)
   end do
end do

call nc_close_file(filter_input_file%ncid)

deallocate(spatial_array, variable_array)
deallocate(dart_varids, time_array)
deallocate(col_index)
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
real(r8), allocatable :: block_array(:, :, :)
real(r8), allocatable :: spatial_array(:)
real(r8), allocatable :: variable_array(:, :, :)

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

function assign_grid_files_array(nblocks) result(grid_files)
   
integer,          intent(in) :: nblocks
type(file_type), allocatable :: grid_files(:)

character(len=4) :: block_name
integer          :: iblock

allocate(grid_files(nblocks))

! Want more control over pathnames from namelist 
do iblock = 1, nblocks
   block_name = zero_fill(integer_to_string(iblock-1), 4)
   grid_files(iblock)%file_path = "../grid_files/grid_g" // block_name // ".nc"
   !!!grid_files(iblock)%file_path = "../NEW_RESTART_FILES/restartOut//grid_g" // block_name // ".nc"
end do

end function assign_grid_files_array

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
