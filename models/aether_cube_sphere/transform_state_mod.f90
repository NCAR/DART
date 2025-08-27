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

type(file_type), allocatable :: ions_files(:), grid_files(:)
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

ions_files = assign_block_files_array(nblocks, restart_ensemble_member, restart_directory, &
   restart_file_prefix, restart_file_middle, restart_file_suffix)

grid_files = assign_grid_files_array(nblocks)

end subroutine initialize_transform_state_mod

!---------------------------------------------------------------

subroutine model_to_dart()

integer :: iblock, dimid, length, ncols, dart_dimid(3), varid, xtype, nDimensions, nAtts
integer :: ix, iy, iz, icol, ncstatus
integer :: ntimes(nblocks), nzs(nblocks), nxs(nblocks), nys(nblocks)
integer :: ion_ntimes(nblocks), ion_nzs(nblocks), ion_nxs(nblocks), ion_nys(nblocks)
integer :: haloed_nxs(nblocks), haloed_nys(nblocks)
integer :: final_ntimes, final_nzs
integer, allocatable :: filter_ions_ids(:)
integer :: filter_time_id, filter_alt_id, filter_lat_id, filter_lon_id
integer :: grid_time_id,   grid_alt_id,   grid_lat_id,   grid_lon_id
integer :: dimids(NF90_MAX_VAR_DIMS)
real(r8) :: blat, blon, del, half_del
character(len=NF90_MAX_NAME) :: name, attribute
integer,  allocatable         :: col_index(:, :, :)
! The time variable in the block files is a double
real(r8), allocatable, dimension(:) :: time_array
! File for reading in variables from block file; These can probably be R4
real(r4), allocatable :: spatial_array(:), variable_array(:, :, :)
real(r4), allocatable :: block_array(:, :, :), block_lats(:, :, :), block_lons(:, :, :)

! Get grid spacing from number of points across each face
call get_grid_delta(np, del, half_del)

!==================================== get info from grid file block ===================

do iblock = 1, nblocks
   ! Open the grid files, read only
   grid_files(iblock)%ncid = nc_open_file_readonly(grid_files(iblock)%file_path)
   ! There doesn't seem to be a helper procedure corresponding to nf90_inquire in
   ! netcdf_utilities_mod so this uses the external function directly from the netcdf library
   ncstatus = nf90_inquire(grid_files(iblock)%ncid, &
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
      ncstatus = nf90_inquire_dimension(grid_files(iblock)%ncid, dimid, name, length)

      if (trim(name) == 'time') then
         ! Don't care about times in the grid files
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
final_nzs = nzs(1)

! Compute the number of columns (without haloes)
ncols = sum(nxs(1:nblocks) * nys(1:nblocks))
write(*, *) 'ncols is ', ncols


!==================================== end of get info from grid file block ===================

! Start with ion files, add in neutrals later, test with old aether_restart files first

do iblock = 1, nblocks
   ! Open the ions block files, read only, and get the metadata
   ions_files(iblock)%ncid = nc_open_file_readonly(ions_files(iblock)%file_path)
   ncstatus = nf90_inquire(ions_files(iblock)%ncid, &
      ions_files(iblock)%nDimensions, ions_files(iblock)%nVariables, &
      ions_files(iblock)%nAttributes,  ions_files(iblock)%unlimitedDimId, &
      ions_files(iblock)%formatNum)

   ! Loop through each of the dimensions to find the metadata values
   do dimid = 1, ions_files(iblock)%nDimensions
      ncstatus = nf90_inquire_dimension(ions_files(iblock)%ncid, dimid, name, length)

      if (trim(name) == 'time') then
         !JLA Error if more than one time???
         ion_ntimes(iblock) = length
      else if (trim(name) == 'x') then
         ion_nxs(iblock)         = length-2*nhalos
      else if (trim(name) == 'y') then
         ion_nys(iblock)        = length-2*nhalos
      else if (trim(name) == 'z') then
         ion_nzs(iblock) = length
      end if
   end do

end do

! Do some consistency checks to make sure the files have the same number of variables, levels and times
! Comparison is to grid files for times and levels, but just among block files for number of variables
if(any(ion_ntimes - ion_ntimes(1) .ne. 0)) then
   write(*, *) 'inconsistent ntimes in ion files'
   stop
endif

final_ntimes = ion_ntimes(1)

if(any(ion_nzs - final_nzs .ne. 0)) then
   write(*, *) 'inconsistent nunber of vertical levels in ion files'
   stop
endif
if(any(ions_files(:)%nVariables - ions_files(1)%nVariables .ne. 0)) then
   write(*, *) 'inconsistent number of variables in ion files'
   stop
endif

! Allocate ncols size temporary storage
allocate(spatial_array(ncols), variable_array(ncols, final_nzs, final_ntimes), time_array(final_ntimes))

! Initialize the filter input file that will be created
filter_input_file = assign_filter_file(dart_ensemble_member, filter_directory, &
   filter_input_prefix, '.nc')
! Create the filter netcdf file
filter_input_file%ncid = nc_create_file(filter_input_file%file_path)

! Create dimensions in filter_input_file; save for use during variable definition
ncstatus = nf90_def_dim(filter_input_file%ncid, 'time', NF90_UNLIMITED, dart_dimid(3))
ncstatus = nf90_def_dim(filter_input_file%ncid, 'z',    nzs(1),         dart_dimid(2))
ncstatus = nf90_def_dim(filter_input_file%ncid, 'col',  ncols,          dart_dimid(1))

!=========================================================
! Create the variables from the grid files first. Should have only time, lat, lon, alt

do varid = 1, 4
   ncstatus = nf90_inquire_variable(grid_files(1)%ncid, varid, name, xtype, nDimensions, dimids, nAtts)
   if (trim(name) == 'time') then
      ncstatus = nf90_def_var(filter_input_file%ncid, name, xtype, dart_dimid(3), filter_time_id)
      grid_time_id = varid
   else if (trim(name) == 'Altitude') then
      ! Rename the 'z' variable as 'alt' so there isn't a dimension and a variable with the same name
      ncstatus = nf90_def_var(filter_input_file%ncid, 'alt', xtype, dart_dimid(2), filter_alt_id)
      ncstatus = nf90_put_att(filter_input_file%ncid, filter_alt_id, 'units', 'm')
      ncstatus = nf90_put_att(filter_input_file%ncid, filter_alt_id, 'long_name', &
         'height above mean sea level')
      grid_alt_id = varid
   else if (trim(name) == 'Latitude') then
      ncstatus = nf90_def_var(filter_input_file%ncid, 'lat', xtype, dart_dimid(1), filter_lat_id)
      ncstatus = nf90_put_att(filter_input_file%ncid, filter_lat_id, 'units', 'degrees_north')
      ncstatus = nf90_put_att(filter_input_file%ncid, filter_lat_id, 'long_name', 'latitude')
      grid_lat_id = varid
   else if (trim(name) == 'Longitude') then
      ncstatus = nf90_def_var(filter_input_file%ncid, 'lon', xtype, dart_dimid(1), filter_lon_id)
      ncstatus = nf90_put_att(filter_input_file%ncid, filter_lon_id, 'units', 'degrees_east')
      ncstatus = nf90_put_att(filter_input_file%ncid, filter_lon_id, 'long_name', 'longitude')
      grid_lon_id = varid
   else
      write(*, *) 'Unexpected variable name in grid file', trim(name)
      stop
   end if

end do

!=========================================================

! Now get the other data fields from the block files (soon to be ions and neutrals)
! The ions and neutrals files have time and all their physical variables, but not latitude, longitude, or altitude
! aether_restarts being used for tests also have lat, lon, alt

! Pointers to the different data fields in the filter nc file
allocate(filter_ions_ids(ions_files(1)%nVariables))

! The filter_input_file is still in define mode. Create all of the variables before entering data mode.
do varid = 1, ions_files(1)%nVariables
   ncstatus = nf90_inquire_variable(ions_files(1)%ncid, varid, name, xtype, nDimensions, dimids, nAtts)
   ! Only time should occur once we switch to ions and neutrals files
   if (trim(name) /= 'time' .and. trim(name) /= 'z' .and. trim(name) /= 'lat' .and. trim(name) /= 'lon') then
      ncstatus = nf90_def_var(filter_input_file%ncid, name, xtype, dart_dimid, filter_ions_ids(varid))
      ! Note that the filter_ions_id maps from the ids for fields in the ions files

      ! Add the units, same in all files so just get from the first
      ncstatus = nf90_get_att(ions_files(1)%ncid, varid, 'units', attribute)
      ncstatus = nf90_put_att(filter_input_file%ncid, filter_ions_ids(varid), 'units', attribute)
   end if
end do

! End of define mode for filter nc file, ready to add data
call nc_end_define_mode(filter_input_file%ncid)

!=========== Block to get lat lon and alt from grid files =================

! The col_index array will keep track of mapping from x and y for each block to final columns
allocate(col_index(nblocks, maxval(nys), maxval(nxs)))

! Allocate storage for the latitude and longitude from the blocks
allocate(block_lats(final_nzs, maxval(haloed_nys), maxval(haloed_nxs)), &
         block_lons(final_nzs, maxval(haloed_nys), maxval(haloed_nxs)), &
         block_array(final_nzs, maxval(haloed_nys), maxval(haloed_nxs)))

! Loop through all the blocks for this variable
do iblock = 1, nblocks
   ! Get the latitude and longitude full arrays 
   ncstatus = nf90_get_var(grid_files(iblock)%ncid, grid_lat_id, block_lats)
   ncstatus = nf90_get_var(grid_files(iblock)%ncid, grid_lon_id, block_lons)

   ! Compute the col_index for each of the horizontal locations in this block 
   do ix = 1, nxs(iblock)
      do iy = 1, nys(iblock)
         blat = block_lats(1, nhalos + iy, nhalos + ix);
         blon = block_lons(1, nhalos + iy, nhalos + ix);
         col_index(iblock, iy, ix) = lat_lon_to_col_index(blat, blon, del, half_del, np)
      end do
   end do
end do

! Only need altitude from 1 block
ncstatus = nf90_get_var(grid_files(1)%ncid, grid_alt_id, block_array)
ncstatus = nf90_put_var(filter_input_file%ncid, &
   filter_alt_id, block_array(:,1,1))

! Loop through blocks to get lat values
do iblock = 1, nblocks
   ! Get lat values for this block
   ncstatus = nf90_get_var(grid_files(iblock)%ncid, grid_lat_id, block_array)
   do iy = 1, nys(iblock)
      do ix = 1, nxs(iblock)
         icol = col_index(iblock, iy, ix)
         spatial_array(icol) = block_array(1, nhalos+iy, nhalos+ix) * RAD2DEG
      end do
   end do
end do
ncstatus = nf90_put_var(filter_input_file%ncid, filter_lat_id, spatial_array)

! Loop through blocks to get lon values
do iblock = 1, nblocks
   ! Get lon values for this block
   ncstatus = nf90_get_var(grid_files(iblock)%ncid, grid_lon_id, block_array)
   do iy = 1, nys(iblock)
      do ix = 1, nxs(iblock)
         icol = col_index(iblock, iy, ix)
         spatial_array(icol) = block_array(1, nhalos+iy, nhalos+ix) * RAD2DEG
      end do
   end do
end do
ncstatus = nf90_put_var(filter_input_file%ncid, filter_lon_id, spatial_array)


!=========== End of Block to get lat lon and alt from grid files =================

! Will get full spatial field for one variable at a time
do varid = 1, ions_files(1)%nVariables
   ! Get metadata for this variable from first block file 
   ncstatus = nf90_inquire_variable(ions_files(1)%ncid, &
      varid, name, xtype, nDimensions, dimids, nAtts)

   if(trim(name) == 'time') then
      ! Time must be the same in all files, so just deal with it from the first one
      ncstatus = nf90_get_var(ions_files(1)%ncid, varid, time_array)
      ncstatus = nf90_put_var(filter_input_file%ncid, filter_time_id, time_array)
   else
      ! Loop through all the blocks for this variable
      do iblock = 1, nblocks
         ! Read into the full 3Dblock array
         ncstatus = nf90_get_var(ions_files(iblock)%ncid, varid, block_array)
         
         do iy = 1, nys(iblock)
            do ix = 1, nxs(iblock)
               icol = col_index(iblock, iy, ix)
               do iz = 1, nzs(1)
                  variable_array(icol, iz, 1) = block_array(iz, nhalos+iy, nhalos+ix)
               end do
            end do
         end do

      end do
      ncstatus = nf90_put_var(filter_input_file%ncid, filter_ions_ids(varid), variable_array)
   end if

end do

call nc_close_file(filter_input_file%ncid)
do iblock = 1, nblocks
   ! Close the grid files and ions files
   call nc_close_file(grid_files(iblock)%ncid)
   call nc_close_file(ions_files(iblock)%ncid)
end do

deallocate(block_lats, block_lons, block_array, spatial_array, variable_array)
deallocate(time_array, col_index, filter_ions_ids)

end subroutine model_to_dart

!---------------------------------------------------------------

subroutine dart_to_model

real(r8) :: blat, blon, del, half_del
integer :: iblock, dimid, length, ncols, varid, xtype, nDimensions, nAtts
integer :: ix, iy, iz, icol, ncstatus, filter_varid
integer :: ntimes(nblocks), nzs(nblocks), nxs(nblocks), nys(nblocks)
integer :: final_nzs
integer :: haloed_nxs(nblocks), haloed_nys(nblocks)
integer :: filter_alt_id, filter_lat_id, filter_lon_id
integer :: grid_alt_id,   grid_lat_id,   grid_lon_id
integer :: dimids(NF90_MAX_VAR_DIMS)
character(len=NF90_MAX_NAME) :: name, attribute
integer,  allocatable         :: col_index(:, :, :)
! File for reading in variables from block file; These can probably be R4
real(r4), allocatable :: variable_array(:, :, :)
real(r4), allocatable :: block_array(:, :, :), block_lats(:, :, :), block_lons(:, :, :)

! JUST DEAL WITH ONE TIME LEVEL?

! Get grid spacing from number of points across each face
call get_grid_delta(np, del, half_del)

!==================================== get info from grid file block ===================

do iblock = 1, nblocks
   ! Open the grid files, read only
   grid_files(iblock)%ncid = nc_open_file_readonly(grid_files(iblock)%file_path)
   ! Get the info for this block
   ncstatus = nf90_inquire(grid_files(iblock)%ncid, &
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
      ncstatus = nf90_inquire_dimension(grid_files(iblock)%ncid, dimid, name, length)

      if (trim(name) == 'time') then
         ! Don't care about times in the grid files
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
final_nzs = nzs(1)

! Compute the number of columns (without haloes)
ncols = sum(nxs(1:nblocks) * nys(1:nblocks))
write(*, *) 'ncols is ', ncols

!==================================== end of get info from grid file block ===================

! Find the latitude and longitude information from the grid files and get the column mapping
ncstatus = nf90_inq_varid(grid_files(1)%ncid, 'Altitude',  grid_alt_id)
ncstatus = nf90_inq_varid(grid_files(1)%ncid, 'Latitude',  grid_lat_id)
ncstatus = nf90_inq_varid(grid_files(1)%ncid, 'Longitude', grid_lon_id)

!=========== Block to get lat lon and alt from grid files =================

! The col_index array will keep track of mapping from x and y for each block to final columns
! Allocate storage for the latitude and longitude from the blocks
allocate(col_index(nblocks, maxval(nys), maxval(nxs)), &
         variable_array(ncols, final_nzs, 1), &
         block_lats(final_nzs,  maxval(haloed_nys), maxval(haloed_nxs)), &
         block_lons(final_nzs,  maxval(haloed_nys), maxval(haloed_nxs)), &
         block_array(final_nzs, maxval(haloed_nys), maxval(haloed_nxs)))

! Loop through all the blocks for this variable
do iblock = 1, nblocks
   ! Get the latitude and longitude full arrays 
   ncstatus = nf90_get_var(grid_files(iblock)%ncid, grid_lat_id, block_lats)
   ncstatus = nf90_get_var(grid_files(iblock)%ncid, grid_lon_id, block_lons)

   ! Compute the col_index for each of the horizontal locations in this block 
   do ix = 1, nxs(iblock)
      do iy = 1, nys(iblock)
         blat = block_lats(1, nhalos + iy, nhalos + ix);
         blon = block_lons(1, nhalos + iy, nhalos + ix);
         col_index(iblock, iy, ix) = lat_lon_to_col_index(blat, blon, del, half_del, np)
      end do
   end do
end do

!=========== End of Block to get lat lon and alt from grid files =================

! The ions block files need to be open read write
do iblock = 1, nblocks
   ions_files(iblock)%ncid = nc_open_file_readwrite(ions_files(iblock)%file_path)
end do

! Get file name for filter_output_file that will be read
filter_output_file = assign_filter_file(dart_ensemble_member, filter_directory, &
   filter_output_prefix, '.nc')
filter_output_file%ncid = nc_open_file_readonly(filter_output_file%file_path)

!=========== Block to loop through ions fields and replace with valued from filter_output
ncstatus = nf90_inquire(ions_files(1)%ncid, ions_files(1)%nDimensions, &
   ions_files(1)%nVariables, ions_files(1)%nAttributes,  ions_files(1)%unlimitedDimId, &
   ions_files(1)%formatNum)

! Will get full spatial field for one variable at a time
do varid = 1, ions_files(1)%nVariables
   ! Get metadata for this variable from first block file 
   ncstatus = nf90_inquire_variable(ions_files(1)%ncid, &
      varid, name, xtype, nDimensions, dimids, nAtts)
   write(*, *) 'var loop ', varid, name
   if(trim(name) .ne. 'time' .and. trim(name) .ne. 'Altitude' .and. trim(name) .ne. 'Latitude' &
      .and. trim(name) .ne. 'Longitude') then
      ! See if this variable is also in the filter output file
      ncstatus = nf90_inq_varid(filter_output_file%ncid, trim(name), filter_varid)
      ! Check on failed ncstatus. 0 is successful but should use the proper name
      if(ncstatus == 0) then
         write(*, *) 'fetching variable ', trim(name)
         ! Read this field from filter_output file 
         ncstatus = nf90_get_var(filter_input_file%ncid, filter_varid, variable_array)

         ! Loop through all the blocks for this variable
         do iblock = 1, nblocks
            do iy = 1, nys(iblock)
               do ix = 1, nxs(iblock)
                  icol = col_index(iblock, iy, ix)
                  do iz = 1, nzs(1)
                     block_array(iz, nhalos+iy, nhalos+ix) = variable_array(icol, iz, 1)
                  end do
               end do
            end do
            ! Write into the full file for this block
            ncstatus = nf90_put_var(ions_files(iblock)%ncid, varid, block_array)
         end do
      endif
   end if
end do

deallocate(col_index, variable_array, block_lats, block_lons, block_array)

! Close the netcdf files
call nc_close_file(filter_output_file%ncid)
do iblock = 1, nblocks
   ! Close the grid files and ions files
   call nc_close_file(grid_files(iblock)%ncid)
   call nc_close_file(ions_files(iblock)%ncid)
end do

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
