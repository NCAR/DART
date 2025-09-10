module transform_state_mod

use netcdf
use types_mod,                  only : r4, r8, varnamelength, RAD2DEG
use netcdf_utilities_mod,       only : nc_open_file_readonly, nc_open_file_readwrite, &
                                       nc_close_file, nc_create_file, nc_end_define_mode
use utilities_mod,              only : find_namelist_in_file, check_namelist_read,  &
                                       error_handler, E_ERR, string_to_integer

use cube_sphere_grid_tools_mod, only : lat_lon_to_col_index, get_grid_delta

implicit none
private

public :: initialize_transform_state_mod,  model_to_dart, dart_to_model, &
          get_ensemble_range_from_command_line

integer  :: iunit, io

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'aether_cube_sphere/transform_state_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

type :: file_type
   character(len=256) :: file_path
   integer            :: ncid, unlimitedDimId, nDimensions, nVariables, nAttributes, formatNum
end type file_type

! It would be nice to get this information from the Aether input files, not possible for now
integer            :: np, nblocks, nhalos
namelist /transform_state_nml/ np, nblocks, nhalos

! Temporary switch between scalar and horizontal f10.7
logical  :: scalar_f10_7 = .false.

contains

!---------------------------------------------------------------

subroutine initialize_transform_state_mod()

call find_namelist_in_file('input.nml', 'transform_state_nml', iunit)
read(iunit, nml = transform_state_nml, iostat = io)
call check_namelist_read(iunit, io, 'transform_state_nml')

end subroutine initialize_transform_state_mod

!---------------------------------------------------------------

subroutine model_to_dart(aether_block_file_dir, dart_file_dir, ensemble_number)

character(len=*), intent(in) :: aether_block_file_dir, dart_file_dir
integer,          intent(in) :: ensemble_number

integer  :: iblock, dimid, length, ncols, dart_dimid(3), varid, xtype, nDimensions, nAtts
integer  :: ix, iy, iz, icol, ncstatus
integer  :: ntimes(nblocks), nxs(nblocks), nys(nblocks)
integer  :: ions_ntimes(nblocks), ions_nxs(nblocks), ions_nys(nblocks)
integer  :: neutrals_ntimes(nblocks), neutrals_nxs(nblocks), neutrals_nys(nblocks)
integer  :: haloed_nxs(nblocks), haloed_nys(nblocks)
integer  :: final_nzs, ions_final_nzs, neutrals_final_nzs
integer  :: filter_time_id, filter_alt_id, filter_lat_id, filter_lon_id 
integer  :: electron_varid, f10_7_varid
integer  :: grid_alt_id,   grid_lat_id,   grid_lon_id
integer  :: dimids(NF90_MAX_VAR_DIMS)
real(r8) :: blat, blon, del, half_del, f10_7_val
real(r8) :: time_array(1)
logical  :: add_to_electrons
character(len = 4) :: ensemble_string
character(len=NF90_MAX_NAME) :: name, attribute
integer,         allocatable :: col_index(:, :, :), filter_ions_ids(:), filter_neutrals_ids(:)
! The time variable in the block files is a double
! File for reading in variables from block file; These can be R4
real(r4),        allocatable :: spatial_array(:), variable_array(:, :, :), electron_array(:, :)
real(r4),        allocatable :: block_array(:, :, :), block_lats(:, :, :), block_lons(:, :, :)
type(file_type), allocatable :: ions_files(:), neutrals_files(:), grid_files(:)
type(file_type)              :: filter_file

! Open the grid and ions and neutrals files here for now with fixed directory names
ions_files = assign_block_file_names(nblocks, aether_block_file_dir, &
   'ions', ensemble_number)

neutrals_files = assign_block_file_names(nblocks, aether_block_file_dir, &
   'neutrals', ensemble_number)

grid_files = assign_block_file_names(nblocks, aether_block_file_dir, 'grid')

! Get grid spacing from number of points across each face
call get_grid_delta(np, del, half_del)

!======================== Get info on x, y and z dimensions from grid files

do iblock = 1, nblocks
   ! Open the grid files, read only
   grid_files(iblock)%ncid = nc_open_file_readonly(grid_files(iblock)%file_path)
end do

call get_aether_block_dimensions(grid_files, nblocks, nhalos, nxs, nys, final_nzs)
! Get the full dimension size with the halos for all blocks
haloed_nxs = nxs + 2*nhalos
haloed_nys = nys + 2*nhalos

! Compute the number of columns (without haloes)
ncols = sum(nxs(1:nblocks) * nys(1:nblocks))
write(*, *) 'ncols is ', ncols

!=============================== Check that ions files are consistent with grids =========

do iblock = 1, nblocks
   ! Open the ions block files, read only, and get the metadata
   ions_files(iblock)%ncid = nc_open_file_readonly(ions_files(iblock)%file_path)
end do

call get_aether_block_dimensions(ions_files, nblocks, nhalos, ions_nxs, ions_nys, ions_final_nzs)

! Check for inconsistent number of vertical levels in ion and grid files
if(ions_final_nzs .ne. final_nzs) &
   call error_handler(E_ERR, 'model_to_dart', &
      'Number of altitudes in grid and ions files differs', source, revision, revdate)

! Make sure ions and grid files have same horizontal sizes
if(any(ions_nxs .ne. nxs) .or. any(ions_nys .ne. nys)) &
   call error_handler(E_ERR, 'model_to_dart', &
      'Number of Latitudes and Longitudes in grid and ions files differ', source, revision, revdate)

!=============================== Check that neutrals files are consistent with grids =========

do iblock = 1, nblocks
   ! Open the neutrals block files, read only, and get the metadata
   neutrals_files(iblock)%ncid = nc_open_file_readonly(neutrals_files(iblock)%file_path)
end do

call get_aether_block_dimensions(neutrals_files, nblocks, nhalos, neutrals_nxs, neutrals_nys, &
   neutrals_final_nzs)

! Check for inconsistent number of vertical levels in ion and grid files
if(neutrals_final_nzs .ne. final_nzs) &
   call error_handler(E_ERR, 'model_to_dart', &
      'Number of altitudes in grid and neutrals files differs', source, revision, revdate)

! Make sure neutrals and grid files have same horizontal sizes
if(any(neutrals_nxs .ne. nxs) .or. any(neutrals_nys .ne. nys)) &
   call error_handler(E_ERR, 'model_to_dart', &
      'Number of Latitudes and Longitudes in grid and neutrals files differ', source, revision, revdate)

!==================================== Write dimensions in the filter nc file ========

! Initialize the filter file that will be created
ensemble_string = zero_fill(integer_to_string(ensemble_number), 4)
filter_file%file_path = trim(dart_file_dir) // 'filter_input_' // &
   ensemble_string // '.nc'

! Create the filter netcdf file
filter_file%ncid = nc_create_file(filter_file%file_path)

! Create dimensions in filter_file; save for use during variable definition
ncstatus = nf90_def_dim(filter_file%ncid, 'time', NF90_UNLIMITED, dart_dimid(3))
ncstatus = nf90_def_dim(filter_file%ncid, 'z',    final_nzs,      dart_dimid(2))
ncstatus = nf90_def_dim(filter_file%ncid, 'col',  ncols,          dart_dimid(1))

!=========================================================
! Create the variables from the grid files; lat, lon, alt

do varid = 1, 4
   ncstatus = nf90_inquire_variable(grid_files(1)%ncid, varid, name, xtype, nDimensions, dimids, nAtts)
   if (trim(name) == 'time') then
      ncstatus = nf90_def_var(filter_file%ncid, name, xtype, dart_dimid(3), filter_time_id)
   else if (trim(name) == 'Altitude') then
      ! Rename the 'z' variable as 'alt' so there isn't a dimension and a variable with the same name
      ncstatus = nf90_def_var(filter_file%ncid, 'alt', xtype, dart_dimid(2), filter_alt_id)
      ncstatus = nf90_put_att(filter_file%ncid, filter_alt_id, 'units', 'm')
      ncstatus = nf90_put_att(filter_file%ncid, filter_alt_id, 'long_name', &
         'height above mean sea level')
      grid_alt_id = varid
   else if (trim(name) == 'Latitude') then
      ncstatus = nf90_def_var(filter_file%ncid, 'lat', xtype, dart_dimid(1), filter_lat_id)
      ncstatus = nf90_put_att(filter_file%ncid, filter_lat_id, 'units', 'degrees_north')
      ncstatus = nf90_put_att(filter_file%ncid, filter_lat_id, 'long_name', 'latitude')
      grid_lat_id = varid
   else if (trim(name) == 'Longitude') then
      ncstatus = nf90_def_var(filter_file%ncid, 'lon', xtype, dart_dimid(1), filter_lon_id)
      ncstatus = nf90_put_att(filter_file%ncid, filter_lon_id, 'units', 'degrees_east')
      ncstatus = nf90_put_att(filter_file%ncid, filter_lon_id, 'long_name', 'longitude')
      grid_lon_id = varid
   else
      call error_handler(E_ERR, 'model_to_dart', &
         'Unexpected variable name in grid file ' // trim(name), source, revision, revdate)
   end if
end do

!=========================================================

! Allocate storage 

! Pointers to the different data fields in the filter nc file
allocate(filter_ions_ids(ions_files(1)%nVariables))
! Pointers to the different data fields in the filter nc file
allocate(filter_neutrals_ids(neutrals_files(1)%nVariables))
! Allocate ncols size temporary storage
allocate(spatial_array(ncols), variable_array(ncols, final_nzs, 1), electron_array(ncols, final_nzs))
! The col_index array will keep track of mapping from x and y for each block to final columns
allocate(col_index(nblocks, maxval(nys), maxval(nxs)))

! Allocate storage for the latitude and longitude from the blocks
allocate(block_lats (final_nzs, maxval(haloed_nys), maxval(haloed_nxs)), &
         block_lons (final_nzs, maxval(haloed_nys), maxval(haloed_nxs)), &
         block_array(final_nzs, maxval(haloed_nys), maxval(haloed_nxs)))

!=========================================================

! Get the metadata for variable fields from the ions block files
! The ions and neutrals files have time and all their physical variables, but not latitude, longitude, or altitude

! Illegal value of filter file index for default
filter_ions_ids = -99

! The filter_file is still in define mode. Create all of the variables before entering data mode.
do varid = 1, ions_files(1)%nVariables
   ncstatus = nf90_inquire_variable(ions_files(1)%ncid, varid, name, xtype, nDimensions, dimids, nAtts)
   ! Only time should actually occur in ions and neutrals files
   if (trim(name) /= 'time' .and. trim(name) /= 'z' .and. trim(name) /= 'lat' .and. trim(name) /= 'lon') then
      ncstatus = nf90_def_var(filter_file%ncid, name, xtype, dart_dimid, filter_ions_ids(varid))

      ! Add the units, same in all files so just get from the first
      ncstatus = nf90_get_att(ions_files(1)%ncid, varid, 'units', attribute)
      ncstatus = nf90_put_att(filter_file%ncid, filter_ions_ids(varid), 'units', attribute)
   end if
end do

!=========================================================

! Get the metadata for variable fields from the neutrals block files
! The ions and neutrals files have time and all their physical variables, but not latitude, longitude, or altitude

! Illegal value of filter file index for default
filter_neutrals_ids = -99

! The filter_file is still in define mode. Create all of the variables before entering data mode.
do varid = 1, neutrals_files(1)%nVariables
   ncstatus = nf90_inquire_variable(neutrals_files(1)%ncid, varid, name, xtype, nDimensions, dimids, nAtts)
   ! Only time should occur once we switch to ions and neutrals files
   if (trim(name) /= 'time' .and. trim(name) /= 'z' .and. trim(name) /= 'lat' .and. trim(name) /= 'lon') then
      ncstatus = nf90_def_var(filter_file%ncid, name, xtype, dart_dimid, filter_neutrals_ids(varid))

      ! Add the units, same in all files so just get from the first
      ncstatus = nf90_get_att(neutrals_files(1)%ncid, varid, 'units', attribute)
      ncstatus = nf90_put_att(filter_file%ncid, filter_neutrals_ids(varid), 'units', attribute)
   end if
end do

!=========================================================

! Add a derived vertical total electron content field
! xtype is currently set to value of last field from neutrals file
ncstatus = nf90_def_var(filter_file%ncid, 'ION_E', xtype, dart_dimid, electron_varid)
ncstatus = nf90_put_att(filter_file%ncid, electron_varid, 'units', '/m3')

if(scalar_f10_7) then
   ! Add a scalar F10.7 
   ncstatus = nf90_def_var(filter_file%ncid, 'F10.7', xtype, dart_dimid(3), f10_7_varid)
   ncstatus = nf90_put_att(filter_file%ncid, f10_7_varid, 'units', 'sfu: W/m^2/Hz')
   ncstatus = nf90_put_att(filter_file%ncid, f10_7_varid, 'long_name', 'Solar Radio Flux at 10.7 cm')
else
   ! Add a two-dimensional F10.7
   ! WARNING: QUANTITY AS PS UNTIL FURTHER STUDY
   ncstatus = nf90_def_var(filter_file%ncid, 'F10.7', xtype, &
      dart_dimid(1:3:2), f10_7_varid)
   ncstatus = nf90_put_att(filter_file%ncid, f10_7_varid, 'units', 'sfu: W/m^2/Hz')
   ncstatus = nf90_put_att(filter_file%ncid, f10_7_varid, 'long_name', 'Solar Radio Flux at 10.7 cm')
endif

! End of define mode for filter nc file, ready to add data
call nc_end_define_mode(filter_file%ncid)

!=========== Block to get lat lon and alt data from grid files =================

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
ncstatus = nf90_put_var(filter_file%ncid, filter_alt_id, block_array(:,1,1))

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
ncstatus = nf90_put_var(filter_file%ncid, filter_lat_id, spatial_array)

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
ncstatus = nf90_put_var(filter_file%ncid, filter_lon_id, spatial_array)

!=========== Copy data from ions files =================

! Electron density is sum of all the ion densities; sum them up
electron_array = 0.0_r8

! Will get full spatial field for one variable at a time
do varid = 1, ions_files(1)%nVariables
   ! Get metadata for this variable from first block file 
   ncstatus = nf90_inquire_variable(ions_files(1)%ncid, &
      varid, name, xtype, nDimensions, dimids, nAtts)

   ! See if this is a density; if so, needs to be added into electrons
   ncstatus = nf90_get_att(ions_files(1)%ncid, varid, 'units', attribute)
   add_to_electrons = trim(attribute) == '/m3'

   if(trim(name) == 'time') then
      ! Time must be the same in all files, so just deal with it from the first one
      ncstatus = nf90_get_var(ions_files(1)%ncid, varid, time_array)
      ncstatus = nf90_put_var(filter_file%ncid, filter_time_id, time_array)
   else
      ! Loop through all the blocks for this variable
      do iblock = 1, nblocks
         ! Read into the full 3Dblock array
         ncstatus = nf90_get_var(ions_files(iblock)%ncid, varid, block_array)
         
         do iy = 1, nys(iblock)
            do ix = 1, nxs(iblock)
               icol = col_index(iblock, iy, ix)
               do iz = 1, final_nzs
                  variable_array(icol, iz, 1) = block_array(iz, nhalos+iy, nhalos+ix)
                  ! Add into electrons if it is a density
                  if(add_to_electrons) electron_array(icol, iz) = &
                     electron_array(icol, iz) + variable_array(icol, iz, 1)
               end do
            end do
         end do

      end do
      ncstatus = nf90_put_var(filter_file%ncid, filter_ions_ids(varid), variable_array)
   end if

end do

!=========== Copy data from neutrals files =================

! Will get full spatial field for one variable at a time
do varid = 1, neutrals_files(1)%nVariables
   ! Get metadata for this variable from first block file 
   ncstatus = nf90_inquire_variable(neutrals_files(1)%ncid, &
      varid, name, xtype, nDimensions, dimids, nAtts)

   ! Already got time from ions files
   if(trim(name) .ne. 'time') then
      ! Loop through all the blocks for this variable
      do iblock = 1, nblocks
         ! Read into the full 3Dblock array
         ncstatus = nf90_get_var(neutrals_files(iblock)%ncid, varid, block_array)
         
         do iy = 1, nys(iblock)
            do ix = 1, nxs(iblock)
               icol = col_index(iblock, iy, ix)
               do iz = 1, final_nzs
                  variable_array(icol, iz, 1) = block_array(iz, nhalos+iy, nhalos+ix)
               end do
            end do
         end do

      end do
      ncstatus = nf90_put_var(filter_file%ncid, filter_neutrals_ids(varid), variable_array)
   end if

end do

!===================== Add in the additional variables ==================
! Write out the electron density field
variable_array(:, :, 1) = electron_array
ncstatus = nf90_put_var(filter_file%ncid, electron_varid, variable_array)

if(scalar_f10_7) then
   ! Add in f10.7 as a zero_dimensional field
   ncstatus = nf90_put_var(filter_file%ncid, f10_7_varid, 1.0_r8 * ensemble_number)
else
   ! Add in f10.7 as a two-dimensional field
   f10_7_val = 1.0_r8 * ensemble_number
   variable_array(:, 1, 1) = f10_7_val
   ncstatus = nf90_put_var(filter_file%ncid, f10_7_varid, variable_array(:, 1, 1))
endif

!===============================================================

! Close files and release storage
call nc_close_file(filter_file%ncid)
do iblock = 1, nblocks
   ! Close the grid files and ions and neutrals files
   call nc_close_file(grid_files(iblock)%ncid)
   call nc_close_file(ions_files(iblock)%ncid)
   call nc_close_file(neutrals_files(iblock)%ncid)
end do

deallocate(block_lats, block_lons, block_array, spatial_array, variable_array)
deallocate(electron_array, col_index, filter_ions_ids, filter_neutrals_ids)

end subroutine model_to_dart

!---------------------------------------------------------------

subroutine dart_to_model(dart_file_dir, aether_block_file_dir, ensemble_number)

character(len=*), intent(in) :: dart_file_dir, aether_block_file_dir
integer,          intent(in) :: ensemble_number

real(r8) :: blat, blon, del, half_del, f10_7_scalar
integer  :: iblock, dimid, length, ncols, varid, xtype, nDimensions, nAtts
integer  :: ix, iy, iz, icol, ncstatus, filter_varid
integer  :: nxs(nblocks), nys(nblocks), final_nzs
integer  :: ions_nxs(nblocks), ions_nys(nblocks), ions_final_nzs
integer  :: neutrals_nxs(nblocks), neutrals_nys(nblocks), neutrals_final_nzs
integer  :: haloed_nxs(nblocks), haloed_nys(nblocks)
integer  :: grid_alt_id,   grid_lat_id,   grid_lon_id
integer  :: dimids(NF90_MAX_VAR_DIMS)
character(len = 4) :: ensemble_string
character(len=NF90_MAX_NAME) :: name, attribute
integer,         allocatable :: col_index(:, :, :)
! File for reading in variables from block file; These can be R4
real(r4),        allocatable :: variable_array(:, :, :)
real(r4),        allocatable :: block_array(:, :, :), block_lats(:, :, :), block_lons(:, :, :)
type(file_type), allocatable :: ions_files(:), neutrals_files(:), grid_files(:)
type(file_type)              :: filter_file

! Open the grid and ions and neutrals files here for now with fixed directory names
ions_files = assign_block_file_names(nblocks, aether_block_file_dir, &
   'ions', ensemble_number)

neutrals_files = assign_block_file_names(nblocks, aether_block_file_dir, &
   'neutrals', ensemble_number)

grid_files = assign_block_file_names(nblocks, aether_block_file_dir, 'grid')

! Get grid spacing from number of points across each face
call get_grid_delta(np, del, half_del)

!======================== Get info on x, y and z dimensions from grid files
do iblock = 1, nblocks
   ! Open the grid files, read only
   grid_files(iblock)%ncid = nc_open_file_readonly(grid_files(iblock)%file_path)
end do

call get_aether_block_dimensions(grid_files, nblocks, nhalos, nxs, nys, final_nzs)
! Get the full dimension size with the halos for all blocks
haloed_nxs = nxs + 2*nhalos
haloed_nys = nys + 2*nhalos

! Compute the number of columns (without haloes)
ncols = sum(nxs(1:nblocks) * nys(1:nblocks))
write(*, *) 'ncols is ', ncols

!==================================== Open and check dimensions for ions files

! The ions block files need to be open read write
do iblock = 1, nblocks
   ions_files(iblock)%ncid = nc_open_file_readwrite(ions_files(iblock)%file_path)
end do

! Check that ions files are consistent with grids
call get_aether_block_dimensions(ions_files, nblocks, nhalos, ions_nxs, ions_nys, ions_final_nzs)

! Check for inconsistent number of vertical levels in ion and grid files
if(ions_final_nzs .ne. final_nzs) &
   call error_handler(E_ERR, 'model_to_dart', &
      'Number of altitudes in grid and ions files differs', source, revision, revdate)

! Make sure ions and grid files have same horizontal sizes
if(any(ions_nxs .ne. nxs) .or. any(ions_nys .ne. nys)) &
   call error_handler(E_ERR, 'model_to_dart', &
      'Number of Latitudes and Longitudes in grid and ion files differ', source, revision, revdate)

!==================================== Open and check dimensions for neutrals files

! The neutrals block files need to be open read write
do iblock = 1, nblocks
   neutrals_files(iblock)%ncid = nc_open_file_readwrite(neutrals_files(iblock)%file_path)
end do

! Check that neutrals files are consistent with grids
call get_aether_block_dimensions(neutrals_files, nblocks, nhalos, neutrals_nxs, neutrals_nys, neutrals_final_nzs)

! Check for inconsistent number of vertical levels in ion and grid files
if(neutrals_final_nzs .ne. final_nzs) &
   call error_handler(E_ERR, 'model_to_dart', &
      'Number of altitudes in grid and neutrals files differs', source, revision, revdate)

! Make sure neutrals and grid files have same horizontal sizes
if(any(neutrals_nxs .ne. nxs) .or. any(neutrals_nys .ne. nys)) &
   call error_handler(E_ERR, 'model_to_dart', &
      'Number of Latitudes and Longitudes in grid and neutrals files differ', source, revision, revdate)

!=========== Get lat lon and alt from grid files =================

! Find the latitude and longitude information from the grid files and get the column mapping
ncstatus = nf90_inq_varid(grid_files(1)%ncid, 'Altitude',  grid_alt_id)
ncstatus = nf90_inq_varid(grid_files(1)%ncid, 'Latitude',  grid_lat_id)
ncstatus = nf90_inq_varid(grid_files(1)%ncid, 'Longitude', grid_lon_id)

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

! Open the filter file that will be read
ensemble_string = zero_fill(integer_to_string(ensemble_number), 4)
filter_file%file_path = trim(dart_file_dir) // 'filter_output_' // &
   ensemble_string // '.nc'
! Open the filter netcdf file
filter_file%ncid = nc_open_file_readonly(filter_file%file_path)

!==========================================================================
! Loop through ions fields and replace with values from filter_file

ncstatus = nf90_inquire(ions_files(1)%ncid, ions_files(1)%nDimensions, &
   ions_files(1)%nVariables, ions_files(1)%nAttributes,  ions_files(1)%unlimitedDimId, &
   ions_files(1)%formatNum)

! Get full spatial field for one variable at a time
do varid = 1, ions_files(1)%nVariables
   ! Get metadata for this variable from first block file 
   ncstatus = nf90_inquire_variable(ions_files(1)%ncid, &
      varid, name, xtype, nDimensions, dimids, nAtts)
   if(trim(name) .ne. 'time' .and. trim(name) .ne. 'Altitude' .and. trim(name) .ne. 'Latitude' &
      .and. trim(name) .ne. 'Longitude') then
      ! See if this variable is also in the filter output file
      ncstatus = nf90_inq_varid(filter_file%ncid, trim(name), filter_varid)
      ! Check on failed ncstatus. 0 is successful but should use the proper name
      if(ncstatus == 0) then
         ! Read this field from filter file 
         ncstatus = nf90_get_var(filter_file%ncid, filter_varid, variable_array)

! CAN WE UPDATE THE HALOS TOO WHEN WRITING BACK???
         ! Loop through all the blocks for this variable
         block_array = 0.0_r8
         do iblock = 1, nblocks
            do iy = 1, nys(iblock)
               do ix = 1, nxs(iblock)
                  icol = col_index(iblock, iy, ix)
                  do iz = 1, final_nzs
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

!==========================================================================
! Loop through neutrals fields and replace with values from filter

ncstatus = nf90_inquire(neutrals_files(1)%ncid, neutrals_files(1)%nDimensions, &
   neutrals_files(1)%nVariables, neutrals_files(1)%nAttributes,  neutrals_files(1)%unlimitedDimId, &
   neutrals_files(1)%formatNum)

! Get full spatial field for one variable at a time
do varid = 1, neutrals_files(1)%nVariables
   ! Get metadata for this variable from first block file 
   ncstatus = nf90_inquire_variable(neutrals_files(1)%ncid, &
      varid, name, xtype, nDimensions, dimids, nAtts)
   if(trim(name) .ne. 'time' .and. trim(name) .ne. 'Altitude' .and. trim(name) .ne. 'Latitude' &
      .and. trim(name) .ne. 'Longitude') then
      ! See if this variable is also in the filter output file
      ncstatus = nf90_inq_varid(filter_file%ncid, trim(name), filter_varid)
      ! Check on failed ncstatus. 0 is successful but should use the proper name
      if(ncstatus == 0) then
         ! Read this field from filter file 
         ncstatus = nf90_get_var(filter_file%ncid, filter_varid, variable_array)

         ! Loop through all the blocks for this variable
         block_array = 0.0_r8
         do iblock = 1, nblocks
            do iy = 1, nys(iblock)
               do ix = 1, nxs(iblock)
                  icol = col_index(iblock, iy, ix)
                  do iz = 1, final_nzs
                     block_array(iz, nhalos+iy, nhalos+ix) = variable_array(icol, iz, 1)
                  end do
               end do
            end do
            ! Write into the full file for this block
            ncstatus = nf90_put_var(neutrals_files(iblock)%ncid, varid, block_array)
         end do
      endif
   end if
end do

!==============================================================================

! Need more information about where F10.7 will be in Aether input files to complete copy back
ncstatus = nf90_inq_varid(filter_file%ncid, 'F10.7', filter_varid)
if(scalar_f10_7) then
   ! Read a scalar f10_7 value
   ncstatus = nf90_get_var(filter_file%ncid, filter_varid, f10_7_scalar)
   write(*, *) 'reading a scalar f10.7 from filter file', f10_7_scalar
else
   ! Read a column sized f10_7 value
   ncstatus = nf90_get_var(filter_file%ncid, filter_varid, variable_array(:, 1, 1))
   ! Average the updated value over all the columns
   f10_7_scalar = sum(variable_array(:, 1, 1)) / ncols
   write(*, *) 'reading column f10.7 ', f10_7_scalar
endif
! Write the updated F10.7 to the appropriate Aether file
!!! NEED MORE INFO TO IMPLEMENT

!==============================================================================

! Free storage and close files
deallocate(col_index, variable_array, block_lats, block_lons, block_array)

! Close the netcdf files
call nc_close_file(filter_file%ncid)
do iblock = 1, nblocks
   ! Close the grid files and ions and neutrals files
   call nc_close_file(grid_files(iblock)%ncid)
   call nc_close_file(ions_files(iblock)%ncid)
   call nc_close_file(neutrals_files(iblock)%ncid)
end do

end subroutine dart_to_model

!--------------------------------------------------------------------

subroutine get_aether_block_dimensions(files, nblocks, nhalos, nxs, nys, nzs) 

integer,         intent(in)    :: nblocks, nhalos
type(file_type), intent(inout) :: files(nblocks)
integer,         intent(out)   :: nxs(nblocks), nys(nblocks), nzs

integer :: iblock, b_nzs(nblocks), ncstatus, dimid, length
character(len=NF90_MAX_NAME) :: name

do iblock = 1, nblocks
   ! Get info about the block file
   ncstatus = nf90_inquire(files(iblock)%ncid, files(iblock)%nDimensions, &
      files(iblock)%nVariables, files(iblock)%nAttributes, &
       files(iblock)%unlimitedDimId, files(iblock)%formatNum)

   ! Verify that a single time level exists
   ncstatus = nf90_inq_dimid(files(iblock)%ncid, 'time', dimid)
   ncstatus = nf90_inquire_dimension(files(iblock)%ncid, dimid, name, length)
   if(length .ne. 1 .or. ncstatus .ne. 0) &
      call error_handler(E_ERR, 'get_aether_block_dimensions', &
         'Number of times in input block files should be 1', source, revision, revdate)

   ! Get the length of x dimension
   ncstatus = nf90_inq_dimid(files(iblock)%ncid, 'x', dimid)
   ncstatus = nf90_inquire_dimension(files(iblock)%ncid, dimid, name, length)
   if(ncstatus .ne. 0) &
      call error_handler(E_ERR, 'get_aether_block_dimensions', &
         'input block files must have x dimension', source, revision, revdate)
   nxs(iblock)         = length-2*nhalos

   ! Get the length of y dimension
   ncstatus = nf90_inq_dimid(files(iblock)%ncid, 'y', dimid)
   ncstatus = nf90_inquire_dimension(files(iblock)%ncid, dimid, name, length)
   if(ncstatus .ne. 0) &
      call error_handler(E_ERR, 'get_aether_block_dimensions', &
         'input block files must have y dimension', source, revision, revdate)
   nys(iblock)         = length-2*nhalos

   ! Get the length of z dimension
   ncstatus = nf90_inq_dimid(files(iblock)%ncid, 'z', dimid)
   ncstatus = nf90_inquire_dimension(files(iblock)%ncid, dimid, name, length)
   if(ncstatus .ne. 0) &
      call error_handler(E_ERR, 'get_aether_block_dimensions', &
         'input block files must have z dimension', source, revision, revdate)
   b_nzs(iblock) = length

end do

! Make sure all blocks have same number of vertical levels
if(any(b_nzs - b_nzs(1) .ne. 0)) then
   call error_handler(E_ERR, 'model_to_dart', &
         'block files have different lengths for z dimension', source, revision, revdate)
else
   nzs = b_nzs(1)
endif

! Make sure all block files have the same number of attributes
if(any(files(:)%nAttributes .ne. files(1)%nAttributes)) &
   call error_handler(E_ERR, 'model_to_dart', &
         'All blocks must have same nunber of variables', source, revision, revdate)

end subroutine get_aether_block_dimensions

!---------------------------------------------------------------

subroutine get_ensemble_range_from_command_line(start_ensemble, end_ensemble)

integer, intent(out) :: start_ensemble, end_ensemble

! Gets the first and last ensemble members to be converted from command line

character(len=4) :: start_ensemble_string, end_ensemble_string
integer          :: nargs

nargs = command_argument_count()

if (nargs /= 2) &
   call error_handler(E_ERR, 'get_ensemble_range_from_command_line', &
      'starting and ending ensemble members must be in command line argument')

call get_command_argument(1, start_ensemble_string)
call get_command_argument(2, end_ensemble_string)

! Convert these to integer values
read(start_ensemble_string, *) start_ensemble
read(end_ensemble_string,   *) end_ensemble

! Not prepared to deal with more than 4 digit ensemble count
if(start_ensemble > 9999 .or. end_ensemble > 9999) &
   call error_handler(E_ERR, 'get_ensemble_range_from_command_line', &
      'Ensemble numbers on command line must be less than 10000')

end subroutine get_ensemble_range_from_command_line

!---------------------------------------------------------------

function assign_block_file_names(nblocks, directory, &
   file_prefix, ensemble_number) result(block_files)

integer,          intent(in)             :: nblocks
character(len=*), intent(in)             :: directory
character(len=*), intent(in)             :: file_prefix
integer,          intent(in), optional   :: ensemble_number
type(file_type),  allocatable            :: block_files(:)

character(len=256) :: file
character(len=4)   :: block_num, ensemble_string
integer            :: iblock

! Storage for each block
allocate(block_files(nblocks))

if(present(ensemble_number)) then
   ensemble_string = zero_fill(integer_to_string(ensemble_number - 1), 4)
endif

do iblock = 1, nblocks
   block_num = zero_fill(integer_to_string(iblock - 1), 4)
   file = trim(directory) // trim(file_prefix) 
   ! Add in ensemble member if needed
   if(present(ensemble_number)) then
      file = trim(file) // '_m' // ensemble_string
   endif
   block_files(iblock)%file_path = trim(file) //  '_g' // block_num // '.nc'
end do


end function assign_block_file_names

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
   call error_handler(E_ERR, 'zero_fill', &
      'Input string is longer than desired output => ensemble size too large', &
      source, revision, revdate)
else if (difference_of_string_lengths > 0) then
   do string_index = 1, difference_of_string_lengths
      filled_string(string_index:string_index) = '0'
   end do
end if

filled_string(difference_of_string_lengths+1:desired_length) = trim(string)

end function zero_fill

!---------------------------------------------------------------

end module transform_state_mod
