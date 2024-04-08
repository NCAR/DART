program create_geometry_file

use               netcdf
use             sort_mod, only : index_sort
use            types_mod, only : i8, r8, PI, RAD2DEG, DEG2RAD
use         location_mod, only : location_type, get_dist, get_location, set_location, VERTISHEIGHT
use        utilities_mod, only : initialize_utilities, finalize_utilities, find_namelist_in_file, &
                                 check_namelist_read
use netcdf_utilities_mod, only : nc_open_file_readonly, nc_get_variable_size, nc_get_variable, &
                                 nc_create_file, nc_end_define_mode, nc_close_file
use  transform_state_mod, only : file_type, integer_to_string, zero_fill
use quad_utils_mod,       only : in_quad

implicit none

integer  :: iunit, io

! There will always be exactly 8 vertices on a cube sphere, regardless of how many blocks
integer, parameter               :: nvertex_columns = 8
integer, parameter               :: nside_faces = 4
integer, parameter               :: nvertex_neighbors = 3
integer, parameter               :: nquad_neighbors = 4
integer, dimension(nvertex_columns, 3) :: matrix_of_corner_indices
logical                          :: is_a_vertex
integer, dimension(3)            :: corner_indices

integer, dimension(3) :: variable_size

integer :: nzs_per_centers_block, nys_per_centers_block, nxs_per_centers_block, &
           nzs_per_corners_block, nys_per_corners_block, nxs_per_corners_block

integer :: truncated_nys_per_centers_block, truncated_nxs_per_centers_block, &
           truncated_nys_per_side_corners_block, truncated_nxs_per_side_corners_block, &
           truncated_nys_per_top_bottom_corners_block, truncated_nxs_per_top_bottom_corners_block

integer :: center_ncols, corner_ncols

real(r8), allocatable, dimension(:, :, :) :: center_altitude_block

real(r8), allocatable, dimension(:, :, :) :: center_lats_block, center_lons_block, & 
                                             corner_lats_block, corner_lons_block

real(r8), allocatable, dimension(:, :, :) :: center_lats_truncated, center_lons_truncated, &
                                             corner_lats_side_truncated, corner_lons_side_truncated, &
                                             corner_lats_top_bottom_truncated, corner_lons_top_bottom_truncated

real(r8), allocatable, dimension(:)       :: center_latitude, center_longitude, &
                                             quad_latitude, quad_longitude

integer, allocatable, dimension(:, :)     :: quad_neighbor_indices

integer, dimension(nvertex_columns, nvertex_neighbors) :: vertex_neighbor_indices

real(r8), dimension(nvertex_columns)            :: vertex_latitude, vertex_longitude

real(r8), dimension(3)                    :: location

! Declare variables for the geometry_file that will be output by this program
real(r8)                         :: min_alt
real(r8)                         :: max_lon
real(r8)                         :: min_lon
real(r8)                         :: max_lat
real(r8)                         :: min_lat

type :: vertex
   type(location_type) :: loc
   integer, dimension(3) :: column_index_of_neighbors
end type vertex

type :: quad
   type(location_type) :: loc
   logical  :: spans_prime_meridian
   logical  :: encloses_pole
   integer, dimension(4) :: column_index_of_neighbors
end type quad

! This is an intentionally unecessary nesting of a derived type just so that the lat and lon of each
! center_location is accessed using the same %loc attibute syntax
type :: center_location
   type(location_type) :: loc
end type center_location

type(vertex), dimension(nvertex_columns) :: vertexs
type(quad), allocatable, dimension(:) :: quads
type(center_location), allocatable, dimension(:) :: center_locations

integer                          :: ivertex
integer                          :: iquad

integer :: icol, ix, iy, nx, ny, inorth_pole_quad_column, isouth_pole_quad_column
real(r8), allocatable, dimension(:) :: distances_from_corners_to_centers

character(len=256) :: restart_directory, grid_directory, filter_directory
namelist /directory_nml/ restart_directory, grid_directory, filter_directory

character(len=256) :: grid_centers_file_prefix, grid_centers_file_suffix, grid_corners_file_prefix, grid_corners_file_suffix
namelist /grid_nml/ grid_centers_file_prefix, grid_centers_file_suffix, grid_corners_file_prefix, grid_corners_file_suffix

integer :: nblocks, nhalos
character(len=256) :: restart_file_prefix, restart_file_middle, restart_file_suffix, filter_input_prefix, filter_input_suffix
namelist /transform_state_nml/ nblocks, nhalos,restart_file_prefix, restart_file_middle, restart_file_suffix, filter_input_prefix, filter_input_suffix

! Assign initialized variables
max_lon = 360.0
min_lon = 0.0
max_lat = 90.0
min_lat = -90.0
corner_indices = [0, 0, 0]
location(:) = 0.0
vertex_latitude(:) = 0.0
vertex_longitude(:) = 0.0
vertex_neighbor_indices(:, :) = 0

call initialize_utilities(progname='create_geometry_file')

call initialize_program()

call assign_triangles_and_quads()

call sort_neighbors()

call check_if_point_in_quad()

call output_triangles_and_quads_to_netcdf()

call finalize_utilities('create_geometry_file')

contains

subroutine initialize_program()

   call find_namelist_in_file('input.nml', 'directory_nml', iunit)
   read(iunit, nml = directory_nml, iostat = io)
   call check_namelist_read(iunit, io, 'directory_nml')

   call find_namelist_in_file('input.nml', 'grid_nml', iunit)
   read(iunit, nml = grid_nml, iostat = io)
   call check_namelist_read(iunit, io, 'grid_nml')

   call find_namelist_in_file('input.nml', 'transform_state_nml', iunit)
   read(iunit, nml = transform_state_nml, iostat = io)
   call check_namelist_read(iunit, io, 'transform_state_nml')

end subroutine initialize_program

subroutine assign_triangles_and_quads()

   type(file_type), allocatable, dimension(:) :: center_files
   type(file_type), allocatable, dimension(:) :: corner_files

   integer :: iblock

   center_files = assign_grid_files_array(nblocks, grid_directory, grid_centers_file_prefix, grid_centers_file_suffix)
   corner_files = assign_grid_files_array(nblocks, grid_directory, grid_corners_file_prefix, grid_corners_file_suffix)

   do iblock = 1, nblocks

      ! Center files

      center_files(iblock)%ncid = nc_open_file_readonly(center_files(iblock)%file_path)

      if (iblock == 1) then
         ! Only call the nc_get_variable_size subroutine once for the center_files because the 
         ! Longitude and Latitude fields must be the same shape
         call nc_get_variable_size(center_files(iblock)%ncid, 'Longitude', variable_size)

         nzs_per_centers_block = variable_size(1)
         nys_per_centers_block = variable_size(2)
         nxs_per_centers_block = variable_size(3)

         truncated_nys_per_centers_block = nys_per_centers_block-2*nhalos
         truncated_nxs_per_centers_block = nxs_per_centers_block-2*nhalos

         allocate(center_altitude_block(nzs_per_centers_block, nys_per_centers_block, nxs_per_centers_block))
         center_altitude_block(:, :, :) = 0
         
         allocate(center_lons_block(nzs_per_centers_block, nys_per_centers_block, nxs_per_centers_block))
         center_lons_block(:, :, :) = 0
         
         allocate(center_lats_block(nzs_per_centers_block, nys_per_centers_block, nxs_per_centers_block))
         center_lats_block(:, :, :) = 0

         allocate(center_lons_truncated(nblocks, truncated_nys_per_centers_block, truncated_nxs_per_centers_block))
         center_lons_truncated(:, :, :) = 0

         allocate(center_lats_truncated(nblocks, truncated_nys_per_centers_block, truncated_nxs_per_centers_block))
         center_lats_truncated(:, :, :) = 0
      end if

      call nc_get_variable(center_files(iblock)%ncid, 'Altitude', center_altitude_block)

      call nc_get_variable(center_files(iblock)%ncid, 'Longitude', center_lons_block)
      center_lons_truncated(iblock, :, :) = center_lons_block(1, nhalos+1:nys_per_centers_block-nhalos, nhalos+1:nxs_per_centers_block-nhalos)

      call nc_get_variable(center_files(iblock)%ncid, 'Latitude', center_lats_block)
      center_lats_truncated(iblock, :, :) = center_lats_block(1, nhalos+1:nys_per_centers_block-nhalos, nhalos+1:nxs_per_centers_block-nhalos)

      ! Corner files      

      corner_files(iblock)%ncid = nc_open_file_readonly(corner_files(iblock)%file_path)

      if (iblock == 1) then
         ! Only call the nc_get_variable_size subroutine once for the corner_files because the 
         ! Longitude Corners and Latitude Corners fields must be the same shape
         call nc_get_variable_size(corner_files(iblock)%ncid, 'Longitude Corners', variable_size)

         nzs_per_corners_block = variable_size(1)
         nys_per_corners_block = variable_size(2)
         nxs_per_corners_block = variable_size(3)

         truncated_nys_per_side_corners_block = nys_per_corners_block-(2*nhalos+2)
         truncated_nxs_per_side_corners_block = nxs_per_corners_block-(2*nhalos+1)

         truncated_nys_per_top_bottom_corners_block = nys_per_corners_block-2*nhalos
         truncated_nxs_per_top_bottom_corners_block = nxs_per_corners_block-2*nhalos

         allocate(corner_lons_block(nzs_per_corners_block, nys_per_corners_block, nxs_per_corners_block))
         corner_lons_block(:, :, :) = 0

         allocate(corner_lats_block(nzs_per_corners_block, nys_per_corners_block, nxs_per_corners_block))
         corner_lats_block(:, :, :) = 0

         allocate(corner_lons_side_truncated(nblocks, truncated_nys_per_side_corners_block, truncated_nxs_per_side_corners_block))
         corner_lons_side_truncated(:, :, :) = 0

         allocate(corner_lats_side_truncated(nblocks, truncated_nys_per_side_corners_block, truncated_nxs_per_side_corners_block))
         corner_lats_side_truncated(:, :, :) = 0

         allocate(corner_lons_top_bottom_truncated(nblocks, truncated_nys_per_top_bottom_corners_block, truncated_nxs_per_top_bottom_corners_block))
         corner_lons_top_bottom_truncated(:, :, :) = 0

         allocate(corner_lats_top_bottom_truncated(nblocks, truncated_nys_per_top_bottom_corners_block, truncated_nxs_per_top_bottom_corners_block))
         corner_lats_top_bottom_truncated(:, :, :) = 0
      end if

      call nc_get_variable(corner_files(iblock)%ncid, 'Longitude Corners', corner_lons_block)
      call nc_get_variable(corner_files(iblock)%ncid, 'Latitude Corners', corner_lats_block)

      if (iblock <= nside_faces) then
         corner_lons_side_truncated(iblock, :, :) = corner_lons_block(1, nhalos+2:nys_per_corners_block-(nhalos+1), nhalos+2:nxs_per_corners_block-nhalos)
         corner_lats_side_truncated(iblock, :, :) = corner_lats_block(1, nhalos+2:nys_per_corners_block-(nhalos+1), nhalos+2:nxs_per_corners_block-nhalos)
      else
         corner_lons_top_bottom_truncated(iblock, :, :) = corner_lons_block(1, nhalos+1:nys_per_corners_block-nhalos, nhalos+1:nxs_per_corners_block-nhalos)
         corner_lats_top_bottom_truncated(iblock, :, :) = corner_lats_block(1, nhalos+1:nys_per_corners_block-nhalos, nhalos+1:nxs_per_corners_block-nhalos)
      end if

   end do

   center_ncols = nblocks*truncated_nys_per_centers_block*truncated_nxs_per_centers_block
   corner_ncols = nside_faces*truncated_nys_per_side_corners_block*truncated_nxs_per_side_corners_block + &
                  (nblocks-nside_faces)*truncated_nys_per_top_bottom_corners_block*truncated_nxs_per_top_bottom_corners_block

   allocate(center_locations(center_ncols))
   allocate(center_latitude(center_ncols))
   allocate(center_longitude(center_ncols))

   min_alt = center_altitude_block(1, 1, 1)

   icol = 0
   do iblock = 1, nblocks
      do iy = 1, truncated_nys_per_centers_block
         do ix = 1, truncated_nxs_per_centers_block
            icol = icol + 1
            ! lon and lat are private attributes of the location type stored internally as radians
            ! and set_location requires the arguments to be in degrees. So, when invoking 
            ! set_location, the aether coordinates in radians must be converted to degrees for the
            ! arguments before they are converted back to radians in the function
            center_locations(icol)%loc = set_location(max(min_lon, min(max_lon, center_lons_truncated(iblock, iy, ix)*RAD2DEG)), max(min_lat, min(max_lat, center_lats_truncated(iblock, iy, ix)*RAD2DEG)), min_alt, VERTISHEIGHT)
            
            location = get_location(center_locations(icol)%loc)
            center_longitude(icol) = location(1)
            center_latitude(icol) = location(2)

         end do
      end do
   end do

   matrix_of_corner_indices = create_matrix_of_corner_indices(nblocks, nvertex_columns, &
                                                              truncated_nys_per_top_bottom_corners_block, &
                                                              truncated_nxs_per_top_bottom_corners_block)

   allocate(distances_from_corners_to_centers(center_ncols))
   allocate(quads(corner_ncols-nvertex_columns))
   allocate(quad_latitude(corner_ncols-nvertex_columns))
   allocate(quad_longitude(corner_ncols-nvertex_columns))
   allocate(quad_neighbor_indices(corner_ncols-nvertex_columns, nquad_neighbors))

   quad_latitude(:) = 0.0
   quad_longitude(:) = 0.0
   quad_neighbor_indices(:, :) = 0

   ivertex = 0
   iquad = 0
   do iblock = 1, nblocks
      if (iblock <= nside_faces) then
         ny = truncated_nys_per_side_corners_block
         nx = truncated_nxs_per_side_corners_block
      else
         ny = truncated_nys_per_top_bottom_corners_block
         nx = truncated_nxs_per_top_bottom_corners_block
      end if
      do iy = 1, ny
         do ix = 1, nx
            corner_indices(1) = iblock
            corner_indices(2) = iy
            corner_indices(3) = ix
            is_a_vertex = is_corner_a_vertex(nvertex_columns, matrix_of_corner_indices, corner_indices)
            if (is_a_vertex .eqv. .true.) then
               ivertex = ivertex + 1
               if (iblock <= nside_faces) then
                  vertexs(ivertex)%loc = set_location(max(min_lon, min(max_lon, corner_lons_side_truncated(iblock, iy, ix)*RAD2DEG)), max(min_lat, min(max_lat, corner_lats_side_truncated(iblock, iy, ix)*RAD2DEG)), min_alt, VERTISHEIGHT)
               else
                  vertexs(ivertex)%loc = set_location(max(min_lon, min(max_lon, corner_lons_top_bottom_truncated(iblock, iy, ix)*RAD2DEG)), max(min_lat, min(max_lat, corner_lats_top_bottom_truncated(iblock, iy, ix)*RAD2DEG)), min_alt, VERTISHEIGHT)
               end if

               location = get_location(vertexs(ivertex)%loc)
               vertex_longitude(ivertex) = location(1)
               vertex_latitude(ivertex) = location(2)
               
            else
               iquad = iquad + 1
               if (iblock <= nside_faces) then
                  quads(iquad)%loc = set_location(max(min_lon, min(max_lon, corner_lons_side_truncated(iblock, iy, ix)*RAD2DEG)), max(min_lat, min(max_lat, corner_lats_side_truncated(iblock, iy, ix)*RAD2DEG)), min_alt, VERTISHEIGHT)
               else
                  quads(iquad)%loc = set_location(max(min_lon, min(max_lon, corner_lons_top_bottom_truncated(iblock, iy, ix)*RAD2DEG)), max(min_lat, min(max_lat, corner_lats_top_bottom_truncated(iblock, iy, ix)*RAD2DEG)), min_alt, VERTISHEIGHT)
               end if

               location = get_location(quads(iquad)%loc)
               quad_longitude(iquad) = location(1)
               quad_latitude(iquad) = location(2)

            end if
         end do
      end do
   end do

   do ivertex = 1, nvertex_columns
      distances_from_corners_to_centers(:) = 0.0

      do icol = 1, center_ncols
         distances_from_corners_to_centers(icol) = get_dist(vertexs(ivertex)%loc, center_locations(icol)%loc)
      end do

      call find_indices_of_n_smallest_elements_of_vector(distances_from_corners_to_centers, vertex_neighbor_indices(ivertex, :))

   end do

   do iquad = 1, corner_ncols-nvertex_columns
      distances_from_corners_to_centers(:) = 0.0

      do icol = 1, center_ncols
         distances_from_corners_to_centers(icol) = get_dist(quads(iquad)%loc, center_locations(icol)%loc)
      end do

      call find_indices_of_n_smallest_elements_of_vector(distances_from_corners_to_centers, quad_neighbor_indices(iquad, :))

   end do

end subroutine assign_triangles_and_quads

subroutine output_triangles_and_quads_to_netcdf()

   type(file_type) :: geometry_file
   integer :: dim_id_center_altitudes
   integer :: dim_id_vertex_columns, dim_id_quad_columns, dim_id_center_columns
   integer :: dim_id_vertex_neighbors, dim_id_quad_neighbors
   integer :: var_id_center_altitude
   integer :: var_id_center_longitude, var_id_center_latitude
   integer :: var_id_vertex_longitude, var_id_vertex_latitude, var_id_vertex
   integer :: var_id_quad_longitude, var_id_quad_latitude, var_id_quad
   integer :: xtype

   xtype = 5

   geometry_file%file_path = trim(filter_directory) // 'geometry_file.nc'

   geometry_file%ncid =  nc_create_file(geometry_file%file_path)

   geometry_file%ncstatus = nf90_def_dim(geometry_file%ncid, 'center_altitudes', nzs_per_centers_block, dim_id_center_altitudes)

   geometry_file%ncstatus = nf90_def_dim(geometry_file%ncid, 'vertex_columns', nvertex_columns, dim_id_vertex_columns)

   geometry_file%ncstatus = nf90_def_dim(geometry_file%ncid, 'vertex_neighbors', nvertex_neighbors, dim_id_vertex_neighbors)

   geometry_file%ncstatus = nf90_def_dim(geometry_file%ncid, 'quad_columns', (corner_ncols-nvertex_columns), dim_id_quad_columns)

   geometry_file%ncstatus = nf90_def_dim(geometry_file%ncid, 'quad_neighbors', nquad_neighbors, dim_id_quad_neighbors)

   geometry_file%ncstatus = nf90_def_dim(geometry_file%ncid, 'center_columns', center_ncols, dim_id_center_columns)

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'center_altitude', xtype, dim_id_center_altitudes, var_id_center_altitude)

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'vertex_longitude', xtype, dim_id_vertex_columns, var_id_vertex_longitude)

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'vertex_latitude', xtype, dim_id_vertex_columns, var_id_vertex_latitude)

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'quad_longitude', xtype, dim_id_quad_columns, var_id_quad_longitude)

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'quad_latitude', xtype, dim_id_quad_columns, var_id_quad_latitude)

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'center_longitude', xtype, dim_id_center_columns, var_id_center_longitude)

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'center_latitude', xtype, dim_id_center_columns, var_id_center_latitude)

   xtype = 4

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'vertex_neighbor_indices', xtype, [dim_id_vertex_columns, dim_id_vertex_neighbors], var_id_vertex)

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'quad_neighbor_indices', xtype, [dim_id_quad_columns, dim_id_quad_neighbors], var_id_quad)

   ! Add the indices of the north and south pole quad
   geometry_file%ncstatus = nf90_put_att(geometry_file%ncid, NF90_GLOBAL, 'index_of_north_pole_quad_column', inorth_pole_quad_column)
   geometry_file%ncstatus = nf90_put_att(geometry_file%ncid, NF90_GLOBAL, 'index_of_south_pole_quad_column', isouth_pole_quad_column)

   call nc_end_define_mode(geometry_file%ncid)

   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_center_altitude, center_altitude_block(:, 1, 1))

   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_center_longitude, center_longitude)
   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_center_latitude, center_latitude)

   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_vertex_longitude, vertex_longitude)
   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_vertex_latitude, vertex_latitude)
   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_vertex, vertex_neighbor_indices)

   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_quad_longitude, quad_longitude)
   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_quad_latitude, quad_latitude)
   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_quad, quad_neighbor_indices)

   call nc_close_file(geometry_file%ncid)

end subroutine output_triangles_and_quads_to_netcdf

function assign_grid_files_array(nblocks, grid_directory, grid_file_prefix, grid_file_suffix) result(grid_files)

   integer, intent(in)          :: nblocks
   character(len=*), intent(in) :: grid_directory
   character(len=*), intent(in) :: grid_file_prefix
   character(len=*), intent(in) :: grid_file_suffix
   type(file_type), allocatable, dimension(:) :: grid_files
   character(len=4) :: block_name
   integer :: iblock

   print *, 'nblocks'
   print *, nblocks

   allocate(grid_files(nblocks))

   do iblock = 1, nblocks
      block_name = zero_fill(integer_to_string(iblock-1), 4)
      grid_files(iblock)%file_path = trim(grid_directory) // trim(grid_file_prefix) // &
                                      block_name // trim(grid_file_suffix)
      print *, grid_files(iblock)%file_path
   end do

end function assign_grid_files_array

subroutine find_indices_of_n_smallest_elements_of_vector(vector, smallest_indices)
   real(r8), dimension(:), intent(in) :: vector
   integer, dimension(:), intent(inout) :: smallest_indices
   integer :: iindex, nindices
   logical, allocatable, dimension(:) :: vector_mask

   nindices = size(smallest_indices)
   allocate(vector_mask(size(vector)))
   vector_mask(:) = .true.

   do iindex=1, nindices
      ! This peculiar (iindex:iindex) array indexing is required because minloc returns an array of dim 1
      smallest_indices(iindex:iindex) = minloc(vector, mask=vector_mask)
      vector_mask(smallest_indices(iindex)) = .false.
   end do

end subroutine find_indices_of_n_smallest_elements_of_vector

function create_matrix_of_corner_indices(nblocks, nvertex_columns, truncated_nys_per_corners_block, &
                                         truncated_nxs_per_corners_block) result(matrix_of_corner_indices)
   integer, intent(in) :: nblocks
   integer, intent(in) :: nvertex_columns
   integer, intent(in) :: truncated_nys_per_corners_block
   integer, intent(in) :: truncated_nxs_per_corners_block
   integer, dimension(nvertex_columns, 3) :: matrix_of_corner_indices

   matrix_of_corner_indices(:, :) = 0

   matrix_of_corner_indices(1, :) = [nblocks-1, 1, 1]
   matrix_of_corner_indices(2, :) = [nblocks-1, 1, truncated_nxs_per_corners_block]
   matrix_of_corner_indices(3, :) = [nblocks-1, truncated_nys_per_corners_block, 1]
   matrix_of_corner_indices(4, :) = [nblocks-1, truncated_nys_per_corners_block, truncated_nxs_per_corners_block]
   matrix_of_corner_indices(5, :) = [nblocks, 1, 1]
   matrix_of_corner_indices(6, :) = [nblocks, 1, truncated_nxs_per_corners_block]
   matrix_of_corner_indices(7, :) = [nblocks, truncated_nys_per_corners_block, 1]
   matrix_of_corner_indices(8, :) = [nblocks, truncated_nys_per_corners_block, truncated_nxs_per_corners_block]

end function create_matrix_of_corner_indices

function is_corner_a_vertex(nvertex_columns, matrix_of_corner_indices, corner_indices) result(is_a_vertex)
   integer, intent(in)                          :: nvertex_columns
   integer, dimension(nvertex_columns, 3), intent(in) :: matrix_of_corner_indices
   integer, dimension(3), intent(in)            :: corner_indices
   logical                                      :: is_a_vertex
   integer                                      :: ivertex

   is_a_vertex = .false.
   do ivertex = 1, nvertex_columns
      if (all(matrix_of_corner_indices(ivertex, :) == corner_indices)) then
         is_a_vertex = .true.
         exit
      end if
   end do

end function is_corner_a_vertex

subroutine check_if_point_in_quad()

   real(r8), dimension(4) :: x_neighbors, y_neighbors
   logical  :: inside
   integer  :: num_outside, num_inside, num_poles
   integer  :: icolumn, ineighbor
   real(r8) :: tolerance

   tolerance = 0.00000314159_r8

   num_outside = 0
   num_inside = 0
   num_poles = 0

   do icolumn=1, (corner_ncols-nvertex_columns)
      x_neighbors(:) = 0.0
      y_neighbors(:) = 0.0
      do ineighbor=1, nquad_neighbors
         x_neighbors(ineighbor) = center_longitude(quad_neighbor_indices(icolumn, ineighbor))
         y_neighbors(ineighbor) = center_latitude(quad_neighbor_indices(icolumn, ineighbor))
      end do

      inside = in_quad(quad_longitude(icolumn), quad_latitude(icolumn), x_neighbors, y_neighbors)

      if (inside) then
         num_inside = num_inside + 1
      else if (abs(-90.0-quad_latitude(icolumn)) < tolerance) then
         isouth_pole_quad_column = icolumn
         num_poles = num_poles + 1
      else if (abs(90.0-quad_latitude(icolumn)) < tolerance) then
         inorth_pole_quad_column = icolumn
         num_poles = num_poles + 1
      else
         num_outside = num_outside + 1
      end if
   end do

end subroutine check_if_point_in_quad

subroutine sort_neighbors()

   real(r8), dimension(nquad_neighbors)    :: atan2_results
   real(r8), dimension(nquad_neighbors)    :: offset_lons, offset_lats
   integer,  dimension(nquad_neighbors)    :: sorted_array
   real(r8)  :: central_lon, central_lat
   integer(i8), dimension(4) :: indices
   integer(i8) :: nneighbors
   integer :: ineighbor
   logical :: in_first_quadrant, in_fourth_quadrant

   ! Need to declare an additional integer equal to 4, since the interface defined for index sort,
   ! when passed r8 values, requires an i8 length
   nneighbors = nquad_neighbors
   atan2_results(:) = 0

   do iquad = 1, (corner_ncols-nvertex_columns)

      central_lon = quad_longitude(iquad)
      central_lat = quad_latitude(iquad)

      ! Assign the offset_lons and offset_lats elements
      do ineighbor = 1, nquad_neighbors
         offset_lats(ineighbor) = center_latitude(quad_neighbor_indices(iquad, ineighbor))
         offset_lons(ineighbor) = center_longitude(quad_neighbor_indices(iquad, ineighbor)) 
      end do

      ! Check if the quad spans the prime meridian. It's definitely possible to perform this check
      ! with only a single logical variable but using two makes the code more readable.
      in_first_quadrant = .false.
      in_fourth_quadrant = .false.

      do ineighbor = 1, nquad_neighbors
         if (offset_lons(ineighbor) >= 0.0 .and. offset_lons(ineighbor) <= 90.0) then
            in_first_quadrant = .true.
         else if (offset_lons(ineighbor) >= 270.0 .and. offset_lons(ineighbor) <= 360.0) then
            in_fourth_quadrant = .true.
         end if
      end do

      if (in_first_quadrant .and. in_fourth_quadrant) then
         central_lon = central_lon + 180.0
         if (central_lon >= 360.0) then
            central_lon = central_lon - 360.0
         end if
         ! The quad spans the prime meridian
         do ineighbor = 1, nquad_neighbors
            offset_lons(ineighbor) = offset_lons(ineighbor) + 180.0
            if (offset_lons(ineighbor) >= 360.0) then
               offset_lons(ineighbor) = offset_lons(ineighbor) - 360.0
            end if
         end do
      end if

      do ineighbor = 1, nquad_neighbors
         offset_lons(ineighbor) = (central_lon-offset_lons(ineighbor))*DEG2RAD
         offset_lats(ineighbor) = (central_lat-offset_lats(ineighbor))*DEG2RAD
         atan2_results(ineighbor) = atan2(offset_lats(ineighbor), offset_lons(ineighbor))      
      end do
      
      call index_sort(atan2_results, indices, nneighbors)

      do ineighbor = 1, nneighbors
         sorted_array(ineighbor) = quad_neighbor_indices(iquad, indices(ineighbor))
      end do

      quad_neighbor_indices(iquad, :) = sorted_array(:)

   end do

end subroutine sort_neighbors

end program create_geometry_file
