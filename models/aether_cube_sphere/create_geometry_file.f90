program create_geometry_file

use               netcdf
use            types_mod, only : r8, RAD2DEG
use         location_mod, only : location_type, get_dist, get_location, set_location, VERTISHEIGHT
use        utilities_mod, only : initialize_utilities, finalize_utilities
use netcdf_utilities_mod, only : nc_open_file_readonly, nc_get_variable_size, nc_get_variable, &
                                 nc_create_file, nc_end_define_mode, nc_close_file
use  transform_state_mod, only : read_namelist, file_type, integer_to_string, zero_fill, nblocks, &
                                 nhalos, dart_directory, grid_directory, grid_centers_file_prefix, &
                                 grid_centers_file_suffix, grid_corners_file_prefix, &
                                 grid_corners_file_suffix

implicit none

! There will always be exactly 8 vertices on a cube sphere, regardless of how many blocks
integer, parameter               :: nvertices = 8
integer, parameter               :: nside_faces = 4
integer, parameter               :: nvertex_triangle_neighbors = 3
integer, parameter               :: ncube_face_quad_neighbors = 4
integer, dimension(nvertices, 3) :: matrix_of_corner_indices
logical                          :: is_a_cube_vertex
integer, dimension(3)            :: corner_indices

! Declare variables for the center and corner netcdf files

type(file_type), allocatable, dimension(:) :: center_files, corner_files
integer, dimension(3) :: variable_size

integer :: nzs_per_centers_block, nys_per_centers_block, nxs_per_centers_block, &
           nzs_per_corners_block, nys_per_corners_block, nxs_per_corners_block

integer :: truncated_nys_per_centers_block, truncated_nxs_per_centers_block, &
           truncated_nys_per_side_corners_block, truncated_nxs_per_side_corners_block, &
           truncated_nys_per_top_bottom_corners_block, truncated_nxs_per_top_bottom_corners_block

integer :: center_ncols, corner_ncols

real(r8), allocatable, dimension(:, :, :) :: center_lats_block, center_lons_block, & 
                                             corner_lats_block, corner_lons_block

real(r8), allocatable, dimension(:, :, :) :: center_lats_truncated, center_lons_truncated, &
                                             corner_lats_side_truncated, corner_lons_side_truncated, &
                                             corner_lats_top_bottom_truncated, corner_lons_top_bottom_truncated

real(r8), allocatable, dimension(:)       :: center_lats_vector, center_lons_vector, &
                                             cube_face_quad_lats_vector, cube_face_quad_lons_vector

integer, allocatable, dimension(:, :)     :: cube_face_quad_neighbor_indices

integer, dimension(nvertices, nvertex_triangle_neighbors) :: vertex_triangle_neighbor_indices

real(r8), dimension(nvertices)            :: vertex_triangle_lats_vector, vertex_triangle_lons_vector

real(r8), dimension(3)                    :: location

! Declare variables for the geometry_file that will be output by this program
real(r8)                         :: alt
real(r8)                         :: max_lon
real(r8)                         :: min_lon
real(r8)                         :: max_lat
real(r8)                         :: min_lat

type :: vertex_triangle
   type(location_type) :: loc
   integer, dimension(3) :: column_index_of_neighbors
end type vertex_triangle

type :: cube_face_quad
   type(location_type) :: loc
   logical  :: spans_prime_meridian
   logical  :: encloses_pole
   integer, dimension(4) :: column_index_of_neighbors
end type cube_face_quad

! This is an intentionally unecessary nesting of a derived type just so that the lat and lon of each
! center_location is accessed using the same %loc attibute syntax
type :: center_location
   type(location_type) :: loc
end type center_location

type(vertex_triangle), dimension(nvertices) :: vertex_triangles
type(cube_face_quad), allocatable, dimension(:) :: cube_face_quads
type(center_location), allocatable, dimension(:) :: center_locations

integer                          :: ivertex_triangle
integer                          :: icube_face_quad

integer :: icol, ix, iy, ncols, nx, ny
real(r8), allocatable, dimension(:) :: distances_from_corners_to_centers

! Assign initialized variables
alt = 96343.08
max_lon = 360.0
min_lon = 0.0
max_lat = 90.0
min_lat = -90.0
corner_indices = [0, 0, 0]
location(:) = 0.0
vertex_triangle_lats_vector(:) = 0.0
vertex_triangle_lons_vector(:) = 0.0
vertex_triangle_neighbor_indices(:, :) = 0

call initialize_utilities(progname='create_geometry_file')

call read_namelist('transform_state_nml')

call assign_triangles_and_quads()

call output_triangles_and_quads_to_netcdf()

call finalize_utilities('create_geometry_file')

contains

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
         
         allocate(center_lons_block(nzs_per_centers_block, nys_per_centers_block, nxs_per_centers_block))
         center_lons_block(:, :, :) = 0
         
         allocate(center_lats_block(nzs_per_centers_block, nys_per_centers_block, nxs_per_centers_block))
         center_lats_block(:, :, :) = 0

         allocate(center_lons_truncated(nblocks, truncated_nys_per_centers_block, truncated_nxs_per_centers_block))
         center_lons_truncated(:, :, :) = 0

         allocate(center_lats_truncated(nblocks, truncated_nys_per_centers_block, truncated_nxs_per_centers_block))
         center_lats_truncated(:, :, :) = 0
      end if

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
   allocate(center_lats_vector(center_ncols))
   allocate(center_lons_vector(center_ncols))

   icol = 0
   do iblock = 1, nblocks
      do iy = 1, truncated_nys_per_centers_block
         do ix = 1, truncated_nxs_per_centers_block
            icol = icol + 1
            ! lon and lat are private attributes of the location type stored internally as radians
            ! and set_location requires the arguments to be in degrees. So, when invoking 
            ! set_location, the aether coordinates in radians must be converted to degrees for the
            ! arguments before they are converted back to radians in the function
            center_locations(icol)%loc = set_location(max(min_lon, min(max_lon, center_lons_truncated(iblock, iy, ix)*RAD2DEG)), max(min_lat, min(max_lat, center_lats_truncated(iblock, iy, ix)*RAD2DEG)), alt, VERTISHEIGHT)
            
            location = get_location(center_locations(icol)%loc)
            center_lons_vector(icol) = location(1)
            center_lats_vector(icol) = location(2)

         end do
      end do
   end do

   matrix_of_corner_indices = create_matrix_of_corner_indices(nblocks, nvertices, &
                                                              truncated_nys_per_top_bottom_corners_block, &
                                                              truncated_nxs_per_top_bottom_corners_block)

   allocate(distances_from_corners_to_centers(center_ncols))
   allocate(cube_face_quads(corner_ncols-nvertices))
   allocate(cube_face_quad_lats_vector(corner_ncols-nvertices))
   allocate(cube_face_quad_lons_vector(corner_ncols-nvertices))
   allocate(cube_face_quad_neighbor_indices(corner_ncols-nvertices, ncube_face_quad_neighbors))

   cube_face_quad_lats_vector(:) = 0.0
   cube_face_quad_lons_vector(:) = 0.0
   cube_face_quad_neighbor_indices(:, :) = 0

   ivertex_triangle = 0
   icube_face_quad = 0
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
            is_a_cube_vertex = is_corner_a_cube_vertex(nvertices, matrix_of_corner_indices, corner_indices)
            if (is_a_cube_vertex .eqv. .true.) then
               ivertex_triangle = ivertex_triangle + 1
               if (iblock <= nside_faces) then
                  vertex_triangles(ivertex_triangle)%loc = set_location(max(min_lon, min(max_lon, corner_lons_side_truncated(iblock, iy, ix)*RAD2DEG)), max(min_lat, min(max_lat, corner_lats_side_truncated(iblock, iy, ix)*RAD2DEG)), alt, VERTISHEIGHT)
               else
                  vertex_triangles(ivertex_triangle)%loc = set_location(max(min_lon, min(max_lon, corner_lons_top_bottom_truncated(iblock, iy, ix)*RAD2DEG)), max(min_lat, min(max_lat, corner_lats_top_bottom_truncated(iblock, iy, ix)*RAD2DEG)), alt, VERTISHEIGHT)
               end if

               location = get_location(vertex_triangles(ivertex_triangle)%loc)
               vertex_triangle_lons_vector(ivertex_triangle) = location(1)
               vertex_triangle_lats_vector(ivertex_triangle) = location(2)
               
            else
               icube_face_quad = icube_face_quad + 1
               if (iblock <= nside_faces) then
                  cube_face_quads(icube_face_quad)%loc = set_location(max(min_lon, min(max_lon, corner_lons_side_truncated(iblock, iy, ix)*RAD2DEG)), max(min_lat, min(max_lat, corner_lats_side_truncated(iblock, iy, ix)*RAD2DEG)), alt, VERTISHEIGHT)
               else
                  cube_face_quads(icube_face_quad)%loc = set_location(max(min_lon, min(max_lon, corner_lons_top_bottom_truncated(iblock, iy, ix)*RAD2DEG)), max(min_lat, min(max_lat, corner_lats_top_bottom_truncated(iblock, iy, ix)*RAD2DEG)), alt, VERTISHEIGHT)
               end if

               location = get_location(cube_face_quads(icube_face_quad)%loc)
               cube_face_quad_lons_vector(icube_face_quad) = location(1)
               cube_face_quad_lats_vector(icube_face_quad) = location(2)

            end if
         end do
      end do
   end do

   do ivertex_triangle = 1, nvertices
      distances_from_corners_to_centers(:) = 0.0

      do icol = 1, center_ncols
         distances_from_corners_to_centers(icol) = get_dist(vertex_triangles(ivertex_triangle)%loc, center_locations(icol)%loc)
      end do

      call find_indices_of_n_smallest_elements_of_vector(distances_from_corners_to_centers, vertex_triangle_neighbor_indices(ivertex_triangle, :))

   end do

   do icube_face_quad = 1, corner_ncols-nvertices
      distances_from_corners_to_centers(:) = 0.0

      do icol = 1, center_ncols
         distances_from_corners_to_centers(icol) = get_dist(cube_face_quads(icube_face_quad)%loc, center_locations(icol)%loc)
      end do

      call find_indices_of_n_smallest_elements_of_vector(distances_from_corners_to_centers, cube_face_quad_neighbor_indices(icube_face_quad, :))

   end do

end subroutine assign_triangles_and_quads

subroutine output_triangles_and_quads_to_netcdf()

   type(file_type) :: geometry_file
   integer :: dim_id_vertex_columns, dim_id_cube_face_quad_columns, dim_id_center_columns
   integer :: dim_id_vertex_neighbors, dim_id_cube_face_quad_neighbors
   integer :: var_id_center_longitude, var_id_center_latitude
   integer :: var_id_vertex_triangle_longitude, var_id_vertex_triangle_latitude, var_id_vertex_triangle
   integer :: var_id_cube_face_quad_longitude, var_id_cube_face_quad_latitude, var_id_cube_face_quad
   integer :: xtype

   xtype = 5

   geometry_file%file_path = trim(dart_directory) // 'geometry_file.nc'

   geometry_file%ncid =  nc_create_file(geometry_file%file_path)

   geometry_file%ncstatus = nf90_def_dim(geometry_file%ncid, 'vertex_columns', nvertices, dim_id_vertex_columns)

   geometry_file%ncstatus = nf90_def_dim(geometry_file%ncid, 'vertex_neighbors', nvertex_triangle_neighbors, dim_id_vertex_neighbors)

   geometry_file%ncstatus = nf90_def_dim(geometry_file%ncid, 'cube_face_quad_columns', (corner_ncols-nvertices), dim_id_cube_face_quad_columns)

   geometry_file%ncstatus = nf90_def_dim(geometry_file%ncid, 'cube_face_quad_neighbors', ncube_face_quad_neighbors, dim_id_cube_face_quad_neighbors)

   geometry_file%ncstatus = nf90_def_dim(geometry_file%ncid, 'center_columns', center_ncols, dim_id_center_columns)

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'vertex_triangle_longitude', xtype, dim_id_vertex_columns, var_id_vertex_triangle_longitude)

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'vertex_triangle_latitude', xtype, dim_id_vertex_columns, var_id_vertex_triangle_latitude)

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'cube_face_quad_longitude', xtype, dim_id_cube_face_quad_columns, var_id_cube_face_quad_longitude)

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'cube_face_quad_latitude', xtype, dim_id_cube_face_quad_columns, var_id_cube_face_quad_latitude)

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'center_longitude', xtype, dim_id_center_columns, var_id_center_longitude)

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'center_latitude', xtype, dim_id_center_columns, var_id_center_latitude)

   xtype = 4

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'vertex_triangle_neighbor_indices', xtype, [dim_id_vertex_columns, dim_id_vertex_neighbors], var_id_vertex_triangle)

   geometry_file%ncstatus = nf90_def_var(geometry_file%ncid, 'cube_face_quad_neighbor_indices', xtype, [dim_id_cube_face_quad_columns, dim_id_cube_face_quad_neighbors], var_id_cube_face_quad)

   call nc_end_define_mode(geometry_file%ncid)

   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_center_longitude, center_lons_vector)
   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_center_latitude, center_lats_vector)

   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_vertex_triangle_longitude, vertex_triangle_lons_vector)
   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_vertex_triangle_latitude, vertex_triangle_lats_vector)
   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_vertex_triangle, vertex_triangle_neighbor_indices)

   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_cube_face_quad_longitude, cube_face_quad_lons_vector)
   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_cube_face_quad_latitude, cube_face_quad_lats_vector)
   geometry_file%ncstatus = nf90_put_var(geometry_file%ncid, var_id_cube_face_quad, cube_face_quad_neighbor_indices)

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

function create_matrix_of_corner_indices(nblocks, nvertices, truncated_nys_per_corners_block, &
                                         truncated_nxs_per_corners_block) result(matrix_of_corner_indices)
   integer, intent(in) :: nblocks
   integer, intent(in) :: nvertices
   integer, intent(in) :: truncated_nys_per_corners_block
   integer, intent(in) :: truncated_nxs_per_corners_block
   integer, dimension(nvertices, 3) :: matrix_of_corner_indices

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

function is_corner_a_cube_vertex(nvertices, matrix_of_corner_indices, corner_indices) result(is_a_cube_vertex)
   integer, intent(in)                          :: nvertices
   integer, dimension(nvertices, 3), intent(in) :: matrix_of_corner_indices
   integer, dimension(3), intent(in)            :: corner_indices
   logical                                      :: is_a_cube_vertex
   integer                                      :: ivertex

   is_a_cube_vertex = .false.
   do ivertex = 1, nvertices
      if (all(matrix_of_corner_indices(ivertex, :) == corner_indices)) then
         is_a_cube_vertex = .true.
         exit
      end if
   end do

end function

end program create_geometry_file
