! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program test_aether_grid

use utilities_mod,             only : initialize_utilities, finalize_utilities

use assim_model_mod,            only : static_init_assim_model

use types_mod,                  only : r8, i8, DEG2RAD, RAD2DEG

use location_mod,               only : location_type, get_location

use utilities_mod,              only : error_handler, E_ERR, E_MSG

use cube_sphere_grid_tools_mod, only : is_point_in_triangle, is_point_in_quad, grid_to_lat_lon, &
                                       lat_lon_to_xyz, col_index_to_lat_lon, lat_lon_to_grid,   &
                                       get_bounding_box, lat_lon_to_col_index, get_grid_delta

! Need basic grid description read from template file by model_mod
use model_mod,                  only : np, ncenter_altitudes, get_state_index, get_state_meta_data

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "test_aether_grid"


call initialize_utilities(source)

call static_init_assim_model()

call test_grid_box()

call finalize_utilities

!------------------------------------------------------------

contains


! Testing subroutine for grid definition and interpolation tools
! This is not designed to be run with more than one process

subroutine test_grid_box
                              
integer             :: i, j, num_bound_points, qty, lon_count, lat_count
integer             :: grid_face(4), grid_lat_ind(4), grid_lon_ind(4)
integer             :: my_face, my_level, my_qty, my_lon_ind, my_lat_ind
integer             :: test_face, test_lat_ind, test_lon_ind, col_index, test_col_index
integer             :: num_test_lats, num_test_lons
integer(i8)         :: state_index
real(r8)            :: del, half_del       ! Grid row spacing and half of that
real(r8)            :: pt_lon_d, pt_lat_d, pt_lon, pt_lat
real(r8)            :: qxyz(4, 3), pxyz(3), grid_pt_lat(4), grid_pt_lon(4)
real(r8)            :: lon_lat_hgt(3), my_lat, my_lon, base_dist, dist_sum
logical             :: inside

type(location_type) :: location

! Error message strings                
character(len=512) :: string1, string2
          
! Get the rest of the geometry, np read by static_init_model
call get_grid_delta(np, del, half_del)

! Test that grid_to_lat_lon and lat_lon_to_grid are inverses of each other
do my_face = 0, 5
   do my_lat_ind = 1, np
      do my_lon_ind = 1, np
         call grid_to_lat_lon(my_face, my_lat_ind, my_lon_ind, del, half_del, pt_lat, pt_lon)
         call lat_lon_to_grid(pt_lat, pt_lon, del, half_del, test_face, test_lat_ind, test_lon_ind)
         if(my_face /= test_face .or. my_lat_ind /= test_lat_ind .or. my_lon_ind /= test_lon_ind) then
            write(string1, *) 'Test failed: lat_lon_to_grid is not inverse of grid_to_lat_lon'
            write(string2, *) my_face, test_face, my_lat_ind, test_lat_ind, my_lon_ind, test_lon_ind
            call error_handler(E_ERR, 'test_grid_box', string1, source, text2=string2)
         endif

         ! Test that col_index_to_lat_lon and lat_lon_to_col_index are inverses of each other
         col_index = my_lon_ind + (my_lat_ind - 1) * np + my_face * np*np
         call col_index_to_lat_lon(col_index, np, del, half_del, pt_lat, pt_lon)
         test_col_index = lat_lon_to_col_index(pt_lat, pt_lon, del, half_del, np)
         if(col_index /= test_col_index) then
            write(string1, *) 'Test failed: lat_lon_to_col_index is not inverse of col_index_to_lat_lon'
            write(string2, *) my_face, my_lat_ind, my_lon_ind, col_index, test_col_index
            call error_handler(E_ERR, 'test_grid_box', string1, source, text2=string2)
         endif
      enddo
   enddo
enddo

! Test points for the following:
! 1. Does the bounding box found contain the observed point?
! 2. Are the computed vertex latitudes and longitudes the same as those in the Aether grid files?

! Largest edges are on the quads in the center of a face
! Get distance, base_dist, along side of center quad
do i = 1, 2
   ! Traverse half of the rows (each np across), plus halfway across the next row
   state_index = (np/2)*np + np/2 + i - 1
   call get_state_meta_data(state_index, location, qty)
   lon_lat_hgt = get_location(location)
   qxyz(i, 1:3) = lat_lon_to_xyz(DEG2RAD*lon_lat_hgt(2), DEG2RAD*lon_lat_hgt(1))
enddo
base_dist = sqrt(sum((qxyz(1, :) - qxyz(2, :))**2))

! Loop through many longitude and latitude points for testing
num_test_lons = 3600
do lon_count = 0, num_test_lons
   pt_lon_d = lon_count * (360.0_r8 / num_test_lons)
   num_test_lats = 1800
      do lat_count = -900, 900
      pt_lat_d = lat_count * (180.0_r8 / num_test_lats)

      ! Convert to radians
      pt_lon = DEG2RAD * pt_lon_d
      pt_lat = DEG2RAD * pt_lat_d

      ! Get the x, y, z coords for this point
      pxyz = lat_lon_to_xyz(pt_lat, pt_lon);

      call get_bounding_box(pt_lat, pt_lon, del, half_del, np, &
         grid_face, grid_lat_ind, grid_lon_ind, grid_pt_lat, grid_pt_lon, num_bound_points)

      do i = 1, num_bound_points
         ! Convert to x, y, z coords to check for whether points are in tri/quad
         qxyz(i, 1:3) = lat_lon_to_xyz(grid_pt_lat(i), grid_pt_lon(i));
      enddo

      ! Get latitude longitude of bounding points from get_state_meta_data as confirmation test
      do i = 1, num_bound_points
         state_index = get_state_index(grid_face(i), grid_lat_ind(i), grid_lon_ind(i), 1, 1)
         call get_state_meta_data(state_index, location, qty)
         lon_lat_hgt = get_location(location)
         ! Deal with Aether file round off
         if(abs(lon_lat_hgt(1) - 360.0_r8) < 0.0001) lon_lat_hgt(1) = 0.0_r8
         if(abs(RAD2DEG*grid_pt_lat(i) - lon_lat_hgt(2)) > 0.0001_r8 .or. &
            abs(RAD2DEG*grid_pt_lon(i) - lon_lat_hgt(1)) > 0.0001_r8) then
            write(string1, *) 'Test failed: Aether files grid points inconsistent with get_state_meta_data'
            write(string2, *) grid_pt_lat(i), grid_pt_lon(i), lon_lat_hgt(2), lon_lat_hgt(1)
            call error_handler(E_ERR, 'test_grid_box', string1, source, text2=string2)
         endif
      enddo

      if(num_bound_points == 3) then
         ! See if the point is inside a local approximately tangent triangle
         inside = is_point_in_triangle(qxyz(1, :), qxyz(2, :), qxyz(3, :), pxyz);
      else
         ! Or quadrilateral
         inside = is_point_in_quad(qxyz, pxyz);
      endif

      if(.not. inside) then
         write(string1, *) 'Test failed: Point is not inside the triangle or quadrilateral'
         write(string2, *) pt_lat, pt_lon, num_bound_points
         call error_handler(E_ERR, 'test_grid_box', string1, source, text2=string2)
      endif

      ! Also check on distance to vertices; this greatly reduces the possibility that
      ! bounding boxes that are bigger than they should be are being found
      dist_sum = 0.0_r8
      do i = 1, num_bound_points
         ! Compute sum of distances between point and each of the vertices
         dist_sum = dist_sum + sqrt(sum((qxyz(i, :) - pxyz)**2))
      end do

      if(num_bound_points == 4) then
         ! For quad, sum should be less than 3.5 times the baseline
         if(dist_sum / base_dist > 3.5_r8) then
            write(string1, *) 'Test failed: Ratio of sum of distances to vertices is too large for quad'
            ! Additional info that could be helpful
            !!!do i = 1, num_bound_points
               !!!write(*, *) 'grid ', i, grid_pt_lat(i), grid_pt_lon(i)
               !!!write(*, *) 'grid xyz ', i, qxyz(i, :) 
            !!!enddo
            write(string2, *) 'point ', pt_lat, pt_lon, 'point xyz ', pxyz
            call error_handler(E_ERR, 'test_grid_box', string1, source, text2=string2)
         endif
      elseif(num_bound_points == 3) then
         ! For triangle, sum should be less than 3 times the baseline
         if(dist_sum / base_dist > 3.0_r8) then
            write(string1, *) 'Test failed: ratio of sum of distances to vertices is too large for triangle'
            call error_handler(E_ERR, 'test_grid_box', string1, source)
         endif
      endif
   enddo
enddo

!-------------------
! Block that loops through all state variables and confirms that the algorithms for mapping
! state vector (face/lon/lat) indices and get_state_meta_data correcty match up.
do my_qty = 1, 2
   do my_level = 1, ncenter_altitudes
      do my_face = 0, 5
         do my_lat_ind = 1, np
            do my_lon_ind = 1, np
               state_index = get_state_index(my_face, my_lat_ind, my_lon_ind, &
                  my_level, my_qty)
               call get_state_meta_data(state_index, location, qty)
               lon_lat_hgt = get_location(location)

               ! Want to compare the lat lon directly from code to that from get_state_meta_data
               call grid_to_lat_lon(my_face, my_lat_ind, my_lon_ind, del, half_del, my_lat, my_lon)

               ! ROUNDOFF FROM AETHER
               if(abs(lon_lat_hgt(1) - 360.0_r8) < 0.0001) lon_lat_hgt(1) = 0.0_r8

               ! Check that things are consistent            
               if(abs(RAD2DEG*my_lat - lon_lat_hgt(2)) > 0.0001_r8 .or. &
                  abs(RAD2DEG*my_lon - lon_lat_hgt(1)) > 0.0001_r8) then
                  write(string1, *) 'Test Failed: Grid points not appropriately mapping'
                  write(string2, *) my_face, my_qty, my_level, my_lat_ind, my_lon_ind
                  call error_handler(E_ERR, 'test_grid_box', string1, source)
               endif

            enddo
         enddo
      enddo
   enddo
enddo

write(string1, *) 'ALL TESTS PASSED'
call error_handler(E_MSG, 'test_grid_box', string1, source)

end subroutine test_grid_box

!-----------------------------------------------------------------------

end program test_aether_grid
