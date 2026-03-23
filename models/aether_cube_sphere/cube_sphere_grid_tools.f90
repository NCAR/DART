! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module cube_sphere_grid_tools_mod

! This module provides tools that know about the horizontal geometry and relation to 
! storage patterns for the Aether cube sphere grid. 

use types_mod,     only : r8, PI

use utilities_mod, only : error_handler, E_ERR

implicit none
private

public :: lat_lon_to_col_index, col_index_to_lat_lon, get_bounding_box,      &
          is_point_in_triangle, is_point_in_quad, lat_lon_to_xyz,            & 
          grid_to_lat_lon, lat_lon_to_grid, get_face, fix_face, get_corners, &
          get_grid_delta

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "cube_sphere_grid_tools"

contains

!------------------------------------------------------------------

! Compute the spacing between grid rows on a face

subroutine get_grid_delta(np, del, half_del)

integer,  intent(in)  :: np
real(r8), intent(out) :: del, half_del

real(r8) :: cube_side

! Cube side is divided into np-1 interior intervals of width 2sqrt(1/3) / np and
! two exterior intervals of half  width, sqrt(1/3) / np 
cube_side = 2.0_r8 * sqrt(1.0_r8 / 3.0_r8)
! Get the spacing of the grid points
del = cube_side / np  
half_del = del / 2.0_r8

end subroutine get_grid_delta

!------------------------------------------------------------------

! Given the face (from 0 to 5) the number of lons/lats on a face (np) and
! the i and j indices of the grid point, returns the latitude and longitude of the point

subroutine grid_to_lat_lon(face, lat_ind, lon_ind, del, half_del, lat, lon)

integer,  intent(in)  :: face, lat_ind, lon_ind
real(r8), intent(in)  :: del, half_del
real(r8), intent(out) :: lat, lon

real(r8) :: x, y, blon, blat, rot_angle
real(r8) :: vect(3), rot_vect(3), RZ(3, 3)

! Get the x and y positions on the face
x = sqrt(1.0_r8/3.0_r8) - (half_del + del * (lon_ind - 1))
if(face == 5) then 
   y =  sqrt(1.0_r8/3.0_r8) - (half_del + del * (lat_ind - 1))
else
   y = -sqrt(1.0_r8/3.0_r8) + (half_del + del * (lat_ind - 1))
endif

! These are the faces tangent to the equator
if(face < 4) then
   blon = atan2(sqrt(1.0_r8 / 3.0_r8), x)
   blat = atan2(y, sqrt(1.0_r8/3.0_r8 + x**2.0_r8))
   blon = blon - PI/4.0_r8

   ! Above is for face 0; add PI/2 for each additional face tangent to equator
   lon = blon + PI/2.0_r8 * face; 
   lat = blat;
elseif(face == 4 .or. face == 5) then
   ! Face 4 is tangent to south pole
   lon = atan2(y, x)
   lat = atan2(sqrt(1.0_r8/3.0_r8), sqrt(x**2.0_r8 + y**2.0_r8))

   if(face == 4) lat = -lat

   !  Get ready for rotation
   vect = lat_lon_to_xyz(lat, lon)

   ! Then rotate 45 degrees around Z
   rot_angle = -PI/4.0_r8;

   ! Create the rotation matrix
   RZ(1, 1:3) = [cos(rot_angle),  sin(rot_angle), 0.0_r8]
   RZ(2, 1:3) = [-sin(rot_angle), cos(rot_angle), 0.0_r8]
   RZ(3, 1:3) = [0.0_r8,          0.0_r8,         1.0_r8]
   rot_vect = matmul(RZ, vect)

   lat = asin(rot_vect(3))
   lon = atan2(rot_vect(2), rot_vect(1))
   ! There are inconsistent treatments of the value near longitude
   ! 0 in the grid files for Aether. Some points have a value near or just less
   ! than 2PI, other points have values just greater than 0. This code 
   ! avoids values near to 2PI and converts to near 0 instead.
   if(lon < 0.0_r8) lon =  lon + 2.0_r8*PI
   if(lon >= 2.0_r8*PI) lon = 0.0_r8

endif

end subroutine grid_to_lat_lon

!-----------------------------------------------------------------------

! Convert from latitude longitude, to the face and indices on the face
! of a correpsonding grid point

subroutine lat_lon_to_grid(lat, lon, del, half_del, face, lat_ind, lon_ind)

real(r8), intent(in)  :: lat, lon, del, half_del
integer,  intent(out) :: face, lat_ind, lon_ind

real(r8) :: len(2)

! Get the face and the length along the two imbedded cube faces for the point
call get_face(lat, lon, face, len);

! Figure out which interval this is in along each cube face; This gives 0 to np grid indices
lon_ind = nint((len(1) + half_del) / del)
lat_ind = nint((len(2) + half_del) / del)

end subroutine lat_lon_to_grid

!-----------------------------------------------------------------------

! For points past the edge of a face, finds the corresponding points on the adjacent face

subroutine fix_face(face, lat_grid, lon_grid, np, f_face, f_lat_grid, f_lon_grid, edge, corner)

integer, intent(in)  :: face, lat_grid, lon_grid, np
integer, intent(out) :: f_face, f_lat_grid, f_lon_grid
logical, intent(out) :: edge, corner

integer :: left_neighbor(6),   right_neighbor(6)
integer :: bottom_neighbor(6), top_neighbor(6)
integer :: left_lon_grid(6),   right_lon_grid(6)
integer :: left_lat_grid(6),   right_lat_grid(6)
integer :: bottom_lon_grid(6), top_lon_grid(6)
integer :: bottom_lat_grid(6), top_lat_grid(6)

! Default is not a corner or an edge
corner = .false.
edge   = .false.

! Just return if no edge
if(lon_grid > 0 .and. lon_grid < np + 1 .and. lat_grid > 0 .and. lat_grid < np + 1) then
   f_face = face
   f_lon_grid = lon_grid
   f_lat_grid = lat_grid
   return
endif

! Return illegal value for face information if on a corner
if((lat_grid == 0 .or. lat_grid == np + 1) .and. (lon_grid == 0 .or. lon_grid == np + 1)) then
   corner = .true.
   f_face     = -99
   f_lon_grid = -99
   f_lat_grid = -99
   return
endif

! Otherwise, on an edge
edge = .true.
! Deal with each side of faces separately
if(lon_grid == 0) then
   ! On left edge 
   left_neighbor = [3, 0, 1, 2, 0, 0]
   f_face = left_neighbor(face + 1)
   left_lon_grid = [np, np, np, np, lat_grid, np+1-lat_grid]
   f_lon_grid = left_lon_grid(face + 1)
   left_lat_grid = [lat_grid, lat_grid, lat_grid, lat_grid, 1, np]
   f_lat_grid = left_lat_grid(face + 1)
elseif(lon_grid == np + 1) then
   ! On right edge
   right_neighbor = [1, 2, 3, 0, 2, 2]
   f_face = right_neighbor(face + 1)
   right_lon_grid = [1, 1, 1, 1, np+1-lat_grid, lat_grid]
   f_lon_grid = right_lon_grid(face + 1)
   right_lat_grid = [lat_grid, lat_grid, lat_grid, lat_grid, 1, np]
   f_lat_grid = right_lat_grid(face + 1)
elseif(lat_grid == 0) then
   ! On bottom edge 
   bottom_neighbor = [4, 4, 4, 4, 3, 1]
   f_face = bottom_neighbor(face + 1)
   bottom_lon_grid = [1, lon_grid, np, np+1-lon_grid, np+1-lon_grid, lon_grid]
   f_lon_grid = bottom_lon_grid(face + 1)
   bottom_lat_grid = [lon_grid, np, np+1-lon_grid, 1, 1, np]
   f_lat_grid = bottom_lat_grid(face + 1)
elseif(lat_grid == np + 1) then
   ! On top edge
   top_neighbor = [5, 5, 5, 5, 1, 3]
   f_face = top_neighbor(face + 1)
   top_lon_grid = [1, lon_grid, np, np+1-lon_grid, lon_grid, np+1-lon_grid]
   f_lon_grid = top_lon_grid(face + 1)
   top_lat_grid = [np+1-lon_grid, 1, lon_grid, np, 1, np]
   f_lat_grid = top_lat_grid(face + 1)
endif

end subroutine fix_face

!-----------------------------------------------------------------------

! Returns which face contains (lat, lon_in) and the length from the edge of the point
! along each of the great circle axes.

subroutine get_face(lat, lon_in, face, len)

real(r8), intent(in)  :: lat, lon_in
integer,  intent(out) :: face
real(r8), intent(out) :: len(2)

integer  :: side, rside, rside2
real(r8) :: inv_sqrt_2, rlon, rlon2, gama, gamb, lon
real(r8) :: vec(3), rot_vec(3), rot_vec2(3), rot(3, 3), rot2(3, 3), lon_grid(2), lon_grid_m(2)

! Range adjustment due to Aether roundoff around 2 PI
lon = lon_in
if(lon >= 2.0_r8*PI) lon = 0.0_r8

! Convert lat lon to x y z on unit sphere
vec = lat_lon_to_xyz(lat, lon)

! Get the longitudes for this point in the two rotated spaces

!====================================================================
! Following code shows individual rotations;
! Can collapse these to a single rotation vector for efficiency
! Rotation 90 degrees around y to put pole on equator
!RY = [cosd(90) 0 -sind(90); 0 1 0; sind(90) 0 cosd(90)];
! Rotate 45 degrees around x
!RX = [1 0 0; 0 cosd(45) sind(45); 0 -sind(45) cosd(45)];
! Then rotate 45 degrees around Z
!RZ = [cosd(45) sind(45) 0; -sind(45) cosd(45) 0; 0 0 1];
! Get longitude in the two rotated spaces
!rot_vec = RZ * RX * RY * vec;
!====================================================================
inv_sqrt_2 = 1.0_r8 / sqrt(2.0_r8);
rot(1, 1:3) = [0.5_r8,     0.5_r8,      -inv_sqrt_2]
rot(2, 1:3) = [0.5_r8,     0.5_r8,      inv_sqrt_2 ]
rot(3, 1:3) = [inv_sqrt_2, -inv_sqrt_2, 0.0_r8     ]
rot_vec = matmul(rot, vec)
! Compute the longitude in the rotated space
rlon = atan2(rot_vec(2), rot_vec(1))
if(rlon < 0.0_r8) rlon = rlon + 2.0_r8*PI

!====================================================================
! Can collapse these to a single rotation vector for efficiency
! Rotation 90 degrees around y to put pole on equator
!RY = [cosd(90) 0 -sind(90); 0 1 0; sind(90) 0 cosd(90)];
! Rotate -45 degrees around x
!RX2 = [1 0 0; 0 cosd(-45) sind(-45); 0 -sind(-45) cosd(-45)];
! Then rotate 45 degrees around Z
!RZ = [cosd(45) sind(45) 0; -sind(45) cosd(45) 0; 0 0 1];
! Get longitude in the two rotated spaces
!rot_vec2 = RZ * RX2 * RY * vec;
!====================================================================
rot2(1, 1:3) = [-0.5_r8,     0.5_r8,        -inv_sqrt_2]
rot2(2, 1:3) = [-0.5_r8,     0.5_r8,        inv_sqrt_2 ]
rot2(3, 1:3) = [inv_sqrt_2,  inv_sqrt_2,    0.0_r8     ]
rot_vec2 = matmul(rot2, vec)
! Compute the longitude in the rotated space
rlon2 = atan2(rot_vec2(2), rot_vec2(1))
if(rlon2 < 0.0_r8) rlon2 = rlon2 + 2.0_r8*PI

! Which non-polar side could we be on, 1 to 4
side = floor(lon / (PI/2.0_r8)) + 1.0_r8
! Which rotated 1 side are we on 
rside = floor(rlon / (PI/2.0_r8)) + 1.0_r8
! Which rotated 2 side
rside2 = floor(rlon2 / (PI/2.0_r8)) + 1.0_r8

! Figure out the face from here (0 to 5, 4 is south, 5 is north)
! These are consistent with the numbering on Aether grid files for the cubed sphere
if    ( side == 1 .and. rside  == 1) then
   face = 0; lon_grid(1) = lon;  lon_grid(2) = rlon
elseif( side == 2 .and. rside2 == 1) then
   face = 1; lon_grid(1) = lon;  lon_grid(2) = rlon2
elseif( side == 3 .and. rside  == 3) then 
   face = 2; lon_grid(1) = lon;  lon_grid(2) = rlon
elseif( side == 4 .and. rside2 == 3) then
   face = 3; lon_grid(1) = lon;  lon_grid(2) = rlon2
elseif(rside == 4 .and. rside2 == 4) then
   face = 4; lon_grid(1) = rlon; lon_grid(2) = rlon2
elseif(rside == 2 .and. rside2 == 2) then
   face = 5; lon_grid(1) = rlon; lon_grid(2) = rlon2
endif

! Can use the fact that the projection is equidistant on the imbedded cube to get what fraction 
! across the imbedded rectangle we are
! Take the longitudes and turn them into a number between -sqrt(1/3) and sqrt(1/3)
lon_grid_m = mod(lon_grid, PI/2.0_r8)

! Use law of sines to go from lon back to position along edge of imbedded cube 
! The triangle of interest has a side of length 2sqrt(1/3) (1/2 of the planar diagonal of the imbedded cube)
! The angles adjacent to this side are the longitude and 45 degrees
! The angle opposite the side of length 2sqrt(1/3) is PI - (longitude + PI/4)
! The side opposite the longitude is how far along the side of the cube
! The cube side is 2sqrt(1/3), so the length along the side is between zero and this value
gama   = PI - (PI/4.0_r8 + lon_grid_m(1))
len(1) = sqrt(2.0_r8/3.0_r8) * sin(lon_grid_m(1)) / sin(gama)

gamb   = PI - (PI/4 + lon_grid_m(2))
len(2) = sqrt(2.0_r8/3.0_r8) * sin(lon_grid_m(2)) / sin(gamb)

! If we are on sides 2 or 3, the lengths need to be modified because the grid storage
! for Aether goes from smallest latitude to largest and the longitudes of the shifted
! poles are going the opposite way
if(face == 2 .or. face == 3) len(2) = 2.0_r8 * sqrt(1.0_r8/3.0_r8) - len(2)

! Same for face 4 (the bottom) but it's the other coordinate that's reversed
if(face == 4) len(1) = 2.0_r8*sqrt(1.0_r8/3.0_r8) - len(1)

end subroutine get_face

!-----------------------------------------------------------------------

! Checks to see if the point under consideration is at a corner
! If it is, return the face, lat_index, and lon_index for each of the three bounding points

subroutine get_corners(face, lat_grid, lon_grid, np, lat, lon, del, half_del, &
   f_face, f_lat_grid, f_lon_grid, num_bound_points)
   
integer,  intent(in)  :: face, lat_grid, lon_grid, np
real(r8), intent(in)  :: lat, lon, del, half_del
integer,  intent(out) :: f_face(4), f_lat_grid(4), f_lon_grid(4), num_bound_points

integer  :: corner, quad, i
integer  :: quad_lon_grid(3, 4), quad_lat_grid(3, 4), quad_face(3, 4)
real(r8) :: pxyz(3), qxyz(4, 3), grid_pt_lat, grid_pt_lon

! Default is to find a triangle
num_bound_points = 3

if(face == 0) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 1
   elseif(lat_grid == 0    .and. lon_grid == np+1) then
      corner = 2
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 5
   else
      corner = 6;
   endif
elseif(face == 1) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 2
   elseif(lat_grid == 0    .and. lon_grid == np+1) then 
      corner = 3
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 6
   else
      corner = 7
   endif
elseif(face == 2) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 3
   elseif(lat_grid == 0    .and. lon_grid == np+1) then
      corner = 4
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 7
   else
      corner = 8
   endif
elseif(face == 3) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 4
   elseif(lat_grid == 0    .and. lon_grid == np+1) then
      corner = 1
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 8
   else
      corner = 5;
   endif
elseif(face == 4) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 1
   elseif(lat_grid == 0    .and. lon_grid == np+1) then
      corner = 4
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 2
   else
      corner = 3
   endif
elseif(face == 5) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 6
   elseif(lat_grid == 0    .and. lon_grid == np+1) then
      corner = 7
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 5
   else
      corner = 8;
   endif
endif

! Harvest the information on the grid points bounding the appropriate corner
! Arrays of info for adjacent quads for bulges (three of them, first index)
quad_lon_grid(1:3, 1:4)  = -99
quad_lat_grid(1:3, 1:4)  = -99
quad_face(1:3, 1:4)      = -99

if(corner == 1) then
   f_face(1:3) =     [3,  0, 4]
   f_lon_grid(1:3) = [np, 1, 1]
   f_lat_grid(1:3) = [1,  1, 1]
   quad_face(1, 1:4) = [3,  0, 0, 3]
   quad_face(2, 1:4) = [0,  0, 4, 4]
   quad_face(3, 1:4) = [3,  3, 4, 4]
   quad_lat_grid(1, 1:4) = [1,  1, 2, 2]
   quad_lat_grid(2, 1:4) = [1,  1, 1, 2]
   quad_lat_grid(3, 1:4) = [1,  1, 1, 1]
   quad_lon_grid(1, 1:4) = [np,   1,  1, np]
   quad_lon_grid(2, 1:4) = [1,    2,  1, 1 ]
   quad_lon_grid(3, 1:4) = [np-1, np, 1, 2 ]
elseif(corner == 2) then
   f_face(1:3) =     [0,  1, 4 ]
   f_lon_grid(1:3) = [np, 1, 1 ]
   f_lat_grid(1:3) = [1,  1, np]
   quad_face(1, 1:4) = [0, 1, 1, 0]
   quad_face(2, 1:4) = [1, 1, 4, 4]
   quad_face(3, 1:4) = [0, 0, 4, 4]
   quad_lat_grid(1, 1:4) = [1, 1, 2,  2   ]
   quad_lat_grid(2, 1:4) = [1, 1, np, np  ]
   quad_lat_grid(3, 1:4) = [1, 1, np, np-1]
   quad_lon_grid(1, 1:4) = [np,   1,  1,  np]
   quad_lon_grid(2, 1:4) = [1,    2,  2,  1 ]
   quad_lon_grid(3, 1:4) = [np-1, np, 1,  1 ]
elseif(corner == 3) then
   f_face(1:3) =     [1,  2, 4 ]
   f_lon_grid(1:3) = [np, 1, np]
   f_lat_grid(1:3) = [1,  1, np]
   quad_face(1, 1:4) = [1, 2, 2, 1]
   quad_face(2, 1:4) = [2, 2, 4, 4]
   quad_face(3, 1:4) = [1, 1, 4, 4]
   quad_lat_grid(1, 1:4) = [1, 1, 2,    2 ]
   quad_lat_grid(2, 1:4) = [1, 1, np-1, np]
   quad_lat_grid(3, 1:4) = [1, 1, np,   np]
   quad_lon_grid(1, 1:4) = [np,   1,  1,  np  ]
   quad_lon_grid(2, 1:4) = [1,    2,  np, np  ]
   quad_lon_grid(3, 1:4) = [np-1, np, np, np-1]
elseif(corner == 4) then
   f_face(1:3) =     [2,  3, 4 ]
   f_lon_grid(1:3) = [np, 1, np]
   f_lat_grid(1:3) = [1,  1, 1 ]
   quad_face(1, 1:4) = [2, 3, 3, 2]
   quad_face(2, 1:4) = [3, 3, 4, 4]
   quad_face(3, 1:4) = [2, 2, 4, 4]
   quad_lat_grid(1, 1:4) = [1, 1, 2, 2]
   quad_lat_grid(2, 1:4) = [1, 1, 1, 1]
   quad_lat_grid(3, 1:4) = [1, 1, 1, 2]
   quad_lon_grid(1, 1:4) = [np,   1,  1,    np]
   quad_lon_grid(2, 1:4) = [1,    2,  np-1, np]
   quad_lon_grid(3, 1:4) = [np-1, np, np,   np]
elseif(corner == 5) then
   f_face(1:3) =     [3,  0,  5 ]
   f_lon_grid(1:3) = [np, 1,  1 ]
   f_lat_grid(1:3) = [np, np, np]
   quad_face(1, 1:4) = [3, 0, 0, 3]
   quad_face(2, 1:4) = [0, 0, 5, 5]
   quad_face(3, 1:4) = [3, 3, 5, 5]
   quad_lat_grid(1, 1:4) = [np-1, np-1, np, np  ]
   quad_lat_grid(2, 1:4) = [np,   np,   np, np-1]
   quad_lat_grid(3, 1:4) = [np,   np,   np, np  ]
   quad_lon_grid(1, 1:4) = [np,   1,  1,  np]
   quad_lon_grid(2, 1:4) = [1,    2,  1,  1 ]
   quad_lon_grid(3, 1:4) = [np-1, np, 1,  2 ]
elseif(corner == 6) then
   f_face(1:3) =     [0,  1, 5 ]
   f_lon_grid(1:3) = [np, 1,  1]
   f_lat_grid(1:3) = [np, np, 1]
   quad_face(1, 1:4) = [0, 1, 1, 0]
   quad_face(2, 1:4) = [1, 1, 5, 5]
   quad_face(3, 1:4) = [0, 0, 5, 5]
   quad_lat_grid(1, 1:4) = [np-1, np-1, np, np]
   quad_lat_grid(2, 1:4) = [np,   np,   1,  1 ]
   quad_lat_grid(3, 1:4) = [np,   np,   1,  2 ]
   quad_lon_grid(1, 1:4) = [np,   1,  1,  np]
   quad_lon_grid(2, 1:4) = [1,    2,  2,  1 ]
   quad_lon_grid(3, 1:4) = [np-1, np, 1,  1 ]
elseif(corner == 7) then
   f_face(1:3) =     [1,  2,  5 ]
   f_lon_grid(1:3) = [np, 1,  np]
   f_lat_grid(1:3) = [np, np, 1 ]
   quad_face(1, 1:4) = [1, 2, 2, 1]
   quad_face(2, 1:4) = [2, 2, 5, 5]
   quad_face(3, 1:4) = [1, 1, 5, 5]
   quad_lat_grid(1, 1:4) = [np-1, np-1, np, np]
   quad_lat_grid(2, 1:4) = [np,   np,   1,   2]
   quad_lat_grid(3, 1:4) = [np,   np,   1,   1]
   quad_lon_grid(1, 1:4) = [np,   1,  1,  np  ]
   quad_lon_grid(2, 1:4) = [1,    2,  np, np  ]
   quad_lon_grid(3, 1:4) = [np-1, np, np, np-1]
elseif(corner == 8) then
   f_face(1:3) =     [2,  3,  5 ]
   f_lon_grid(1:3) = [np, 1,  np]
   f_lat_grid(1:3) = [np, np, np]
   quad_face(1, 1:4) = [2, 3, 3, 2]
   quad_face(2, 1:4) = [3, 3, 5, 5]
   quad_face(3, 1:4) = [2, 2, 5, 5];
   quad_lat_grid(1, 1:4) = [np-1, np-1, np, np  ]
   quad_lat_grid(2, 1:4) = [np,   np,   np, np  ]
   quad_lat_grid(3, 1:4) = [np,   np,   np, np-1]
   quad_lon_grid(1, 1:4) = [np,   1,  1,    np]
   quad_lon_grid(2, 1:4) = [1,    2,  np-1, np]
   quad_lon_grid(3, 1:4) = [np-1, np, np,   np]
endif

! Load up the array for the point
pxyz = lat_lon_to_xyz(lat, lon)

! Get lats and lons of the triangle vertices
do i = 1, 3
   call grid_to_lat_lon(f_face(i), f_lat_grid(i), f_lon_grid(i), &
      del, half_del, grid_pt_lat, grid_pt_lon)
   ! Convert to x, y, z coords to check for whether points are in tris/quads
   qxyz(i, 1:3) = lat_lon_to_xyz(grid_pt_lat, grid_pt_lon)
enddo

! See if the point is in the triangle; if so, all is good
if(is_point_in_triangle(qxyz(1, :), qxyz(2, :), qxyz(3, :), pxyz)) return

! If it's not in the triangle, have to check the three adjacent quads at the corner
! This can happen because edges are really on great circles
num_bound_points = 4

do quad = 1, 3
   ! Compute lat and lon for a quad
   do i = 1, 4
         call grid_to_lat_lon(quad_face(quad, i), quad_lat_grid(quad, i), quad_lon_grid(quad, i), &
            del, half_del, grid_pt_lat, grid_pt_lon)
      ! Convert to x, y, z coords to check for whether points are in tris/quads
      qxyz(i, 1:3) = lat_lon_to_xyz(grid_pt_lat, grid_pt_lon)
   enddo

   ! See if the point is inside this quad
   if(is_point_in_quad(qxyz, pxyz)) then
      f_face = quad_face(quad, 1:4)
      f_lat_grid = quad_lat_grid(quad, 1:4)
      f_lon_grid = quad_lon_grid(quad, 1:4)
      return
   endif
enddo

! Falling of the end should not happen;
call error_handler(E_ERR, 'get_corners', 'Reached end of subroutine get_corners', &
   source, 'This should not be possible')

end subroutine get_corners

!-----------------------------------------------------------------------

! Given the latitude and longitude of a point, returns the face, array indices, latitude
! and longitude of the bounding three or four grid points along the number of points
! (3 triangle; 4 quad). np is the number of grid points across each face of the cube sphere.

subroutine get_bounding_box(lat, lon, del, half_del, np, &
   grid_face, grid_lat_ind, grid_lon_ind, grid_pt_lat, grid_pt_lon, num_bound_points)

integer,  intent(in)  :: np
real(r8), intent(in)  :: lat, lon, del, half_del
integer,  intent(out) :: grid_face(4), grid_lat_ind(4), grid_lon_ind(4), num_bound_points
real(r8), intent(out) :: grid_pt_lat(4), grid_pt_lon(4)

integer  :: face, low_grid(2), hi_grid(2), i, my_pt, corner_index
integer  :: lat_grid(4), lon_grid(4), face1_pts(2), face2_pts(2), face1_count, face2_count
real(r8) :: len(2), qxyz(4, 3), pxyz(3)
logical  :: on_edge, edge, corner

! Get the face and the length along the two imbedded cube faces for the point
call get_face(lat, lon, face, len);

! Figure out which interval this is in along each cube face; This gives 0 to np grid indices
low_grid(1) = floor((len(1) + half_del) / del)
low_grid(2) = floor((len(2) + half_del) / del)
hi_grid = low_grid + 1

! Get the indices for the lat and lon directions: Points go counterclockwise starting from lower left 
! For now assume this is a quad, but will correct below if it is a triangle
lat_grid(1) = low_grid(2); lat_grid(2) = hi_grid(2); lat_grid(3) = lat_grid(2); lat_grid(4) = lat_grid(1)
lon_grid(1) = low_grid(1); lon_grid(2) = lon_grid(1); lon_grid(3) = hi_grid(1); lon_grid(4) = lon_grid(3)

! If points are on the edge map to adjacent faces
on_edge = .false.
do i = 1, 4
   call fix_face(face, lat_grid(i), lon_grid(i), np, &
      grid_face(i), grid_lat_ind(i), grid_lon_ind(i), edge, corner)
   ! If any point is on an edge, on_edge is true
   if(edge) on_edge = .true.
   if(corner) then
      corner_index = i
      exit
   endif
enddo

! If it's at a corner, need to find the triangles in a different fashion
! It is possible that the point initially looks like it is in a corner due to the fact 
! that the edges of the grid are on great circles from the corresponding faces, but the
! grid points at the edge are not connected by these great circles.
if(corner) then
   call get_corners(face, lat_grid(corner_index), lon_grid(corner_index), np, &
      lat, lon, del, half_del, grid_face, grid_lat_ind, grid_lon_ind, num_bound_points)
   if(num_bound_points == 4) corner = .false.
else
   ! If not initially at a corner it's definitely in a quad
   num_bound_points = 4
endif

! Compute the lat and lon corresponding to these point
do i = 1, num_bound_points
   call grid_to_lat_lon(grid_face(i), grid_lat_ind(i), grid_lon_ind(i), &
      del, half_del, grid_pt_lat(i), grid_pt_lon(i))
enddo

! Make on_edge true only if we are on an edge but not at a corner
on_edge = (on_edge .and. .not. corner)

if(on_edge) then
   ! If this is an edge, may need to revise box selection
   ! See if the point is in the box (approximately)
   ! Load up the arrays for the vertex points
   do i = 1, num_bound_points   
      ! Convert to x, y, z coords to check for whether points are in tris/quads
      qxyz(i, 1:3) = lat_lon_to_xyz(grid_pt_lat(i), grid_pt_lon(i))
   enddo

   ! Convert point to xyz
   pxyz = lat_lon_to_xyz(lat, lon);

   if(.not. is_point_in_quad(qxyz, pxyz)) then
      ! Not in this box, need to move 'equatorward'
      ! Find indices (from 1 to 4) of points on the same face
      face1_pts(1:2) = 0; face2_pts(1:2) = 0;
      face1_count = 0;    face2_count = 0;
      do i = 1, 4
         if(grid_face(i) == grid_face(1)) then
            face1_count = face1_count + 1;    face1_pts(face1_count) = i
         else
            face2_count = face2_count + 1;    face2_pts(face2_count) = i
         endif
      enddo

      ! First process points of the first face
      ! Are the latitudes or the longitudes on the edge
      if(grid_lon_ind(face1_pts(1)) == grid_lon_ind(face1_pts(2))) then
         ! Adjust the face1 latitudes
         do i = 1, 2
            my_pt = face1_pts(i)
            if(grid_lat_ind(my_pt) > np/2) then
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) - 1
            else
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) + 1
            endif
         enddo
      else
         ! Adjust the face1 longitudes
         do i = 1, 2
            my_pt = face1_pts(i)
            if(grid_lon_ind(my_pt) > np/2) then
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) - 1
            else
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) + 1
            endif
         enddo
      endif

      ! Do the same thing for face2 
      ! Are the latitudes or the longitudes on the edge
      if(grid_lon_ind(face2_pts(1)) == grid_lon_ind(face2_pts(2))) then
         ! Adjust the face2 latitudes
         do i = 1, 2
            my_pt = face2_pts(i);
            if(grid_lat_ind(my_pt) > np/2) then
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) - 1
            else
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) + 1
            endif
         enddo
      else
         ! Adjust the face2 longitudes
         do i = 1, 2
            my_pt = face2_pts(i);
            if(grid_lon_ind(my_pt) > np/2) then
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) - 1
            else
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) + 1
            endif
         enddo
      endif

      ! Compute the lat and lon corresponding to these point
      do i = 1, num_bound_points
         call grid_to_lat_lon(grid_face(i), grid_lat_ind(i), grid_lon_ind(i), &
            del, half_del, grid_pt_lat(i), grid_pt_lon(i))
      enddo
   endif

endif

end subroutine get_bounding_box

!-----------------------------------------------------------------------

! Given the index of a horizontal column, returns the latitude and longitude in radians

subroutine col_index_to_lat_lon(col_index, np, del, half_del, lat, lon)

integer,  intent(in)  :: col_index, np
real(r8), intent(in)  :: del, half_del
real(r8), intent(out) :: lat, lon

integer :: face, resid, lat_ind, lon_ind

! Which face are we on? np**2 points per face
face = (col_index - 1) / (np**2)
resid = col_index - face * np**2

! Get latitude index
lat_ind = (resid - 1) / np + 1

lon_ind = resid - (lat_ind - 1) * np

! Get the corresponding latitude and longitude
call grid_to_lat_lon(face, lat_ind, lon_ind, del, half_del, lat, lon)

end subroutine col_index_to_lat_lon

!-----------------------------------------------------------------------

! Given a latitude and longitude, return corresponding model column index

function lat_lon_to_col_index(lat, lon, del, half_del, np)

integer              :: lat_lon_to_col_index
real(r8), intent(in) :: lat, lon, del, half_del
integer,  intent(in) :: np 

integer :: face, lat_ind, lon_ind

! Get the face, lat_ind and lon_ind
call lat_lon_to_grid(lat, lon, del, half_del, face, lat_ind, lon_ind)

! Confirm that this duplicates the version that requires initialization;
lat_lon_to_col_index = face * np * np + (lat_ind - 1) * np + lon_ind

end function lat_lon_to_col_index

!-----------------------------------------------------------------------

! Convert latitude and longitude to 3D x,y,z

function lat_lon_to_xyz(lat, lon)          

real(r8)             :: lat_lon_to_xyz(3)
real(r8), intent(in) :: lat, lon

lat_lon_to_xyz(1) = cos(lat) * cos(lon)
lat_lon_to_xyz(2) = cos(lat) * sin(lon)
lat_lon_to_xyz(3) = sin(lat)

end function lat_lon_to_xyz

!-----------------------------------------------------------------------

! Determines if the projection of a point p onto the plane of a triangle with vertices
! v1, v2 and v3 is inside the triangle or not. Computes the areas of each of the triangles
! between p and a pair of vertices. These should sum to the area of the triangle if p
! is inside and be larger than that if p is outside. 

function is_point_in_triangle(v1, v2, v3, p)

logical              :: is_point_in_triangle
real(r8), intent(in) :: v1(3), v2(3), v3(3), p(3)

real(r8) :: a(3), b(3), perp(3), unit_perp(3), p_proj(3)
real(r8) :: offset, len_s1, len_s2, len_s3, len_p1, len_p2, len_p3, at, at1, at2, at3
real(r8) :: area_dif, dif_frac, threshold

! Get the projection of the point p onto the plane containing the triangle
! Start by getting perpendicular vector to plane by cross product
a = v1 - v2
b = v2 - v3
perp(1) = a(2) * b(3) - a(3) * b(2)
perp(2) = a(3) * b(1) - a(1) * b(3)
perp(3) = a(1) * b(2) - a(2) * b(1)
! Get unit vector in direction of perp
unit_perp = perp / sqrt(dot_product(perp, perp))
! Projection of vector from v1 to p on the unit perp vector is how much to move to get to plane
offset = dot_product((p-v1), unit_perp)
p_proj = p - offset * unit_perp

! Compute lengths of the sides
len_s1 = sqrt(dot_product(v1-v2, v1-v2))
len_s2 = sqrt(dot_product(v3-v2, v3-v2))
len_s3 = sqrt(dot_product(v1-v3, v1-v3))

! Compute the lengths from the point p
len_p1 = sqrt(dot_product(p_proj-v1, p_proj-v1))
len_p2 = sqrt(dot_product(p_proj-v2, p_proj-v2))
len_p3 = sqrt(dot_product(p_proj-v3, p_proj-v3))

! Area of triangle
at = heron(len_s1, len_s2, len_s3)

! Compute areas of sub triangles
at1 = heron(len_p1, len_p2, len_s1)
at2 = heron(len_p2, len_p3, len_s2)
at3 = heron(len_p3, len_p1, len_s3)

! Difference between sub triangles and the triangle area
area_dif = at1 + at2 + at3 - at

! Quadrilaterals on the interior of the cube sphere sides are really spherical quads,
! Their sides are great circles. This routine assumes that the triangles composing the quads
! have straight sides in regular space. The algorithm finds points that are inside the 
! spherical quads. These quads actually 'bulge' out compared to the regular sides, so it is possible
! to have points that are inside the spherical quad but just barely outside of the regular
! quads. This threshold is tuned so that these points still show as inside. The tuning is for
! np = 18 (number of points along a grid face is 18). Fewer points might require a larger
! threshold while more points might be okay with a smaller one.
threshold = 0.002_r8

dif_frac = area_dif / at
is_point_in_triangle = abs(dif_frac) < threshold

end function is_point_in_triangle

!-----------------------------------------------------------------------

! Returns true if point p is in quadrilateral with vertices v
   
function is_point_in_quad(v, p)

logical              :: is_point_in_quad
real(r8), intent(in) :: v(4, 3), p(3)

logical :: inside_t(4)

integer :: i
   
! See if the point is inside this quad; it's inside if it's in one or more contained triangles
inside_t(1) = is_point_in_triangle(v(1, :), v(2, :), v(3, :), p)
inside_t(2) = is_point_in_triangle(v(1, :), v(2, :), v(4, :), p)
inside_t(3) = is_point_in_triangle(v(1, :), v(3, :), v(4, :), p)
inside_t(4) = is_point_in_triangle(v(2, :), v(3, :), v(4, :), p)

is_point_in_quad = any(inside_t)
   
end function is_point_in_quad

!-----------------------------------------------------------------------

! Computes Herons formula to get area of triangle from lenghts of sides
! Super accuracy is not needed in the area calculation here

function heron(a, b, c)

real(r8)             :: heron
real(r8), intent(in) :: a, b, c

real(r8) :: s, arg
                                 
s = (a + b + c) /2
arg = (s * (s - a) * (s - b) * (s - c))

! Make sure we don't roundoff to a negative
if(arg <= 0.0_r8) then
   heron = 0.0_r8
else
   heron = sqrt(arg)
endif

end function heron

!-----------------------------------------------------------------------

end module cube_sphere_grid_tools_mod
