% Given the latitude and longitude of a point, returns the face, array indices, latitude
% and longitude of the bounding three or four grid points along the number of points
% (3 triangle; 4 quad). np is the number of grid points across each face of the cube sphere.

function [grid_face, grid_lat_ind, grid_lon_ind, grid_pt_lat, grid_pt_lon, num_bound_points] = ...
   get_bounding_box(lat, lon, np)

% Cube side is divided into np-1 interior intervals of width 2sqrt(1/3) / np and
% two exterior intervals with half the width, sqrt(1/3) / np 
cube_side = 2 * sqrt(1 / 3); 
del = cube_side / np; 
half_del = del / 2;

% Get the face and the length along the two imbedded cube faces for the point
[face, len] = get_face(lat, lon);

% Figure out which interval this is in along each cube face; This gives 0 to np grid indices
low_grid = floor((len + half_del) / del);
hi_grid = low_grid + 1;

% Get the indices for the lat and lon directions: Points go counterclockwise starting from lower left 
% For now assume this is a quad, but will correct below if it is a triangle
lat_grid(1) = low_grid(2); lat_grid(2) = hi_grid(2); lat_grid(3) = lat_grid(2); lat_grid(4) = lat_grid(1);
lon_grid(1) = low_grid(1); lon_grid(2) = lon_grid(1); lon_grid(3) = hi_grid(1); lon_grid(4) = lon_grid(3);

% If points are on the edge map to adjacent faces
on_edge = false;
grid_lon_ind(1:4) = 0; grid_lat_ind(1:4) = 0;
for i = 1:4
   [grid_face(i), grid_lat_ind(i), grid_lon_ind(i), edge, corner] = fix_face(face, lat_grid(i), lon_grid(i), np);
   % If any point is on an edge, on_edge is true
   if(edge) on_edge = true; end
   if(corner)
      corner_index = i;
      break
   end
end

% If it's at a corner, need to find the triangles in a different fashion
% It is possible that the point initially looks like it is in a corner due to the fact 
% that the edges of the grid are on great circles from the corresponding faces, but the
% grid point at the edge are not connected by these great circles.
if(corner)
   [grid_face, grid_lat_ind, grid_lon_ind, num_bound_points] = ...
      get_corners(face, lat_grid(corner_index), lon_grid(corner_index), lat, lon, np);
   if(num_bound_points == 4) corner = false; end
else
   % If not initially at a corner it's definitely in a quad
   num_bound_points = 4;
end

% Compute the lat and lon corresponding to these point
for i = 1:num_bound_points
   [grid_pt_lat(i), grid_pt_lon(i)] = grid_to_lat_lon(grid_face(i), grid_lat_ind(i), grid_lon_ind(i), np);
end

% Make on_edge true only if we are on an edge but not at a corner
on_edge = on_edge && ~corner;

if(on_edge)
   % If this is an edge, may need to revise box selection
   % See if the point is in the box (approximately)
   % Load up the arrays for the vertex points
   for i = 1:num_bound_points   
      % Convert to x, y, z coords to check for whether points are in tris/quads
      qxyz(i, 1:3) = lat_lon_to_xyz(grid_pt_lat(i), grid_pt_lon(i));
   end

   % Convert point to xyz
   pxyz = lat_lon_to_xyz(lat, lon);
   inside_t(1:4) = false;
   % See if the point is inside a quad; it's inside if it's in one or more contained triangles
   [inside_t(1), dif_frac_t(1)] = is_point_in_triangle(qxyz(1, :), qxyz(2, :), qxyz(3, :), pxyz);
   [inside_t(2), dif_frac_t(2)] = is_point_in_triangle(qxyz(1, :), qxyz(2, :), qxyz(4, :), pxyz);
   [inside_t(3), dif_frac_t(3)] = is_point_in_triangle(qxyz(1, :), qxyz(3, :), qxyz(4, :), pxyz);
   [inside_t(4), dif_frac_t(4)] = is_point_in_triangle(qxyz(2, :), qxyz(3, :), qxyz(4, :), pxyz);

   if(~any(inside_t))
      % Not in this box, need to move 'equatorward'
      % Find indices (from 1 to 4) of points on the same face
      face1_pts(1:2) = 0; face2_pts(1:2) = 0;
      face1_count = 0;    face2_count = 0;
      for i = 1:4
         if(grid_face(i) == grid_face(1)) 
            face1_count = face1_count + 1;    face1_pts(face1_count) = i;
         else
            face2_count = face2_count + 1;    face2_pts(face2_count) = i;
         end
      end

      % First process points of the first face
      % Are the latitudes or the longitudes on the edge
      if(grid_lon_ind(face1_pts(1)) == grid_lon_ind(face1_pts(2)))
         % Adjust the face1 latitudes
         for i = 1:2
            my_pt = face1_pts(i);
            if(grid_lat_ind(my_pt) > np/2) 
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) - 1;
            else
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) + 1;
            end
         end
      else
         % Adjust the face1 longitudes
         for i = 1:2
            my_pt = face1_pts(i);
            if(grid_lon_ind(my_pt) > np/2) 
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) - 1;
            else
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) + 1;
            end
         end
      end

      % Do the same thing for face2 
      % Are the latitudes or the longitudes on the edge
      if(grid_lon_ind(face2_pts(1)) == grid_lon_ind(face2_pts(2)))
         % Adjust the face1 latitudes
         for i = 1:2
            my_pt = face2_pts(i);
            if(grid_lat_ind(my_pt) > np/2) 
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) - 1;
            else
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) + 1;
            end
         end
      else
         % Adjust the face1 longitudes
         for i = 1:2
            my_pt = face2_pts(i);
            if(grid_lon_ind(my_pt) > np/2) 
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) - 1;
            else
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) + 1;
            end
         end
      end

      % Compute the lat and lon corresponding to these point
      for i = 1:num_bound_points
         [grid_pt_lat(i), grid_pt_lon(i)] = grid_to_lat_lon(grid_face(i), grid_lat_ind(i), grid_lon_ind(i), np);
      end
   end

end

