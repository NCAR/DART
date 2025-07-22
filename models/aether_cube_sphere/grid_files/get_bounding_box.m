% Given the latitude and longitude of a point, returns the face, array indices, latitude
% and longitude of the bounding three or four grid points along with a flag that indicates
% if the bounding polygon is a triangle (grid corner) or quadrilateral. np is the number
% of grid points across each face.
function [grid_face, grid_lat_ind, grid_lon_ind, grid_pt_lat, grid_pt_lon, on_edge, corner] = ...
   get_bounding_box(lat, lon, np)

% Some geometry about the grid point distribution
% Cube side is divided into np-1 interior intervals of width 2sqrt(1/3) / np and
% two exterior intervals of half  width, sqrt(1/3) / np 
cube_side = 2 * sqrt(1 / 3); 
del = cube_side / np; 
half_del = del / 2;

% Get the face and the length along the two imbedded cube faces for the point
[face, len] = get_face(lat, lon);

% Figure out which interval this is in along each cube face; This gives 0 to np intervals
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
if(corner)
   [grid_face, grid_lat_ind, grid_lon_ind] = get_corners(face, lat_grid(corner_index), lon_grid(corner_index), np);
   num_bound_points = 3;
else
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
      qxyz(i, 1) = cos(grid_pt_lat(i)) .* cos(grid_pt_lon(i));
      qxyz(i, 2) = cos(grid_pt_lat(i)) .* sin(grid_pt_lon(i));
      qxyz(i, 3) = sin(grid_pt_lat(i));
   end

   % Load up the array for the point
   pxyz(1) = cos(lat) .* cos(lon);
   pxyz(2) = cos(lat) .* sin(lon);
   pxyz(3) = sin(lat);
   inside_t(1:4) = false;
   % See if the point is inside a quad; it's inside if it's in one or more contained triangles
   [inside_t(1), dif_frac_t(1)] = is_point_in_triangle(qxyz(1, :), qxyz(2, :), qxyz(3, :), pxyz);
   [inside_t(2), dif_frac_t(2)] = is_point_in_triangle(qxyz(1, :), qxyz(2, :), qxyz(4, :), pxyz);
   [inside_t(3), dif_frac_t(3)] = is_point_in_triangle(qxyz(1, :), qxyz(3, :), qxyz(4, :), pxyz);
   [inside_t(4), dif_frac_t(4)] = is_point_in_triangle(qxyz(2, :), qxyz(3, :), qxyz(4, :), pxyz);

   % Do our own temporary version of whether or not we are in the quad
   if(min(dif_frac_t) > 1e-3)
      change_lat = false;
      change_lon = false;
      % Not really in this box, need to move 'equatorward'
      dif_frac_t
      grid_face
      grid_lat_ind
      grid_lon_ind
      % Figure out if it is the latitude or longitude indices that need to be shifted
      % The ones that don't need to be shifted will be on the edges (indices 1 and np)
      % See if the longitude grid is the one that is already on great circles
      is_one = grid_lon_ind == 1;
      is_np  = grid_lon_ind == 18;
      if(sum(is_one) + sum(is_np) == 4) 
         change_lat = true; 
         fprintf('lat needs changed\n');
         %change the latitudes to move towards the 'equator'
         for i = 1:4
            if(grid_lat_ind(i) >= np/2) 
               grid_lat_ind(i) = grid_lat_ind(i) - 1;
            else
               grid_lat_ind(i) = grid_lat_ind(i) + 1;
            end
         end
       
      else
         change_lon = true;
         fprintf('lon needs changed\n');
         % Change the longitudes to move towards the 'equator'
         for i = 1:4
            if(grid_lon_ind(i) >= np/2) 
               grid_lon_ind(i) = grid_lon_ind(i) - 1;
            else
               grid_lon_ind(i) = grid_lon_ind(i) + 1;
            end
         end
         
      end
      grid_lat_ind
      grid_lon_ind
      % Compute the lat and lon corresponding to these point
      for i = 1:num_bound_points
         [grid_pt_lat(i), grid_pt_lon(i)] = grid_to_lat_lon(grid_face(i), grid_lat_ind(i), grid_lon_ind(i), np);
      end
   end

end









