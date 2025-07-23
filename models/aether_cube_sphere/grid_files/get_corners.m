function [f_face, f_lat_grid, f_lon_grid, num_bound_points] = get_corners(face, lat_grid, lon_grid, lat, lon, np)

% Checks to see if the point under consideration is at a corner
% If it is, return the face, lat_index, and lon_index for each of the three bounding points

% Default is to find a triangle
num_bound_points = 3;

if(face == 0)
   if    (lat_grid == 0    && lon_grid == 0   ) corner = 1;
   elseif(lat_grid == 0    && lon_grid == np+1) corner = 2;
   elseif(lat_grid == np+1 && lon_grid == 0   ) corner = 5;
   else                                         corner = 6;
   end
elseif(face == 1)
   if    (lat_grid == 0    && lon_grid == 0   ) corner = 2;
   elseif(lat_grid == 0    && lon_grid == np+1) corner = 3;
   elseif(lat_grid == np+1 && lon_grid == 0   ) corner = 6;
   else                                         corner = 7;
   end
elseif(face == 2)
   if    (lat_grid == 0    && lon_grid == 0   ) corner = 3;
   elseif(lat_grid == 0    && lon_grid == np+1) corner = 4;
   elseif(lat_grid == np+1 && lon_grid == 0   ) corner = 7;
   else                                         corner = 8;
   end
elseif(face == 3)
   if    (lat_grid == 0    && lon_grid == 0   ) corner = 4;
   elseif(lat_grid == 0    && lon_grid == np+1) corner = 1;
   elseif(lat_grid == np+1 && lon_grid == 0   ) corner = 8;
   else                                         corner = 5;
   end
elseif(face == 4)
   if    (lat_grid == 0    && lon_grid == 0   ) corner = 1;
   elseif(lat_grid == 0    && lon_grid == np+1) corner = 4;
   elseif(lat_grid == np+1 && lon_grid == 0   ) corner = 2;
   else                                         corner = 3;
   end
elseif(face == 5)
   if    (lat_grid == 0    && lon_grid == 0   ) corner = 6;
   elseif(lat_grid == 0    && lon_grid == np+1) corner = 7;
   elseif(lat_grid == np+1 && lon_grid == 0   ) corner = 5;
   else                                         corner = 8;
   end
end

% Harvest the information on the grid points bounding the appropriate corner
% Arrays of info for adjacent quads for bulges (three of them, first index)
quad_lon_grid(1:3, 1:4) = 0;  quad_lat_grid(1:3, 1:4) = 0; quad_face_grid(1:3, 1:4) = -99;

if(corner == 1)
   f_face(1:3) = [3 0 4];
   f_lon_grid(1:3) = [np 1 1];
   f_lat_grid(1:3) = [1  1 1];
   quad_face =     [3  0 0 3 ;    0  0 4 4;     3    3  4 4];
   quad_lat_grid = [1  1 2 2 ;    1  1 1 2;     1    1  1 1];
   quad_lon_grid = [np 1 1 np;    1  2 1 1;     np-1 np 1 2];
elseif(corner == 2)
   f_face(1:3) = [0 1 4];
   f_lon_grid(1:3) = [np 1 1 ];
   f_lat_grid(1:3) = [1  1 np];
   quad_face =     [0  1 1 0 ;    1 1 4  4 ;      0    0  4  4   ];
   quad_lat_grid = [1  1 2 2 ;    1 1 np np;      1    1  np np-1];
   quad_lon_grid = [np 1 1 np;    1 2 2  1 ;      np-1 np 1  1   ];
elseif(corner == 3)
   f_face(1:3) = [1 2 4];
   f_lon_grid(1:3) = [np 1 np];
   f_lat_grid(1:3) = [1  1 np];
   quad_face =     [1  2 2 1 ;  2 2 4    4 ;      1    1  4  4   ];
   quad_lat_grid = [1  1 2 2 ;  1 1 np-1 np;      1    1  np np  ];
   quad_lon_grid = [np 1 1 np;  1 2 np   np;      np-1 np np np-1];
elseif(corner == 4)
   f_face(1:3) = [2 3 4];
   f_lon_grid(1:3) = [np 1 np];
   f_lat_grid(1:3) = [1  1 1 ];
   quad_face = [2 3 3 2;         3  3 4    4 ;     2    2  4  4 ];
   quad_lat_grid = [1 1 2 2;     1  1 1    1 ;     1    1  1  2 ];
   quad_lon_grid = [np 1 1 np;   1  2 np-1 np;     np-1 np np np];
elseif(corner == 5)
   f_face(1:3) = [3 0 5];
   f_lon_grid(1:3) = [np 1  1 ];
   f_lat_grid(1:3) = [np np np];
   quad_face =     [3    0    0  3 ;   0  0  5  5   ;  3    3  5  5 ];
   quad_lat_grid = [np-1 np-1 np np;   np np np np-1;  np   np np np];
   quad_lon_grid = [np 1 1    np   ;   1  2  1  1   ;  np-1 np 1  2 ];
elseif(corner == 6)
   f_face(1:3) = [0 1 5];
   f_lon_grid(1:3) = [np 1  1];
   f_lat_grid(1:3) = [np np 1];
   quad_face =     [0    1    1  0 ;   1  1  5 5;    0    0  5  5];
   quad_lat_grid = [np-1 np-1 np np;   np np 1 1;    np   np 1  2];
   quad_lon_grid = [np   1    1  np;   1  2  2 1;    np-1 np 1  1];
elseif(corner == 7)
   f_face(1:3) = [1 2 5];
   f_lon_grid(1:3) = [np 1  np];
   f_lat_grid(1:3) = [np np 1 ];
   quad_face =     [1    2    2  1 ;   2  2  5   5 ;     1    1  5  5   ];
   quad_lat_grid = [np-1 np-1 np np;   np np 1   2 ;     np   np 1  1   ];
   quad_lon_grid = [np   1    1  np;   1  2  np  np;     np-1 np np np-1];
elseif(corner == 8)
   f_face(1:3) = [2 3 5];
   f_lon_grid(1:3) = [np 1  np];
   f_lat_grid(1:3) = [np np np];
   quad_face =     [2    3    3  2 ;    3  3  5    5 ;    2    2   5  5   ];
   quad_lat_grid = [np-1 np-1 np np;    np np np   np;    np   np  np np-1];
   quad_lon_grid = [np   1    1  np;    1  2  np-1 np;    np-1 np  np np  ];
end

% See if the point is in the triangle; if so, all is good
% Load up the array for the point
pxyz(1) = cos(lat) .* cos(lon);
pxyz(2) = cos(lat) .* sin(lon);
pxyz(3) = sin(lat);

% Get lats and lons of the triangle vertices
% Compute the lat and lon corresponding to these point
for i = 1:3
   [grid_pt_lat(i), grid_pt_lon(i)] = grid_to_lat_lon(f_face(i), f_lat_grid(i), f_lon_grid(i), np);
   % Convert to x, y, z coords to check for whether points are in tris/quads
   qxyz(i, 1) = cos(grid_pt_lat(i)) .* cos(grid_pt_lon(i));
   qxyz(i, 2) = cos(grid_pt_lat(i)) .* sin(grid_pt_lon(i));
   qxyz(i, 3) = sin(grid_pt_lat(i));
end

[inside, dif_frac] = is_point_in_triangle(qxyz(1, :), qxyz(2, :), qxyz(3, :), pxyz);
if(inside) return; end

% If it's not in the triangle, have to check the adjacent quads
num_bound_points = 4;

for quad = 1:3
   % Compute lat and lon for a quad
   for i = 1:4
      [grid_pt_lat(i), grid_pt_lon(i)] = ...
         grid_to_lat_lon(quad_face(quad, i), quad_lat_grid(quad, i), quad_lon_grid(quad, i), np);
      % Convert to x, y, z coords to check for whether points are in tris/quads
      qxyz(i, 1) = cos(grid_pt_lat(i)) .* cos(grid_pt_lon(i));
      qxyz(i, 2) = cos(grid_pt_lat(i)) .* sin(grid_pt_lon(i));
      qxyz(i, 3) = sin(grid_pt_lat(i));
   end

   % See if the point is inside this quad; it's inside if it's in one or more contained triangles
   [inside_t(1), dif_frac_t(1)] = is_point_in_triangle(qxyz(1, :), qxyz(2, :), qxyz(3, :), pxyz);
   [inside_t(2), dif_frac_t(2)] = is_point_in_triangle(qxyz(1, :), qxyz(2, :), qxyz(4, :), pxyz);
   [inside_t(3), dif_frac_t(3)] = is_point_in_triangle(qxyz(1, :), qxyz(3, :), qxyz(4, :), pxyz);
   [inside_t(4), dif_frac_t(4)] = is_point_in_triangle(qxyz(2, :), qxyz(3, :), qxyz(4, :), pxyz);

   if(any(inside_t))
      f_face = quad_face(quad, :);
      f_lat_grid = quad_lat_grid(quad, :);
      f_lon_grid = quad_lon_grid(quad, :);
      return
   end
end

% Falling of the end is not happy
fprintf('UNEXPECTED FAILURE IN GET_CORNERS.M\n');
stop





