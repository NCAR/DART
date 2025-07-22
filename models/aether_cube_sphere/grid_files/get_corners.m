function [f_face, f_lat_grid, f_lon_grid] = get_corners(face, lat_grid, lon_grid, np)

[face lat_grid lon_grid]
% Checks to see if the point under consideration is at a corner
% If it is, return the face, lat_index, and lon_index for each of the three bounding points

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
if(corner == 1)
   f_face(1:3) = [3 0 4];
   f_lon_grid(1:3) = [np 1 1];
   f_lat_grid(1:3) = [1  1 1];
elseif(corner == 2)
   f_face(1:3) = [0 1 4];
   f_lon_grid(1:3) = [np 1 1 ];
   f_lat_grid(1:3) = [1  1 np];
elseif(corner == 3)
   f_face(1:3) = [1 2 4];
   f_lon_grid(1:3) = [np 1 np];
   f_lat_grid(1:3) = [1  1 np];
elseif(corner == 4)
   f_face(1:3) = [2 3 4];
   f_lon_grid(1:3) = [np 1 np];
   f_lat_grid(1:3) = [1  1 1 ];
elseif(corner == 5)
   f_face(1:3) = [3 0 5];
   f_lon_grid(1:3) = [np 1  1 ];
   f_lat_grid(1:3) = [np np np];
elseif(corner == 6)
   f_face(1:3) = [0 1 5];
   f_lon_grid(1:3) = [np 1  1];
   f_lat_grid(1:3) = [np np 1];
elseif(corner == 7)
   f_face(1:3) = [1 2 5];
   f_lon_grid(1:3) = [np 1  np];
   f_lat_grid(1:3) = [np np 1 ];
elseif(corner == 8)
   f_face(1:3) = [2 3 5];
   f_lon_grid(1:3) = [np 1  np];
   f_lat_grid(1:3) = [np np np];
end

