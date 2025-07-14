function [f_face, f_lon_grid, f_lat_grid] = fix_face(face, lon_grid, lat_grid, np)

% For points past the edge of a face, finds the corresponding points on the adjacent face
% Need to do something more special for corner points

[face lon_grid lat_grid np]

% Just return if no edge
if(lon_grid > 0 && lon_grid < np + 1 && lat_grid > 0 && lat_grid < np + 1)
   f_face = face;
   f_lon_grid = lon_grid;
   f_lat_grid = lat_grid;
   return
end

if((lon_grid == 0 | lon_grid == np + 1) && (lat_grid == 0 | lat_grid == np + 1))
   % These are corner points
else
   if(lon_grid == 0)
      % On left edge not corner
      left_neighbor = [3 0 1 2 0 0];
      f_face = left_neighbor(face + 1);
      left_lon_grid = [np np np np lat_grid np+1-lat_grid];
      f_lon_grid = left_lon_grid(face + 1);
      left_lat_grid = [lat_grid lat_grid lat_grid lat_grid 1 np];
      f_lat_grid = left_lat_grid(face + 1);
   elseif(lon_grid == np + 1)
      % On right edge not corner
      right_neighbor = [1 2 3 0 2 2];
      f_face = right_neighbor(face + 1);
      right_lon_grid = [1 1 1 1 np+1-lat_grid lat_grid];
      f_lon_grid = right_lon_grid(face + 1);
      right_lat_grid = [lat_grid lat_grid lat_grid lat_grid 1 np];
      f_lat_grid = right_lat_grid(face + 1);
      [f_face f_lon_grid f_lat_grid]
   elseif(lat_grid == 0)
      % On bottom edge not corner
      bottom_neighbor = [4 4 4 4 3 1];
      f_face = bottom_neighbor(face + 1);
      bottom_lon_grid = [1 lon_grid np np+1-lon_grid np+1-lon_grid lon_grid];
      f_lon_grid = bottom_lon_grid(face + 1);
      bottom_lat_grid = [lon_grid np np+1-lon_grid 1 1 np];
      f_lat_grid = bottom_lat_grid(face + 1);
   elseif(lat_grid == np + 1)
      % On top edge not corner
      top_neighbor = [5 5 5 5 1 3];
      f_face = top_neighbor(face + 1);
      top_lon_grid = [1 lon_grid np np+1-lon_grid lon_grid np+1-lon_grid];
      f_lon_grid = top_lon_grid(face + 1);
      top_lat_grid = [np+1-lon_grid 1 lon_grid np 1 np];
      f_lat_grid = top_lat_grid(face + 1);
   end
   
end
