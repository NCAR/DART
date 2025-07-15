% Initial work on finding the actual 3 or 4 grid points that create triangle or quad that
% contain a given lat/lon on the Aether cubed sphere.

% Parameter that has the number of model grid points along each dimension 
% This does not include halos; the points are offset from the boundaries
np = 18;

% Some geometry about the grid point distribution
% Cube side is divided into np-1 interior intervals of width 2sqrt(1/3) / np and
% two exterior intervals of half  width, sqrt(1/3) / np 
cube_side = 2 * sqrt(1 / 3); 
del = cube_side / np;
half_del = del / 2;

% Precompute the positions of the points along the cube
% These are only used to verify that these grid computations give the same values as the grid files
for i = 1:18
   x(i) = -sqrt(1/3) + half_del + del * (i - 1);
   box_lon(i) = atan2(sqrt(1 / 3), x(i));
end
box_lon = box_lon - pi / 4;

%-------------------------------------------
% Plot a background spherical surface
th = linspace(0,2*pi, 1000) ;
phi = linspace(-pi/2, pi/2, 1000) ;
[T,P] = meshgrid(th,phi) ;

X = cos(P) .* cos(T);
Y = cos(P) .* sin(T);
Z = sin(P);

% Make the grey scale vary by face
for i = 1:size(X, 1)
   for j = 1:size(X, 2)
      lat(i, j) = asin(Z(i, j));
      lon(i, j) = atan2(Y(i, j), X(i, j));
      if(lon(i, j) < 0) lon(i, j) =  lon(i, j) + 2*pi; end
      if(lon(i, j) == 2*pi) lon(i, j) = 0; end
      [face, len] = get_face(lat(i, j), lon(i, j));
      C(i, j) = face;
   end
end

h = surf(X,Y,Z, C);
set(h, 'linestyle', 'none')
axis equal
% Change the color map
for i = 1:size(colormap, 1)
   c_grey(i, 1:3) = 0.6 + (i / size(colormap, 1) * 0.2);
end
colormap(c_grey);
colorbar
hold on
%-------------------------------------------

% Open grid files for each face for comparison
for face = 0:5
   % Generate the grid file name
   fname = strcat('grid_g000', num2str(face), '.nc');
   xt = ncread(fname, 'Longitude');
   glon(face + 1, :, :) = squeeze(xt(1, 3:end-2, 3:end-2));
   yt = ncread(fname, 'Latitude');
   glat(face + 1, :, :) = squeeze(yt(1, 3:end-2, 3:end-2));

   % Get x, y, z for this face for plotting
   x = cos(glat(face + 1, :, :)) .* cos(glon(face + 1, :, :));
   y = cos(glat(face + 1, :, :)) .* sin(glon(face + 1, :, :));
   z = sin(glat(face + 1, :, :));
   % Plot points on each face with distinct colors 
   h = plot3(squeeze(x), squeeze(y), squeeze(z), '.', 'markersize', 24, 'color', 'blue', 'linewidth', 16);
   hold on

end
%-------------------------------------------

% Enter lats and lons in degrees for starters
pt_lon_d = 269.7;
pt_lat_d = -35;
% Convert to radians since this is generally used in algorithm
pt_lon = deg2rad(pt_lon_d);
pt_lat = deg2rad(pt_lat_d);
% Convert to x, y, z for plotting
px = cos(pt_lat) .* cos(pt_lon);
py = cos(pt_lat) .* sin(pt_lon);
pz = sin(pt_lat);
% Plot the point in green
plot3(px, py, pz, '.', 'markersize', 24, 'color', 'green', 'linewidth', 16);

% Get the face, the longitudes on the two projections that don't have a pole in this face
% and the length along the two imbedded cube faces.
[face, len] = get_face(pt_lat, pt_lon); 

% Each of the four bounding points gets this face initially
% Points go counterclockwise starting from lower left 
grid_face(1:4) = face;

% Figure out which interval this is in along each cube face; This gives 0 to np intervals
low_grid = floor((len + half_del) / del);
hi_grid = low_grid + 1;

% The longitude grid values are
lat_grid(1) = low_grid(2); lat_grid(2) = hi_grid(2); lat_grid(3) = lat_grid(2); lat_grid(4) = lat_grid(1);
lon_grid(1) = low_grid(1); lon_grid(2) = lon_grid(1); lon_grid(3) = hi_grid(1); lon_grid(4) = lon_grid(3);

% If points are on the edge map to adjacent faces
f_lon_grid(1:4) = 0; f_lat_grid(1:4) = 0;
for i = 1:4
   [f_face(i), f_lat_grid(i), f_lon_grid(i), corner_detected] = fix_face(face, lat_grid(i), lon_grid(i), np);
   if(corner_detected) 
      corner_index = i;
      break
   end
end

% At this point, either plot for corner triangles or for quads
if(corner_detected) 
   num_bound_points = 3;
   [f_face, f_lat_grid, f_lon_grid] = get_corners(face, lat_grid(corner_index), lon_grid(corner_index), np);
else 
   num_bound_points = 4;
end

% Highlight the bounding points;
for i = 1:num_bound_points
   bp_lon(i) = glon(f_face(i) + 1, f_lat_grid(i), f_lon_grid(i));
   bp_lat(i) = glat(f_face(i) + 1, f_lat_grid(i), f_lon_grid(i));

   % Convert each to x, y, z and plot
   bpx = cos(bp_lat(i)) .* cos(bp_lon(i));
   bpy = cos(bp_lat(i)) .* sin(bp_lon(i));
   bpz = sin(bp_lat(i));
   % Plot the bounding points in red
   plot3(bpx, bpy, bpz, '.', 'markersize', 24, 'color', 'red', 'linewidth', 16);
end

