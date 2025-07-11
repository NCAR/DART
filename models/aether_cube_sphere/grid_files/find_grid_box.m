% Initial work on finding the actual 3 or 4 grid points that create triangle or quad that
% contain a given lat/lon on the Aether cubed sphere.

% Parameter that has the number of model grid points along each dimension 
% This does not include halos; the points are offset from the boundaries
np = 18;

% Some geometry about the grid point distribution
% Cube side is divided into 17 interior intervals of width 2sqrt(1/3) / np and
% two exterior intervals of half  width, sqrt(1/3) / np 
cube_side = 2 * sqrt(1 / 3); 
del = cube_side / np;
half_del = del / 2;

% Precompute the positions of the points along the cube
for i = 1:18
   x(i) = -sqrt(1/3) + half_del + del * (i - 1);
   box_lon(i) = atan2(sqrt(1 / 3), x(i));
end
box_lon = box_lon - pi / 4;

% Open grid files for each face for comparison
for i = 0:5
   % Generate the grid file name
   fname = strcat('grid_g000', num2str(i), '.nc');
   xt = ncread(fname, 'Longitude');
   glon(i + 1, :, :) = squeeze(xt(1, 3:end-2, 3:end-2));
   yt = ncread(fname, 'Latitude');
   glat(i + 1, :, :) = squeeze(yt(1, 3:end-2, 3:end-2));
end


% Enter lats and lons in degrees for starters
pt_lon_d = 50;
pt_lat_d = 89;
pt_lon = deg2rad(pt_lon_d);
pt_lat = deg2rad(pt_lat_d);
% Convert to x, y, z for plotting
px = cos(pt_lat) .* cos(pt_lon);
py = cos(pt_lat) .* sin(pt_lon);
pz = sin(pt_lat);

% Get the face, the longitudes on the two projections that don't have a pole in this face
% and the length along the two imbedded cube faces.

[face, lon_grid, len] = get_face(pt_lat, pt_lon); 

%[face lona lonb lena lenb]

% Figure out which interval this is in along each cube face; This gives 0 to np intervals
low_grid = floor((len + half_del) / del);
% This may not be right for the bottom and top faces; 
hi_grid = low_grid + 1;

% Get x, y, z for this face for plotting
x = cos(glat(face + 1, :, :)) .* cos(glon(face + 1, :, :));
y = cos(glat(face + 1, :, :)) .* sin(glon(face + 1, :, :));
z = sin(glat(face + 1, :, :));
plot3(squeeze(x), squeeze(y), squeeze(z), '.', 'markersize', 24, 'color', 'red', 'linewidth', 16);
hold on

% Plot the point in blue
plot3(px, py, pz, '.', 'markersize', 24, 'color', 'blue', 'linewidth', 16);

% Highlight the bounding points
bp_lon(1) = glon(face + 1, low_grid(2), low_grid(1)); bp_lat(1) = glat(face + 1, low_grid(2), low_grid(1));
bp_lon(2) = glon(face + 1, low_grid(2), hi_grid(1));  bp_lat(2) = glat(face + 1, low_grid(2), hi_grid(1));
bp_lon(3) = glon(face + 1, hi_grid(2),  low_grid(1)); bp_lat(3) = glat(face + 1, hi_grid(2),  low_grid(1));
bp_lon(4) = glon(face + 1, hi_grid(2),  hi_grid(1));  bp_lat(4) = glat(face + 1, hi_grid(2),  hi_grid(1));

% Convert each to x, y, z and plot
for i = 1:4
bpx = cos(bp_lat(i)) .* cos(bp_lon(i));
bpy = cos(bp_lat(i)) .* sin(bp_lon(i));
bpz = sin(bp_lat(i));
plot3(bpx, bpy, bpz, '.', 'markersize', 24, 'color', 'green', 'linewidth', 16);
   
end

% If it is 0 or np, we are on the boundary row
% Have to find neighbors in adjactent facer

% Otherwise, can compute bounding grid point indices





