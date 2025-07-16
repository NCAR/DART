% Initial work on finding the actual 3 or 4 grid points that create triangle or quad that
% contain a given lat/lon on the Aether cubed sphere.

% Parameter that has the number of model grid points along each dimension 
% This does not include halos; the points are offset from the boundaries
np = 18;

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

% Compare computed lat lon to grid files for a given face
face = 4;
for i = 1:18
   for j = 1:18
      [dlat(j, i), dlon(j, i)] = grid_to_lat_lon(face, j, i, np);
   end
end

% Plot the d grids
   % Get x, y, z for this face for plotting
   x = cos(dlat) .* cos(dlon);
   y = cos(dlat) .* sin(dlon);
   z = sin(dlat);
   % Plot points on each face with distinct colors 
   %%%h = plot3(x, y, z, '.', 'markersize', 24, 'color', 'red', 'linewidth', 16);

% Compute the max of the differences
max(max(abs(squeeze(glon(face + 1, :, :)) - dlon)))
max(max(abs(squeeze(glat(face + 1, :, :)) - dlat)))

%-------------------------------------------

% Enter lats and lons in degrees for starters
pt_lon_d = 269.7;
pt_lat_d = -29;
% Convert to radians since this is generally used in algorithm
pt_lon = deg2rad(pt_lon_d);
pt_lat = deg2rad(pt_lat_d);
% Convert to x, y, z for plotting
px = cos(pt_lat) .* cos(pt_lon);
py = cos(pt_lat) .* sin(pt_lon);
pz = sin(pt_lat);
% Plot the point in green
plot3(px, py, pz, '.', 'markersize', 24, 'color', 'green', 'linewidth', 16);


[grid_face, grid_lat_ind, grid_lon_ind, grid_pt_lat, grid_pt_lon, corner_detected] = ...
   get_bounding_box(pt_lat, pt_lon, np);

if(corner_detected)
   num_bound_points = 3;
else
   num_bound_points = 4;
end

% Highlight the bounding points;
for i = 1:num_bound_points
   bp_lon(i) = glon(grid_face(i) + 1, grid_lat_ind(i), grid_lon_ind(i));
   bp_lat(i) = glat(grid_face(i) + 1, grid_lat_ind(i), grid_lon_ind(i));

   % Convert each to x, y, z and plot
   bpx = cos(bp_lat(i)) .* cos(bp_lon(i));
   bpy = cos(bp_lat(i)) .* sin(bp_lon(i));
   bpz = sin(bp_lat(i));
   % Plot the bounding points in red
   plot3(bpx, bpy, bpz, '.', 'markersize', 24, 'color', 'red', 'linewidth', 16);

   % Also plot the lat lon pairs determined automatically
   % Convert each to x, y, z and plot
   cpx = cos(grid_pt_lat(i)) .* cos(grid_pt_lon(i));
   cpy = cos(grid_pt_lat(i)) .* sin(grid_pt_lon(i));
   cpz = sin(grid_pt_lat(i));
   % Plot the bounding points in red
   plot3(cpx, cpy, cpz, 'o', 'markersize', 16, 'color', 'red', 'linewidth', 3);

end











