% Plots aether grid corners on a cube plot
% Plot grid corners from an Aether cube sphere grid file section
x = ncread('grid_g0000.nc', 'Longitude');
lon = squeeze(x(1, 3:end-2, 3:end-2));
y = ncread('grid_g0000.nc', 'Latitude');
lat = squeeze(y(1, 3:end-2, 3:end-2));

% Convert points from lat/lon to x, y, z
x = cos(lat) .* cos(lon);
y = cos(lat) .* sin(lon);
z = sin(lat);

plot3(x, y, z, '.', 'markersize', 24, 'color', 'red', 'linewidth', 16);
hold on

% These have regular longitudes in base space
% Also have regular longitudes in the first rotated space
inv_sqrt_2 = 1 / sqrt(2);
for i = 1:size(x, 1);
   for j = 1:size(x, 2);
      vec = [x(i, j); y(i, j); z(i, j)];
      rot_vec = [0.5 0.5 -inv_sqrt_2; 0.5 0.5 inv_sqrt_2; inv_sqrt_2 -inv_sqrt_2 0] * vec;
      % Compute the longitude in the rotated space
      rlon(i, j) = atan2(rot_vec(2, :), rot_vec(1, :));
      if(rlon(i, j) < 0) rlon(i, j) = rlon(i, j) + 2*pi; end
   end
end

out0 = [lon(1, :)' rlon(:, 1)]

% Equidistant projection?
% Try 37 points total from -1 to 1 along a cube edge
% Corner of inscribed cube in radius 1 sphere is sqrt(1/3)
% Along an edge with constant z = sqrt(1/3), y = sqrt(1/3)
% Total of 37 x points equally spaced from - to + sqrt(1/3)
delx = 2 * sqrt(1/3) / 36;
for i = 1:37
   x = -sqrt(1/3) + delx * (i - 1);
   box_lon(i) = atan2(sqrt(1/3), x);
end

box_lon = box_lon - pi/4;


%------------------------------------------------------------------------

% Plot grid corners from an Aether cube sphere grid file section
x = ncread('grid_g0001.nc', 'Longitude');
lon = squeeze(x(1, 3:end-2, 3:end-2));
y = ncread('grid_g0001.nc', 'Latitude');
lat = squeeze(y(1, 3:end-2, 3:end-2));

% Convert points from lat/lon to x, y, z
x = cos(lat) .* cos(lon);
y = cos(lat) .* sin(lon);
z = sin(lat);

plot3(x, y, z, '.', 'markersize', 24, 'color', 'red', 'linewidth', 16);

% These have regular longitudes in the base space
% And in second rotated space
for i = 1:size(x, 1);
   for j = 1:size(x, 2);
      vec = [x(i, j); y(i, j); z(i, j)];
      rot_vec2 = [-0.5 0.5 -inv_sqrt_2; -0.5 0.5 inv_sqrt_2; inv_sqrt_2 inv_sqrt_2 0] * vec;
      % Compute the longitude in the rotated space
      rlon(i, j) = atan2(rot_vec2(2, :), rot_vec2(1, :));
      if(rlon(i, j) < 0) rlon(i, j) = rlon(i, j) + 2*pi; end
   end
end

out1 = [lon(1, :)' rlon(:, 1)]

%------------------------------------------------------------------------

% Plot grid corners from an Aether cube sphere grid file section
x = ncread('grid_g0002.nc', 'Longitude');
lon = squeeze(x(1, 3:end-2, 3:end-2));
y = ncread('grid_g0002.nc', 'Latitude');
lat = squeeze(y(1, 3:end-2, 3:end-2));

% Convert points from lat/lon to x, y, z
x = cos(lat) .* cos(lon);
y = cos(lat) .* sin(lon);
z = sin(lat);

plot3(x, y, z, '.', 'markersize', 24, 'color', 'red', 'linewidth', 16);

% These have regular longitudes in base space
% Also have regular longitudes in the first rotated space
inv_sqrt_2 = 1 / sqrt(2);
for i = 1:size(x, 1);
   for j = 1:size(x, 2);
      vec = [x(i, j); y(i, j); z(i, j)];
      rot_vec = [0.5 0.5 -inv_sqrt_2; 0.5 0.5 inv_sqrt_2; inv_sqrt_2 -inv_sqrt_2 0] * vec;
      % Compute the longitude in the rotated space
      rlon(i, j) = atan2(rot_vec(2, :), rot_vec(1, :));
      if(rlon(i, j) < 0) rlon(i, j) = rlon(i, j) + 2*pi; end
   end
end

out2 = [lon(1, :)' rlon(:, 1)]

%------------------------------------------------------------------------

% Plot grid corners from an Aether cube sphere grid file section
x = ncread('grid_g0003.nc', 'Longitude');
lon = squeeze(x(1, 3:end-2, 3:end-2));
y = ncread('grid_g0003.nc', 'Latitude');
lat = squeeze(y(1, 3:end-2, 3:end-2));

% Convert points from lat/lon to x, y, z
x = cos(lat) .* cos(lon);
y = cos(lat) .* sin(lon);
z = sin(lat);

plot3(x, y, z, '.', 'markersize', 24, 'color', 'red', 'linewidth', 16);

% These have regular longitudes in the base space
% And in second rotated space
for i = 1:size(x, 1);
   for j = 1:size(x, 2);
      vec = [x(i, j); y(i, j); z(i, j)];
      rot_vec2 = [-0.5 0.5 -inv_sqrt_2; -0.5 0.5 inv_sqrt_2; inv_sqrt_2 inv_sqrt_2 0] * vec;
      % Compute the longitude in the rotated space
      rlon(i, j) = atan2(rot_vec2(2, :), rot_vec2(1, :));
      if(rlon(i, j) < 0) rlon(i, j) = rlon(i, j) + 2*pi; end
   end
end

out3 = [lon(1, :)' rlon(:, 1)]

%------------------------------------------------------------------------

% Plot grid corners from an Aether cube sphere grid file section
x = ncread('grid_g0004.nc', 'Longitude');
lon = squeeze(x(1, 3:end-2, 3:end-2));
y = ncread('grid_g0004.nc', 'Latitude');
lat = squeeze(y(1, 3:end-2, 3:end-2));

% Convert points from lat/lon to x, y, z
x = cos(lat) .* cos(lon);
y = cos(lat) .* sin(lon);
z = sin(lat);

plot3(x, y, z, '.', 'markersize', 24, 'color', 'red', 'linewidth', 16);

% This should have regular longitudes in both of the rotated spaces

% Also have regular longitudes in the first rotated space
inv_sqrt_2 = 1 / sqrt(2);
for i = 1:size(x, 1);
   for j = 1:size(x, 2);
      vec = [x(i, j); y(i, j); z(i, j)];
      rot_vec = [0.5 0.5 -inv_sqrt_2; 0.5 0.5 inv_sqrt_2; inv_sqrt_2 -inv_sqrt_2 0] * vec;
      % Compute the longitude in the rotated space
      rlon(i, j) = atan2(rot_vec(2, :), rot_vec(1, :));
      if(rlon(i, j) < 0) rlon(i, j) = rlon(i, j) + 2*pi; end
   end
end

% And in second rotated space
for i = 1:size(x, 1);
   for j = 1:size(x, 2);
      vec = [x(i, j); y(i, j); z(i, j)];
      rot_vec2 = [-0.5 0.5 -inv_sqrt_2; -0.5 0.5 inv_sqrt_2; inv_sqrt_2 inv_sqrt_2 0] * vec;
      % Compute the longitude in the rotated space
      rlon2(i, j) = atan2(rot_vec2(2, :), rot_vec2(1, :));
      if(rlon2(i, j) < 0) rlon2(i, j) = rlon2(i, j) + 2*pi; end
   end
end

out4 = [rlon(1, :)' rlon2(:, 1)]

%------------------------------------------------------------------------

% Plot grid corners from an Aether cube sphere grid file section
x = ncread('grid_g0005.nc', 'Longitude');
lon = squeeze(x(1, 3:end-2, 3:end-2));
y = ncread('grid_g0005.nc', 'Latitude');
lat = squeeze(y(1, 3:end-2, 3:end-2));

% Convert points from lat/lon to x, y, z
x = cos(lat) .* cos(lon);
y = cos(lat) .* sin(lon);
z = sin(lat);

plot3(x, y, z, '.', 'markersize', 24, 'color', 'red', 'linewidth', 16);

% This should have regular longitudes in both of the rotated spaces

% Also have regular longitudes in the first rotated space
inv_sqrt_2 = 1 / sqrt(2);
for i = 1:size(x, 1);
   for j = 1:size(x, 2);
      vec = [x(i, j); y(i, j); z(i, j)];
      rot_vec = [0.5 0.5 -inv_sqrt_2; 0.5 0.5 inv_sqrt_2; inv_sqrt_2 -inv_sqrt_2 0] * vec;
      % Compute the longitude in the rotated space
      rlon(i, j) = atan2(rot_vec(2, :), rot_vec(1, :));
      if(rlon(i, j) < 0) rlon(i, j) = rlon(i, j) + 2*pi; end
   end
end

% And in second rotated space
for i = 1:size(x, 1);
   for j = 1:size(x, 2);
      vec = [x(i, j); y(i, j); z(i, j)];
      rot_vec2 = [-0.5 0.5 -inv_sqrt_2; -0.5 0.5 inv_sqrt_2; inv_sqrt_2 inv_sqrt_2 0] * vec;
      % Compute the longitude in the rotated space
      rlon2(i, j) = atan2(rot_vec2(2, :), rot_vec2(1, :));
      if(rlon2(i, j) < 0) rlon2(i, j) = rlon2(i, j) + 2*pi; end
   end
end

out5 = [rlon(1, :)' rlon2(:, 1)]

%------------------------------------------------------------------------
 
