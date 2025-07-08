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

%------------------------------------------------------------------------
 
