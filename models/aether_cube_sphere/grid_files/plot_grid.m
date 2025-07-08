
% Plot grid corners from an Aether cube sphere grid file section
x = ncread('grid_g0000.nc', 'Longitude');
lon = squeeze(x(1, 3:end-2, 3:end-2));
y = ncread('grid_g0000.nc', 'Latitude');
lat = squeeze(y(1, 3:end-2, 3:end-2));

plot(lon, lat, 'x', 'markersize', 10, 'color', 'blue');
daspect([1, 1, 1])
hold on

% Plot grid corners from an Aether cube sphere grid file section
x = ncread('grid_g0001.nc', 'Longitude');
lon = squeeze(x(1, 3:end-2, 3:end-2));
y = ncread('grid_g0001.nc', 'Latitude');
lat = squeeze(y(1, 3:end-2, 3:end-2));

plot(lon, lat, 'o', 'markersize', 10, 'color', 'r', 'color', 'r');



% Plot grid corners from an Aether cube sphere grid file section
x = ncread('grid_g0002.nc', 'Longitude');
lon = squeeze(x(1, 3:end-2, 3:end-2));
y = ncread('grid_g0002.nc', 'Latitude');
lat = squeeze(y(1, 3:end-2, 3:end-2));

plot(lon, lat, 'x', 'markersize', 10, 'color', 'r', 'color', 'b');




% Plot grid corners from an Aether cube sphere grid file section
x = ncread('grid_g0003.nc', 'Longitude');
lon = squeeze(x(1, 3:end-2, 3:end-2));
y = ncread('grid_g0003.nc', 'Latitude');
lat = squeeze(y(1, 3:end-2, 3:end-2));

plot(lon, lat, 'o', 'markersize', 10, 'color', 'r', 'color', 'r');


% Plot grid corners from an Aether cube sphere grid file section
x = ncread('grid_g0004.nc', 'Longitude');
lon = squeeze(x(1, 3:end-2, 3:end-2));
y = ncread('grid_g0004.nc', 'Latitude');
lat = squeeze(y(1, 3:end-2, 3:end-2));

plot(lon, lat, 'o', 'markersize', 10, 'color', 'r', 'color', 'g');


% Plot grid corners from an Aether cube sphere grid file section
x = ncread('grid_g0005.nc', 'Longitude');
lon = squeeze(x(1, 3:end-2, 3:end-2));
y = ncread('grid_g0005.nc', 'Latitude');
lat = squeeze(y(1, 3:end-2, 3:end-2));

plot(lon, lat, 'o', 'markersize', 10, 'color', 'r', 'color', 'g');


 
