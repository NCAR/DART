level = 2;
field_name = 'Temperature';
%field_name = 'O2+';
file_name = 'work/increment_0001.nc';

% Routine to plot DART filter file from Aether cube_sphere_grid
lat = ncread(file_name, 'lat');
lon = ncread(file_name, 'lon');
alt = ncread(file_name, 'alt');

field = ncread(file_name, field_name);

% Get color range index
range = max(field(:, level)) - min(field(:, level));
color_range = 256;
col = colormap;

% Loop to get colors
for i = 1:size(lat, 1)
   frac(i) = (field(i, level) - min(field(:, level))) / range;
   color_index(i) = floor(frac(i) * color_range) + 1;
   color_index(i) = min(color_index(i), 256);
   plot(lon(i), lat(i), 'o', 'markersize', 10, 'color', col(color_index(i), :), 'linewidth', 28);
   set(gca, 'fontsize', 16)
   hold on
end



