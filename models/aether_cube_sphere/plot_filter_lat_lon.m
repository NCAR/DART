level = 9;
%field_name = 'O+';
%field_name = 'O2+';
%field_name = 'Temperature';
field_name = 'O';
in_file_name = 'work/filter_input_0001.nc';
out_file_name = 'work/filter_output_0001.nc';

% Routine to plot DART filter file from Aether cube_sphere_grid
lat = ncread(in_file_name, 'lat');
lon = ncread(in_file_name, 'lon');
alt = ncread(in_file_name, 'alt');

in_field = ncread(in_file_name, field_name);
out_field = ncread(out_file_name, field_name);

field = out_field - in_field;

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

figure(2)
field = in_field;
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

figure(3)
field = out_field;
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



