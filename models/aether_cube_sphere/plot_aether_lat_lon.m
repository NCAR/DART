% Plot values from aether neutrals block files on sphere
%nblocks = 24;
%field_name = 'O';
%level = 11;
nblocks = 6;
nhalos = 2;
%field_name = 'O2+';
field_name = 'Temperature';
level = 31;

g_base_name = 'B6_AETHER_INPUT_FILES/restartOut/grid_g';
n_base_name = 'increments/neutrals_inc_m0009_g';
%n_base_name = 'increments/ions_inc_m0009_g';

% Loop through the fields from each block first to get min and max for plot colors
for i = 0:nblocks - 1
   num_prelim = int2str(10000 + i);
   num_end = num_prelim(2:5);
   n_file_name = strcat(n_base_name, num_end, '.nc')

   % Read in the field from the neutrals file
   field = ncread(n_file_name, field_name);
   nx = size(field, 3);
   xs = nhalos + 1;
   xe = nx - nhalos;
   
   if(i == 0) 
      n_min = min(min(field(level, xs:xe, xs:xe)));
      n_max = max(max(field(level, xs:xe, xs:xe)));
   else
      t_min = min(min(field(level, xs:xe, xs:xe)));
      t_max = max(max(field(level, xs:xe, xs:xe)));
      n_min = min(t_min, n_min);
      n_max = max(t_max, n_max);
   end
end

for i = 0:nblocks -1 
   num_prelim = int2str(10000 + i);
   num_end = num_prelim(2:5);
   g_file_name = strcat(g_base_name, num_end, '.nc')
   n_file_name = strcat(n_base_name, num_end, '.nc')

   % Read the grid file lat and lons
   lat = ncread(g_file_name, 'Latitude');
   lon = ncread(g_file_name, 'Longitude');
   
   % Read in the field from the neutrals file
   field = ncread(n_file_name, field_name);

   % Begin by confirming that the grids look like cube sphere
   % Convert points from lat/lon to x, y, z
   x = squeeze(cos(lat(level, :, :)) .* cos(lon(level, :, :)));
   y = squeeze(cos(lat(level, :, :)) .* sin(lon(level, :, :)));
   z = squeeze(sin(lat(level, :, :)));

   %WARNING: NOT MAX OF ALL BLOCKS
   % Get color range index
   range = n_max - n_min;
   color_range = 256;
   col = colormap;

   % Loop to get colors
   for i = nhalos+1:size(x, 1) - nhalos
      for j = nhalos+1:size(x, 2) - nhalos
         frac = (field(level, i, j) - n_min) / range
         color_index = floor(frac * color_range) + 1;
         color_index = min(color_index, 256);
         %plot(x(i, j), y(i, j), z(i, j), 'o', 'markersize', 24, 'color', col(color_index, :), 'linewidth', 16);
         plot(lon(level, i, j), lat(level, i, j), 'o', 'markersize', 10, 'color', col(color_index, :), 'linewidth', 20);
         set(gca, 'fontsize', 16)
         hold on
      end
   end
end
