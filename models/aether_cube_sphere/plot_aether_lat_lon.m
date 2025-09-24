% Plots aether file increments

% Plot values from aether neutrals block files on sphere
nblocks = 6;
nhalos = 2;

% Level to be plotted
level = input('Input level to plot (integer)');

% Ensemble member to be plotted
ens = input('Input ensemble member to plot (integer)');
ens_prelim = int2str(10000 + ens);
ens_final = ens_prelim(2:5);

% Get a valid field_num from user
for i = 1:1000
   field_num = 0;
   % Input field name
   field_num = input('Input 1 for O2+; 2 for O2; 3 for N2+; 4 for Temp.');

   if(field_num == 1)
      file_type = 1;
      field_name = 'O2+'; break
   elseif(field_num == 2)
      file_type = 2;
      field_name = 'O2'; break
   elseif(field_num == 3)
      file_type = 1;
      field_name = 'N2+'; break
   elseif(field_num == 4)
      file_type = 2;
      field_name = 'Temperature'; break
   end
end

% Input file directory (grid files must be here, too
input_base_name = 'TEST_INPUT/';
% Output file directory 
output_base_name = 'TEST_OUTPUT/';

g_base_name = strcat(input_base_name, 'grid_g');
if(file_type == 1)
   input_base_name  = strcat(input_base_name,  'ions_m', ens_final, '_g');
   output_base_name = strcat(output_base_name, 'ions_m', ens_final, '_g');
else
   input_base_name  = strcat(input_base_name,  'neutrals_m', ens_final, '_g');
   output_base_name = strcat(output_base_name, 'neutrals_m', ens_final, '_g');
end

% Loop through the fields from each block first to get min and max for plot colors
for i = 0:nblocks - 1
   num_prelim = int2str(10000 + i);
   num_end = num_prelim(2:5);

   input_file_name = strcat(input_base_name, num_end, '.nc');
   output_file_name = strcat(output_base_name, num_end, '.nc');

   % Read in the fields from input and output files
   in_field = ncread(input_file_name, field_name);
   out_field = ncread(output_file_name, field_name);
   field = out_field - in_field;

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
   g_file_name = strcat(g_base_name, num_end, '.nc');

   input_file_name = strcat(input_base_name, num_end, '.nc');
   output_file_name = strcat(output_base_name, num_end, '.nc');

   % Read the grid file lat and lons
   lat = ncread(g_file_name, 'Latitude');
   lon = ncread(g_file_name, 'Longitude');
   
   % Read in the fields from input and output files
   in_field = ncread(input_file_name, field_name);
   out_field = ncread(output_file_name, field_name);
   field = out_field - in_field;

   % Get color range index
   range = n_max - n_min;
   color_range = 256;
   col = colormap;

   % Loop to get colors
   for i = nhalos+1:size(field, 2) - nhalos
      for j = nhalos+1:size(field, 3) - nhalos
         frac = (field(level, i, j) - n_min) / range;
         color_index = floor(frac * color_range) + 1;
         color_index = min(color_index, 256);
         %plot(x(i, j), y(i, j), z(i, j), 'o', 'markersize', 24, 'color', col(color_index, :), 'linewidth', 16);
         plot(lon(level, i, j), lat(level, i, j), 'o', 'markersize', 10, 'color', col(color_index, :), 'linewidth', 20);
         set(gca, 'fontsize', 16)
         hold on
      end
   end
end
