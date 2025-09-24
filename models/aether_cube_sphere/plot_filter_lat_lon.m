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
      field_name = 'O2+'; break
   elseif(field_num == 2) 
      field_name = 'O2'; break
   elseif(field_num == 3)
      field_name = 'N2+'; break
   elseif(field_num == 4)
      field_name = 'Temperature'; break
   end
end

in_file_name = strcat('work/filter_input_', ens_final, '.nc');
out_file_name = strcat('work/filter_output_', ens_final, '.nc');

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



