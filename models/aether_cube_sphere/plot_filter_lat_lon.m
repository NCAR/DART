% Plots the difference between aether_cube_sphere filter_input and
% filter_output files along with the input and output fields. 
% The user can select level, ensemble member,
% and one of four variables that are in the default test aether
% restart files. The variables are located at a set of vertical columns
% defined by the cube sphere grid. The plotting is quite primitive,
% just putting a symbol with a color scaled by the range of the 
% values on the plot at the given point. 
%
% This script is designed to be run in the aether_cube_sphere directory
% while the filter files are in the work directory below this.

% Level to be plotted
level = input('Input level to plot (integer)');

% Ensemble member to be plotted
ens = input('Input ensemble member to plot (integer)');
ens_prelim = int2str(10000 + ens);
ens_final = ens_prelim(2:5);

% Get a valid field index from user
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

% Create text strings for the input and output filter files
in_file_name = strcat('work/filter_input_', ens_final, '.nc');
out_file_name = strcat('work/filter_output_', ens_final, '.nc');

% Get the metadata latitude, longitude and altitude for each column
lat = ncread(in_file_name, 'lat');
lon = ncread(in_file_name, 'lon');
alt = ncread(in_file_name, 'alt');

% Read the field data from the input and output files
in_field = ncread(in_file_name, field_name);
out_field = ncread(out_file_name, field_name);

% Compute the increments
field = out_field - in_field;

% Get color range index so that scale spans the values in the field
range = max(field(:, level)) - min(field(:, level));
color_range = 256;
col = colormap;

% Loop to plot each point for the increments
for i = 1:size(lat, 1)
   frac(i) = (field(i, level) - min(field(:, level))) / range;
   color_index(i) = floor(frac(i) * color_range) + 1;
   color_index(i) = min(color_index(i), 256);
   plot(lon(i), lat(i), 'o', 'markersize', 10, 'color', col(color_index(i), :), 'linewidth', 28);
   set(gca, 'fontsize', 16)
   hold on
end

% Plot the input field in the same way
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

% Plot the output field in the same way
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

