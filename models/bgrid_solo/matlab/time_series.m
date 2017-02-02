%% time_series

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

fname = input('Input file name');
tlon = getnc(fname, 'TmpI');
num_tlon = size(tlon, 1);
tlat = getnc(fname, 'TmpJ');
num_tlat = size(tlat, 1);
vlon = getnc(fname, 'VelI');
num_vlon = size(vlon, 1);
vlat = getnc(fname, 'VelJ');
num_vlat = size(vlat, 1);
level = getnc(fname, 'level');
num_level = size(level, 1);
times = getnc(fname, 'time');
num_times = size(times, 1);

state_vec = getnc(fname, 'state');

% Select field to plot (ps, t, u, v)
field_num = input('Input field type, 1=ps, 2=t, 3=u, or 4=v')

% Select x and y coordinates
x_coord = 15;


close all;

if field_num < 2 
   plot_num_level = 1;
else
   plot_num_level = num_level;
end

% Loop to do this for all levels 
for field_level = 1:plot_num_level
figure(field_level);
hold on;

% Extract ps or T fields
if field_num < 3
   offset = field_num + field_level - 1;

   field_vec = state_vec(:, offset : num_level + 1 : (num_level + 1) * (num_tlon * num_tlat));

% WARNING: MAKE SURE THAT DOING Y THEN X IN RESHAPE IS OKAY
   field = reshape(field_vec, [num_times, num_tlat num_tlon]);
   plot(field(:, 12:3:18, x_coord));
   
% Otherwise it's on v-grid
else

   base = (num_level + 1) * (num_tlon * num_tlat);
   offset = (field_level - 1) * 2 + (field_num - 2);
   field_vec = state_vec(:, base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);
   field = reshape(field_vec, [num_times num_vlat num_vlon]);
   plot(field(:, 12:3:18, x_coord));

end
end   %level loop

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
