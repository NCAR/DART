%% ens_scatter

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

fname = input('Input file name for truth');
%fname = 'True_State.nc';
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

% Load the ensemble file
ens_fname = input('Input file name for ensemble');
%ens_fname = 'Prior_Diag.nc'
ens_vec = getnc(ens_fname, 'state');

% Ensemble size is
ens_size = size(ens_vec, 2);

% Get ensemble mean
%ens_mean = mean(ens_vec(time_ind, :, :), 2);

% Select field to plot (ps, t, u, v)
field_num = input('Input field type, 1=ps, 2=t, 3=u, or 4=v')

% Get level for free atmosphere fields
if field_num > 1
   field_level = input('Input level');
else
   field_level = 1;
end

% Select x and y coordinates
x_coord = input('Select x coordinate');
y_coord = input('Select y coordinate');


figure(1);
close;
figure(1);
hold on;

% Extract ps or T fields
if field_num < 3
   offset = field_num + field_level - 1;

   field_vec = state_vec(:, offset : num_level + 1 : (num_level + 1) * (num_tlon * num_tlat));



% WARNING: MAKE SURE THAT DOING Y THEN X IN RESHAPE IS OKAY
   field = reshape(field_vec, [num_times, num_tlat num_tlon]);
   plot(field(:, y_coord, x_coord), 'r');
   
% Loop through the ensemble members
   for i = 1 : ens_size
      ens_member = ens_vec(:, i, offset : num_level + 1 : (num_level + 1) * (num_tlon * num_tlat));

      ens = reshape(ens_member, [num_times num_tlat num_tlon]);
      plot(ens(:, y_coord, x_coord));
   end

% Otherwise it's on v-grid
else

   base = (num_level + 1) * (num_tlon * num_tlat);
   offset = (field_level - 1) * 2 + (field_num - 2);
   field_vec = state_vec(:, base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);
   field = reshape(field_vec, [num_times num_vlat num_vlon]);
   plot(field(:, y_coord, x_coord), 'r');

% Loop through the ensemble members
   for i = 1 : ens_size
      ens_member = ens_vec(:, i, base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);
      ens = reshape(ens_member, [num_times num_vlat, num_vlon]);
      plot(ens(:, y_coord, x_coord));
   end

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
