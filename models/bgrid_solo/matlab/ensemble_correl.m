%% ens_error

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

fname = input('Input file name for True state');
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

% Select time
time = input('Input time level for correlations')

% Get the ensemble field for this time level only
ens_1t = ens_vec(time, :, :);

% Extract ps or T key point
if field_num < 3
   offset = field_num + field_level - 1;

   key_ens = ens_1t(1, :, offset : num_level + 1 : (num_level + 1) * (num_tlon * num_tlat));

   ens = reshape(key_ens, [ens_size num_tlat num_tlon]);
   key = ens(:, y_coord, x_coord);

% Otherwise it's on v-grid; extract key point
else

   base = (num_level + 1) * (num_tlon * num_tlat);
   offset = (field_level - 1) * 2 + (field_num - 2);
   key_ens = ens_1t(:, base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);
   ens = reshape(key_ens, [ens_size num_vlat, num_vlon]);
   key = ens(:, y_coord, x_coord)

end

% Compute correlations of key variable with each state var
% Set the size for the cor matrix
cor = ens_1t(1, 1, :);

for i = 1 : size(cor, 3)
   cor_mat = corrcoef(key, ens_1t(1, :, i));
   cor(i) = cor_mat(1, 2);
end





% Loop to plot a bunch of these
for i = 1 : 100

% Select field to plot (ps, t, u, v)
field_num = input('Input field type, 1=ps, 2=t, 3=u, or 4=v')

% Get level for free atmosphere fields
if field_num > 1
   field_level = input('Input level');
else
   field_level = 1;
end

% Extract ps or T fields
if field_num < 3
   offset = field_num + field_level - 1;

   cor_vec = cor(offset : num_level + 1 : (num_level + 1) * (num_tlon * num_tlat));
   fcor = reshape(cor_vec, [num_tlat num_tlon]);

% Otherwise it's on v-grid
else

   base = (num_level + 1) * (num_tlon * num_tlat);
   offset = (field_level - 1) * 2 + (field_num - 2);
   cor_vec = cor(base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);
   fcor = reshape(cor_vec, [num_vlat, num_vlon]);

end


% Plot the correlation
figure(1)
contour_vals = [-0.9 -0.8 -0.7 -0.6 -0.4  0.4 0.5 0.6 0.7 0.8 0.9];
[C, h] = contourf(fcor, contour_vals);
clabel(C, h);
colorbar;

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
