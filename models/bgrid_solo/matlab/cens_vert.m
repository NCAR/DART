%% cens_vert

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Assumes 2 copies of data are ensemble mean and spread
% Should be checked and automated

fname = input('Input file name for true state')
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

state_vec = getnc(fname, 'state');

% Load the ensemble file
ens_fname = input('Input file name for ensemble');
%ens_fname = 'Prior_Diag.nc'
ens_vec = getnc(ens_fname, 'state');

% Ensemble size is
ens_size = size(ens_vec, 2);

% Get a time level from the user
time_ind = input('Input time level');

% Extract state and ensemble for just this time
single_state = state_vec(time_ind, :);

% Get ensemble mean and spread
ens_mean = ens_vec(time_ind, 1, :);
ens_spread = ens_vec(time_ind, 2, :);

% Select field to plot (ps, t, u, v)
field_num = input('Input field type, 2=t, 3=u, or 4=v')

% Extract ps or T fields
if field_num == 2

% Get an array for the vertical cross section
   field_vert(1:num_level, 1:num_tlat) = 0.0;
   ens_vert(1:num_level, 1:num_tlat) = 0.0;
   spread_vert(1:num_level, 1:num_tlat) = 0.0;

   for field_level = 1 : num_level
      offset = field_num + field_level - 1;
      field_vec = single_state(offset : num_level + 1 : (num_level + 1) * (num_tlon * num_tlat));
      field = reshape(field_vec, [num_tlat, num_tlon]);
      field_vert(num_level - field_level + 1, :) = transpose(mean(field, 2));

      ens_vec = ens_mean(offset : num_level + 1 : (num_level + 1) * (num_tlon * num_tlat));
      ens = reshape(ens_vec, [num_tlat num_tlon]);
      ens_vert(num_level - field_level + 1, :) = transpose(mean(ens, 2));

      spread_vec = ens_spread(offset : num_level + 1 : (num_level + 1) * (num_tlon * num_tlat));
      spread = reshape(spread_vec, [num_tlat num_tlon]);
      spread_vert(num_level - field_level + 1, :) = transpose(mean(spread, 2));
   end

% Otherwise it's on v-grid
else
   field_vert(1:num_level, 1:num_vlat) = 0.0;
   ens_vert(1:num_level, 1:num_vlat) = 0.0;
   spread_vert(1:num_level, 1:num_vlat) = 0.0;

   for field_level = 1 : num_level
      base = (num_level + 1) * (num_tlon * num_tlat);
      offset = (field_level - 1) * 2 + (field_num - 2);
%      offset = (field_num - 3) * num_level + field_level;
      field_vec = single_state(base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);
      field = reshape(field_vec, [num_vlat, num_vlon]);
      field_vert(num_level - field_level + 1, :) = transpose(mean(field, 2));

      ens_vec = ens_mean(base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);
      ens = reshape(ens_vec, [num_vlat, num_vlon]);
      ens_vert(num_level - field_level + 1, :) = transpose(mean(ens, 2));

      spread_vec = ens_spread(base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);
      spread = reshape(spread_vec, [num_vlat, num_vlon]);
      spread_vert (num_level - field_level + 1, :) = transpose(mean(spread, 2));
   end

end


figure(1);
[C, h] = contourf(field_vert);
clabel(C, h);

figure(2);
[C, h] = contourf(ens_vert);
clabel(C, h);

% Compute and plot the difference field
ens_err = (ens_vert - field_vert);
figure(3);
[C, h] = contourf(ens_err);
clabel(C, h);

% Compute statistics of the error field
max_err = max(max(ens_err));
min_err = min(min(ens_err));
rms_err = mean(mean(abs(ens_err)));

% Label figure 3 with these statistics
title_string = ['Min = ', num2str(min_err), ' Max =  ', num2str(max_err), '   RMS ERROR = ', num2str(rms_err)];
title (title_string)

% Output the spread plot, too
figure(4);
[C, h] = contourf(spread_vert);
clabel(C, h);

% Compute statistics of the spread field
max_spread = max(max(spread_vert));
min_spread = min(min(spread_vert));
rms_spread = mean(mean(spread_vert));

% Label figure 4 with these statistics
title_string = ['Min = ', num2str(min_spread), ' Max =  ', num2str(max_spread), '   RMS ERROR = ', num2str(rms_spread)];
title (title_string)

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
