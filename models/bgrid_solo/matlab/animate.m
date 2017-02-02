%% animate ... makes a matlab movie object of a field ... interactively gathers input

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Assumes 2 copies of data are ensemble mean and spread
% Should be checked and automated


%fname = 'True_State.nc';
fname = input('Input true state name');
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
%time_ind = input('Input time level');
for time_ind = 1:40

% Extract state and ensemble for just this time
single_state = state_vec(time_ind, :);

% Get ensemble mean and spread
ens_mean = ens_vec(time_ind, 1, :);
ens_spread = ens_vec(time_ind, 2, :);

% Select field to plot (ps, t, u, v)
%field_num = input('Input field type, 1=ps, 2=t, 3=u, or 4=v')
field_num = 2;

% Get level for free atmosphere fields
if field_num > 1
%   field_level = input('Input level');
   field_level = 3;
else
   field_level = 1;
end

% Extract ps or T fields
if field_num < 3
   offset = field_num + field_level - 1;

   field_vec = single_state(offset : num_level + 1 : (num_level + 1) * (num_tlon * num_tlat));
   ens_vec2 = ens_mean(offset : num_level + 1 : (num_level + 1) * (num_tlon * num_tlat));
   spread_vec = ens_spread(offset : num_level + 1 : (num_level + 1) * (num_tlon * num_tlat));

   field = reshape(field_vec, [num_tlat num_tlon]);
   ens = reshape(ens_vec2, [num_tlat num_tlon]);
   spread = reshape(spread_vec, [num_tlat num_tlon]);

% Otherwise it's on v-grid
else

   base = (num_level + 1) * (num_tlon * num_tlat);
   offset = (field_level - 1) * 2 + (field_num - 2);
   field_vec = single_state(base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);
   ens_vec2 = ens_mean(base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);
   spread_vec = ens_spread(base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);

   field = reshape(field_vec, [num_vlat, num_vlon]);
   ens = reshape(ens_vec2, [num_vlat, num_vlon]);
   spread = reshape(spread_vec, [num_vlat, num_vlon]);

end


%figure(time_ind);
figure(1);
[C, h] = contourf(field, 10);
%clabel(C, h);
colorbar;
M(time_ind) = getframe;

%figure(2);
%[C, h] = contourf(ens, 10);
%clabel(C, h);

% Compute and plot the difference field
%ens_err = (ens - field);
%figure(3);
%[C, h] = contourf(ens_err, 10);
%clabel(C, h);

% Compute statistics of the error field
%max_err = max(max(ens_err));
%min_err = min(min(ens_err));
%rms_err = mean(mean(abs(ens_err)));

% Label figure 3 with these statistics
%title_string = ['Min = ', num2str(min_err), ' Max =  ', num2str(max_err), '   RMS ERROR = ', num2str(rms_err)];
%title (title_string)

% Output the spread plot, too
%figure(4);
%[C, h] = contourf(spread, 10);
%clabel(C, h);

% Compute statistics of the spread field
%max_spread = max(max(ens_spread));
%min_spread = min(min(ens_spread));
%rms_spread = mean(mean(ens_spread));

% Label figure 4 with these statistics
%title_string = ['Min = ', num2str(min_spread), ' Max =  ', num2str(max_spread), '   RMS ERROR = ', num2str(rms_spread)];
%title (title_string)


end

Movie(M);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
