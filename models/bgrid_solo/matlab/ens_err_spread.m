%% ens_err_spread

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

fname = 'True_State.nc';
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
state_size = size(state_vec, 2);

% Load the ensemble file
ens_fname = input('Input file name for ensemble');
%ens_fname = 'Prior_Diag.nc'
ens_full_vec = getnc(ens_fname, 'state');

% Ensemble size is
ens_size = size(ens_full_vec, 2);

% Initialize storage for error averaging
rms(1:num_times, 1:4, 1:num_level) = 0.0;
sd_final(1:num_times, 1:4, 1:num_level) = 0.0;


% Loop through all the time levels
for time_ind = 1 : num_times

% Extract state and ensemble for just this time
single_state = state_vec(time_ind, :);

% Get ensemble mean
ens_mean = mean(ens_full_vec(time_ind, :, :), 2);

% Get ensemble spread (S.D.) too
ens_sd(1:state_size) = 0.0;
for i = 1 : state_size
   ens_sd(i) = sqrt(var(ens_full_vec(time_ind, :, i)));
end

% Loop through each of the variable types 
for field_num = 1:4

% Get level for free atmosphere fields
if field_num > 1
   max_level = num_level;
else
   max_level = 1;
end

% Loop through all levels
for field_level = 1 : max_level

% Extract ps or T fields
if field_num < 3
   offset = field_num + field_level - 1;

   field_vec = single_state(offset : num_level + 1 : (num_level + 1) * (num_tlon * num_tlat));
   ens_vec = ens_mean(offset : num_level + 1 : (num_level + 1) * (num_tlon * num_tlat));
   sd_vec = ens_sd(offset : num_level + 1 : (num_level + 1) * (num_tlon * num_tlat));

   field = reshape(field_vec, [num_tlat num_tlon]);
   ens = reshape(ens_vec, [num_tlat num_tlon]);
   sd = reshape(sd_vec, [num_tlat num_tlon]);

% Otherwise it's on v-grid
else

   base = (num_level + 1) * (num_tlon * num_tlat);
   offset = (field_level - 1) * 2 + (field_num - 2);
   field_vec = single_state(base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);
   ens_vec = ens_mean(base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);
   sd_vec = ens_sd(base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);

   field = reshape(field_vec, [num_vlat, num_vlon]);
   ens = reshape(ens_vec, [num_vlat, num_vlon]);
   sd = reshape(sd_vec, [num_vlat, num_vlon]);

end


% Compute and plot the difference field
ens_err = (ens - field);

% Compute statistics of the error field
max_err = max(max(ens_err));
min_err = min(min(ens_err));
rms_err = mean(mean(abs(ens_err)));
sd_mean = mean(mean(sd));
rms(time_ind, field_num, field_level) = rms(time_ind, field_num, field_level) + rms_err;
sd_final(time_ind, field_num, field_level) = sd_final(time_ind, field_num, field_level) + sd_mean;


% End of level loop
end

% End of variable type loop
end

% End of time for loop
end

% Output some graphics of error as function of lead
figure(1);
hold on;
plot(rms(:, 1, 1));
plot(sd_final(:, 1, 1), ':');
title 'Error for ps'

plot_temp(1:num_times, 1:num_level) = 0.0;
figure(2);
hold on;
plot_temp = reshape(rms(:, 2, :), [num_times num_level]);
plot(1:num_times, plot_temp);
plot_temp = reshape(sd_final(:, 2, :), [num_times num_level]);
plot(1:num_times, plot_temp, ':');
title 'Errors for temperature'

figure(3);
hold on;
plot_temp = reshape(rms(:, 3, :), [num_times num_level]);
plot(1:num_times, plot_temp);
plot_temp = reshape(sd_final(:, 3, :), [num_times num_level]);
plot(1:num_times, plot_temp, ':');
title 'Errors for U'

figure(4);
hold on;
plot_temp = reshape(rms(:, 4, :), [num_times num_level]);
plot(1:num_times, plot_temp);
plot_temp = reshape(sd_final(:, 4, :), [num_times num_level]);
plot(1:num_times, plot_temp, ':');
title 'Errors for V'


% Loop for another try
%ensemble;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
