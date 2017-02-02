%% cprior_post_err

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Assumes first two copies are ensemble mean and ensemble spread
% Should be automated

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

% Load the prior ensemble file
ens_fname = input('Input file name for prior ensemble');
%ens_fname = 'Prior_Diag.nc'
ens_full_vec = getnc(ens_fname, 'state');


% Load the posterior ensemble file
post_fname = input('Input file name for posterior ensemble');
post_full_vec = getnc(post_fname, 'state');

% Ensemble size is; assume same for prior posterior
ens_size = size(ens_full_vec, 2);

% Initialize storage for error averaging
rms(1:num_times, 1:4, 1:num_level) = 0.0;
post_rms(1:num_times, 1:4, 1:num_level) = 0.0;


% Loop through all the time levels
for time_ind = 1 : num_times

% Extract state and ensemble for just this time
single_state = state_vec(time_ind, :);

% Get ensemble mean
ens_mean = ens_full_vec(time_ind, 1, :);
post_mean = post_full_vec(time_ind, 1, :);


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
   post_vec = post_mean(offset : num_level + 1 : (num_level + 1) * (num_tlon * num_tlat));

   field = reshape(field_vec, [num_tlat num_tlon]);
   ens = reshape(ens_vec, [num_tlat num_tlon]);
   post = reshape(post_vec, [num_tlat num_tlon]);

% Otherwise it's on v-grid
else

   base = (num_level + 1) * (num_tlon * num_tlat);
   offset = (field_level - 1) * 2 + (field_num - 2);
   field_vec = single_state(base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);
   ens_vec = ens_mean(base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);
   post_vec = post_mean(base + offset : 2 * num_level : base + 2 * num_level * num_vlat * num_vlon);

   field = reshape(field_vec, [num_vlat, num_vlon]);
   ens = reshape(ens_vec, [num_vlat, num_vlon]);
   post = reshape(post_vec, [num_vlat, num_vlon]);

end


% Compute and plot the difference field
ens_err = (ens - field);
post_err = (post - field);

% Compute statistics of the error field
max_err = max(max(ens_err));
min_err = min(min(ens_err));
rms_err = mean(mean(abs(ens_err)));
rms(time_ind, field_num, field_level) = rms(time_ind, field_num, field_level) + rms_err;

% Compute statistics for the posterior error field, too
max_err = max(max(post_err));
min_err = min(min(post_err));
rms_err = mean(mean(abs(post_err)));
post_rms(time_ind, field_num, field_level) = post_rms(time_ind, field_num, field_level) + rms_err;


% End of level loop
end

% End of variable type loop
end

% End of time for loop
end


% Divide sum by number of time levels
rms = rms;
post_rms = post_rms;

% Output some graphics of error as function of lead
figure(1);
plot(rms(:, 1, 1));
title 'Error for ps'
hold on;
plot(post_rms(:, 1, 1), ':');

plot_temp(1:num_times, 1:num_level) = 0.0;
figure(2);
hold on;
plot_temp = reshape(rms(:, 2, :), [num_times num_level]);
plot(1:num_times, plot_temp);
title 'Errors for temperature'
plot_temp = reshape(post_rms(:, 2, :), [num_times num_level]);
plot(1:num_times, plot_temp, ':');

figure(3);
hold on;
plot_temp = reshape(rms(:, 3, :), [num_times num_level]);
plot(1:num_times, plot_temp);
title 'Errors for U'
plot_temp = reshape(post_rms(:, 3, :), [num_times num_level]);
plot(1:num_times, plot_temp, ':');

figure(4);
hold on;
plot_temp = reshape(rms(:, 4, :), [num_times num_level]);
plot(1:num_times, plot_temp);
title 'Errors for V'
plot_temp = reshape(post_rms(:, 4, :), [num_times num_level]);
plot(1:num_times, plot_temp, ':');


% Loop for another try
%ensemble;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
