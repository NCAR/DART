% Data Assimilation Research Testbed -- DART
% Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
% Assumes two copies are ensemble mean followed by ensemble spread
% Should be automated and checked at some point

% Input the working directory
dir_name = input('Input directory; . for current directory')

% Load the ensemble file
tname = input('Input file name for True state');
fname = [dir_name '/' tname];
loc = getnc(fname, 'loc1d');
num_loc = size(loc, 1);
true_times = getnc(fname, 'time');
num_true_times = size(true_times, 1);

state_vec = getnc(fname, 'state');
state_size = size(state_vec, 2);

% Load the ensemble file
ename = input('Input file name for ensemble');
ens_fname = [dir_name '/' ename];
ens_full_vec = getnc(ens_fname, 'state');

% Ensemble size is
ens_size = size(ens_full_vec, 2);

% Num times in ensemble is
ens_times = getnc(ens_fname, 'time');
num_ens_times = size(ens_times, 1);
% Times is min of truth and ensemble
num_times = min(num_true_times, num_ens_times);

% Initialize storage for error averaging
rms(1:num_times) = 0.0;
sd_final(1:num_times) = 0.0;


% Loop through all the time levels
for time_ind = 1 : num_times

% Extract state and ensemble for just this time
single_state = state_vec(time_ind, :);

% Get ensemble mean
ens_mean = ens_full_vec(time_ind, 1, :);

% Get ensemble spread (S.D.) too
ens_sd = ens_full_vec(time_ind, 2, :);

% Compute and plot the difference field
ens_err = (ens_mean(:) - single_state(:));

% Compute statistics of the error field
max_err = max(max(ens_err));
min_err = min(min(ens_err));
rms_err = mean(mean(abs(ens_err)));
sd_mean = mean(mean(ens_sd));
rms(time_ind) = rms_err;
sd_final(time_ind) = sd_mean;

% End of time for loop
end


%Open a file for output of summary statistics from second half of run
fid = fopen([dir_name '/summary.ascii'], 'w');
fprintf(fid, ['Summary data for  ', dir_name, '\n']);

% Output some graphics of error as function of lead
figure(1);
subplot(2, 1, 1);
hold on;
plot(rms(:));
plot(sd_final(:), '--');
title ([dir_name ': Domain Average RMS Error']);
grid on;
subplot(2, 1, 2);
hold on;
plot(rms(num_times / 2 : num_times));
plot(sd_final(num_times / 2 : num_times), '--');
grid on;
% Compute the means over this time interval
mean_rms = mean(rms(num_times / 2 : num_times))
mean_sd = mean(sd_final(num_times / 2: num_times))
legend(strcat('rms mean sd = ', num2str(mean_rms), ' :: ', num2str(mean_sd)));
fprintf(fid, ['P ', num2str(mean_rms), '\n']);
