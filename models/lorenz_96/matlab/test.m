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

end
