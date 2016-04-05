% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$
 
% Assumes two copies are ensemble mean followed by ensemble spread
% Should be automated and checked at some point

% Input the working directory
dir_name = input('Input directory; . for current directory')

% Load the ensemble file
tname = input('Input file name for True state');
fname = [dir_name '/' tname];
lon = getnc(fname, 'lon');
num_lon = size(lon, 1);
lat = getnc(fname, 'lat');
num_lat = size(lat, 1);
level = getnc(fname, 'lev');
num_level = size(level, 1);
true_times = getnc(fname, 'time');
num_true_times = size(true_times, 1);

%state_vec = getnc(fname, 'state');
%state_size = size(state_vec, 2);

% Load the ensemble file
ename = input('Input file name for ensemble');
ens_fname = [dir_name '/' ename];
%ens_full_vec = getnc(ens_fname, 'state');

% Num times in ensemble is
ens_times = getnc(ens_fname, 'time');
num_ens_times = size(ens_times, 1);
% Times is min of truth and ensemble
num_times = min(num_true_times, num_ens_times);

% Initialize storage for error averaging
red_num_level = 38;
rms(1:num_times, 1:6, 1:red_num_level) = 0.0;
sd_final(1:num_times, 1:9, 1:red_num_level) = 0.0;

% Can skip the first three fields (u1, v1, t1) for now 
%field_name = ['u', 'v', 't', 'qnH', 'qnOH', 'qnO']

% Loop through each of the variable types 
for field_num = 1:6

if field_num == 1 field_name = 'u'    ;end
if field_num == 2 field_name = 'v'    ;end
if field_num == 3 field_name = 't'    ;end
if field_num == 4 field_name = 'qnH'  ;end
if field_num == 5 field_name = 'qnOH' ;end
if field_num == 6 field_name = 'qnO'  ;end

ens_vec   = getnc(ens_fname, field_name);
state_vec = getnc(fname, field_name);

% Loop through all the time levels
for time_ind = 1 : num_times
% Loop through all levels
for field_level = 1 : red_num_level

   field =  reshape(state_vec(time_ind,  :, :, field_level),[num_lat,num_lon]); %truth[lat, lon]
   ens =    reshape(ens_vec(time_ind, 11, :, :, field_level),[num_lat,num_lon]);
   sd =     reshape(ens_vec(time_ind, 12, :, :, field_level),[num_lat,num_lon]); 

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

%Open a file for output of summary statistics from second half of run
fid = fopen([dir_name '/summary.ascii'], 'w');
fprintf(fid, ['Summary data for  ', dir_name, '\n']);

% Output some graphics of error as function of lead
figure(1);
subplot(2, 1, 1);
hold on;
plot_temp = reshape(rms(:, 1, :), [num_times red_num_level]);
plot(1:num_times, plot_temp);
plot_temp = reshape(sd_final(:, 1, :), [num_times red_num_level]);
plot(1:num_times, plot_temp, '--');
title ([dir_name ': Errors for U']);
grid on;
subplot(2, 1, 2);
hold on;
plot_temp = reshape(rms(:, 1, :), [num_times red_num_level]);
mean_rms = mean(plot_temp(num_times / 2:num_times, :), 1)
plot(num_times / 2:num_times, plot_temp(num_times/2: num_times, :));
plot_temp = reshape(sd_final(:, 1, :), [num_times red_num_level]);
mean_sd = mean(plot_temp(num_times / 2:num_times, :), 1)
plot(num_times / 2:num_times, plot_temp(num_times/2:num_times, :), '--');
grid on;
legend(num2str(mean_rms), num2str(mean_sd));
fprintf(fid, ['U ', num2str(mean_rms), '\n']);
print_file = strcat(dir_name, '/u_ts.eps')
%print(gcf, '-depsc', print_file);

figure(2);
subplot(2, 1, 1);
hold on;
plot_temp = reshape(rms(:, 2, :), [num_times red_num_level]);
plot(1:num_times, plot_temp);
plot_temp = reshape(sd_final(:, 2, :), [num_times red_num_level]);
plot(1:num_times, plot_temp, '--');
title ([dir_name ': Errors for V']);
grid on;
subplot(2, 1, 2);
hold on;
plot_temp = reshape(rms(:, 2, :), [num_times red_num_level]);
mean_rms = mean(plot_temp(num_times / 2:num_times, :), 1)
plot(num_times / 2:num_times, plot_temp(num_times/2:num_times, :));
plot_temp = reshape(sd_final(:, 2, :), [num_times red_num_level]);
mean_sd = mean(plot_temp(num_times / 2:num_times, :), 1)
plot(num_times/2:num_times, plot_temp(num_times/2:num_times, :), '--');
grid on;
legend(num2str(mean_rms), num2str(mean_sd));
fprintf(fid, ['V ', num2str(mean_rms), '\n']);
print_file = strcat(dir_name, '/v_ts.eps')
%print(gcf, '-depsc', print_file);
