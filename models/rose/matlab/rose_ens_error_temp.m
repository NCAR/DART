% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$
 
% Assumes 2 copies of data are ensemble mean and spread
% Should be checked and automated

%Load the true file%
fname = input('Input true state name');
%fname = 'True_State.nc';
lon = getnc(fname, 'lon');
num_lon = size(lon, 1);
lat = getnc(fname, 'lat');
num_lat = size(lat, 1);
level = getnc(fname, 'lev');
num_level = size(level, 1);

%Load the ensemble file%
ens_fname = input('Input file name for ensemble');
%ens_fname = 'Prior_Diag.nc'

% Select field to plot (u, v, t, ...)
field_name = input('Input field type, 1=u1, 2=v1, 3=t1, 4=u, 5=v, 6=t, 7=qnH, 8=qnOH, 9=qnO')
ens_vec = getnc(ens_fname, field_name);
state_vec = getnc(fname, field_name);

% Get a time level from the user%
time_ind = input('Input time level');

% Get level for free atmosphere fields
field_level = input('Input level');

% Extract state and ensemble for just this time
field = reshape(state_vec( time_ind, :, :, field_level), [num_lat,num_lon]); %truth[lat, lon]
clear state_vec;

ens_spread = ens_vec(time_ind, 12, :, :, field_level);

% Get ensemble mean and spread
ens =   reshape(ens_vec(time_ind, 11, :, :, field_level),[num_lat,num_lon]);  %ens_mean[lat, lon]
spread = reshape(ens_vec(time_ind, 12, :, :, field_level),[num_lat,num_lon]); %ens_spread[lat, lon]
clear ens_vec

figure;
subplot(2, 2, 1);
[C, h] = contourf(lon,lat,field, 10);
%clabel(C, h);
h = colorbar
set(h, 'Fontsize', 16);
title('Truth', 'fontsize', 16);

subplot(2, 2, 2);
[C, h] = contourf(lon,lat,ens, 10);
%clabel(C, h);
h = colorbar
set(h, 'Fontsize', 16);
title('Ensemble Mean Analysis', 'fontsize', 16);

% Compute and plot the difference field
ens_err = (ens - field);
subplot(2, 2, 3);
[C, h] = contourf(lon,lat,ens_err, 10);
%clabel(C, h);
h = colorbar
set(h, 'Fontsize', 16);

% Compute statistics of the error field
max_err = max(max(ens_err));
min_err = min(min(ens_err));
rms_err = mean(mean(abs(ens_err)));

% Label figure 3 with these statistics
title_string = ['RMS ERROR = ', num2str(rms_err)];
title (title_string, 'Fontsize', 16)

% Output the spread plot, too
%figure(2);
subplot(2, 2, 4);
[C, h] = contourf(lon,lat,spread, 10);
clabel(C, h);

% Compute statistics of the spread field
max_spread = max(max(ens_spread));
min_spread = min(min(ens_spread));
rms_spread = mean(mean(ens_spread));

% Label figure 4 with these statistics
title_string = ['Min = ', num2str(min_spread), ' Max =  ', num2str(max_spread), '   RMS ERROR = ', num2str(rms_spread)];
title (title_string)



