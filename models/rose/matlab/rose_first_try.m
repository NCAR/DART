% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$
 
%Load the true file%
%fname = input('Input ens name');
fname = '../work/Prior_Diag.nc';
lon = getnc(fname, 'lon');
num_lon = size(lon, 1);
lat = getnc(fname, 'lat');
num_lat = size(lat, 1);
level = getnc(fname, 'lev');
num_level = size(level, 1);

%field_name = input('Input field type, 1=u1, 2=v1, 3=t1, 4=u, 5=v, 6=t, 7=qnH, 8=qnOH, 9=qnO')
field_name = 't';
state_vec = getnc(fname, field_name);
time_ind = input('Input time level');
field_level = input('Input level');
ens_ind = input('Input ensemble #');

field = reshape(state_vec( time_ind, ens_ind, :, :, field_level), [num_lat,num_lon]); %truth[lat, lon]
clear state_vec;

figure(1);
[C, h] = contourf(lon,lat,field, 10);
%clabel(C, h);
h = colorbar
set(h, 'Fontsize', 16);
title('ENS', 'fontsize', 16);

% Loop for another try
rose_first_try2;

