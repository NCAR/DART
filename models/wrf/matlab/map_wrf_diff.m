%% map_wrf_diff

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Select field to plot (U, V, W, GZ, T, MU, QV, QC, QR)

field_num = input('Input field type, 1=U, 2=V, 3=W, 4=GZ, 5=T, 6=MU, 7=QV, 8=QC, 9=QR: ');

map_proj = {'lambert', 'ups', 'mercator'};

member = input('Input ensemble member: ');

itime = input('Time: ');

fname = 'Prior_Diag';
tlon = getnc(fname, 'XLON');
we = size(tlon, 2);
tlat = getnc(fname, 'XLAT');
sn = size(tlat, 1);
level = getnc(fname, 'level');
bt = size(level, 1);

% Get level for free atmosphere fields
if field_num == 6
   field_level = 1;
else
   field_level = input('Input level: ');
end

start_var = 1;
nx = we + 1;
ny = sn;
var_units = 'U (m/s)';
if field_num > 1
start_var = start_var + bt*(we + 1)*sn;
nx = we;
ny = sn + 1;
var_units = 'V (m/s)';
end
if field_num > 2
start_var = start_var + bt*we*(sn + 1);
nx = we;
ny = sn;
var_units = 'W (m/s)';
end
if field_num > 3
start_var = start_var + (bt + 1)*we*sn;
var_units = 'GZ (m^2/s^2)';
end
if field_num > 4
start_var = start_var + (bt + 1)*we*sn;
var_units = 'T (K)';
end
if field_num > 5
start_var = start_var + bt*we*sn;
var_units = 'MU (Pa)';
end
if field_num > 6
start_var = start_var + we*sn;
var_units = 'QV (kg/kg)';
end
if field_num > 7
start_var = start_var + bt*we*sn*(field_num-7);
var_units = 'QC (kg/kg)';
end
if field_num > 8
var_units = 'QR (kg/kg)';
end

plot_title = [var_units '   Level: ' num2str(field_level) '   Time: ' num2str(itime)];

start_var = start_var + nx*ny*(field_level - 1);
end_var = start_var + nx*ny - 1;

% Extract field

fname = 'True_State';
state_vec_truth = getnc(fname, 'state',[itime -1 start_var],[itime -1 end_var],[1 1 1]);

fname = 'Prior_Diag';
state_vec_prior = getnc(fname, 'state',[itime member start_var],[itime member end_var],[1 1 1]);

fname = 'Posterior_Diag';
state_vec_posterior = getnc(fname, 'state',[itime member start_var],[itime member end_var],[1 1 1]);

field_vec = state_vec_posterior - state_vec_truth;
%field_vec = state_vec_posterior - state_vec_prior;
%field_vec = state_vec_posterior;

field = reshape(field_vec, [nx, ny]);

% Plot field

%nc=5

%colormap = (prism(nc))
if field_num > 2
[C, h] = contourf(tlon,tlat,field');
else
%[C, h] = contourf(field');
[C,h] = contour ( field' , [0.5:1:5] );
hold on
[Cm,hm] = contour ( field' ,- [0.5:1:5] , 'k:');
end
title(plot_title)
colorbar('vert')
clabel(C, h);

% Loop for another try
%map_wrf;
