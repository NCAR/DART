%% map_wrf_diff_time_vect

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

member = input('Input ensemble member: ');

fname = 'Prior_Diag';
tlon = getnc(fname, 'XLON');
we = size(tlon, 2);
tlat = getnc(fname, 'XLAT');
sn = size(tlat, 1);
level = getnc(fname, 'level');
bt = size(level, 1);
true_times = getnc(fname, 'time');
num_true_times = size(true_times, 1)

     stime = input('Initial time : ');
     ftime = input('End time : ');

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
iso = [0.5:1:5];
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
iso = [0.01:0.01:0.1];
end
if field_num > 3
start_var = start_var + (bt + 1)*we*sn;
var_units = 'GZ (m^2/s^2)';
iso = [50:50:300];
end
if field_num > 4
start_var = start_var + (bt + 1)*we*sn;
var_units = 'T (K)';
iso = [0.5:0.5:5];
end
if field_num > 5
start_var = start_var + bt*we*sn;
var_units = 'MU (Pa)';
iso = [100:100:600];
end
if field_num > 6
start_var = start_var + we*sn;
var_units = 'QV (kg/kg)';
iso = [0.0001:0.0001:0.001];
end
if field_num > 7
start_var = start_var + bt*we*sn*(field_num-7);
var_units = 'QC (kg/kg)';
iso = [0.00001:0.00001:0.0001];
end
if field_num > 8
var_units = 'QR (kg/kg)';
iso = [0.00001:0.00001:0.0001];
end

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 0.9*scrsz(4) 0.9*scrsz(4)])

     m = ceil(sqrt(ftime-stime+1));

start_var = start_var + nx*ny*(field_level - 1);
end_var = start_var + nx*ny - 1;

     pane = 1;

for itime = stime:ftime

plot_title = [var_units '   Level: ' num2str(field_level) '   Time: ' num2str(itime)];

% Extract field

fname = 'True_State';
state_vec_truth = getnc(fname, 'state',[itime -1 start_var],[itime -1 end_var],[1 1 1]);

fname = 'Prior_Diag';
state_vec_prior = getnc(fname, 'state',[itime member start_var],[itime member end_var],[1 1 1]);

fname = 'Posterior_Diag';
%fname = 'True_Bdy';
state_vec_posterior = getnc(fname, 'state',[itime member start_var],[itime member end_var],[1 1 1]);

%field_vec = state_vec_prior - state_vec_truth;
field_vec = state_vec_posterior - state_vec_prior;
%field_vec = state_vec_posterior;

field = reshape(field_vec, [nx, ny]);

% Plot field

subplot(m,m,pane);

%nc=5

%colormap = (prism(nc))
if field_num > 2
  [C,h] = contour(tlon,tlat, field', iso );
  hold on
  [Cm,hm] = contour (tlon,tlat, field', -iso, ':');
else
  %[C, h] = contourf(field');
  [C,h] = contour (field', iso);
  hold on
  [Cm,hm] = contour (field', -iso, '--');
end
title(plot_title)
%colorbar('vert')
clabel(C, h);
clabel(Cm, hm);

pane = pane + 1;

end

% Loop for another try
%map_wrf;
