%% map_wrf_spread_time

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

fname = 'Prior_Diag';
tlon = getnc(fname, 'XLON');
we = size(tlon, 2);
tlat = getnc(fname, 'XLAT');
sn = size(tlat, 1);
level = getnc(fname, 'level');
bt = size(level, 1);
ncopy = getnc(fname, 'copy');
ens_size = size(ncopy, 1) - 2;

mean_ind = ens_size + 1;
std_ind = mean_ind + 1;

icopy = std_ind;

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

nx = we + 1;
ny = sn;
var_units = 'U (m/s)';
var_name = 'U';
maxlev = bt;
iso = [0.5:1:5];
if field_num > 1
nx = we;
ny = sn + 1;
var_units = 'V (m/s)';
var_name = 'V';
end
if field_num > 2
nx = we;
ny = sn;
var_units = 'W (m/s)';
var_name = 'W';
maxlev = bt + 1;
iso = [0.01:0.01:0.1];
end
if field_num > 3
var_units = 'GZ (m^2/s^2)';
var_name = 'PH';
iso = [50:50:300];
end
if field_num > 4
var_units = 'T (K)';
var_name = 'T';
maxlev = bt;
iso = [0.5:0.5:5];
end
if field_num > 5
var_units = 'MU (Pa)';
var_name = 'MU';
maxlev = 1;
iso = [100:100:600];
end
if field_num > 6
var_units = 'QV (kg/kg)';
var_name = 'QVAPOR';
maxlev = bt;
iso = [0.0001:0.0001:0.001];
end
if field_num > 7
var_units = 'QC (kg/kg)';
var_name = 'QCLOUD';
iso = [0.00001:0.00001:0.0001];
end
if field_num > 8
var_units = 'QR (kg/kg)';
var_name = 'QRAIN';
iso = [0.00001:0.00001:0.0001];
end

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 0.9*scrsz(4) 0.9*scrsz(4)])

     m = ftime-stime+1;

     pane = 1;

for itime = stime:ftime

plot_title = [var_units '   Level: ' num2str(field_level) '   Time: ' num2str(itime)];

if maxlev > 1
corner_m = [itime icopy field_level -1 -1];
end_point_m = [itime icopy field_level -1 -1];
stride = [1 1 1 1 1];
else
corner_m = [itime icopy -1 -1];
end_point_m = [itime icopy -1 -1];
stride = [1 1 1 1];
end

% Extract field

fname = 'Prior_Diag';
state_vec_prior = getnc(fname, var_name,corner_m,end_point_m,stride);

field = state_vec_prior;

% Plot field

subplot(m,2,pane);

%nc=5

%colormap = (prism(nc))
[C, h] = contourf(field);
%[C,h] = contour (field, iso);
hold on
%[Cm,hm] = contour (field, -iso, '--');

title(plot_title)
%colorbar('vert')
%clabel(C, h);
%clabel(Cm, hm);

pane = pane + 1;

fname = 'Posterior_Diag';
state_vec_posterior = getnc(fname, var_name,corner_m,end_point_m,stride);

field = state_vec_posterior;

% Plot field

subplot(m,2,pane);

%nc=5

%colormap = (prism(nc))
[C, h] = contourf(field);
%[C,h] = contour (field, iso);
hold on
%[Cm,hm] = contour (field, -iso, '--');

title(plot_title)
%colorbar('vert')
%clabel(C, h);
%clabel(Cm, hm);

pane = pane + 1;

end

% Loop for another try
%map_wrf_diff_time;
