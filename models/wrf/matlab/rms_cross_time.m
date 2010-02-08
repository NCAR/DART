%% rms_cross_time

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

fname = 'Prior_Diag';
tlon = getnc(fname, 'XLON');
we = size(tlon, 2);
tlat = getnc(fname, 'XLAT');
sn = size(tlat, 1);
level = getnc(fname, 'level');
bt = size(level, 1);
copy = getnc(fname, 'copy');
ens_size = size(copy, 1) - 2;

true_times = getnc(fname, 'time');
num_true_times = size(true_times, 1);
%num_true_times = 27

dx = 200.;

mean_ind = ens_size + 1;
std_ind = mean_ind + 1;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)])

for field_num = [1:9]

start_var = 1;
nx = we + 1;
ny = sn;
var_units = 'U (m/s)';
maxlev = bt;
maxval = 8.0;
if field_num > 1
start_var = start_var + bt*(we + 1)*sn ;
nx = we;
ny = sn + 1;
var_units = 'V (m/s)';
end
if field_num > 2
start_var = start_var + bt*we*(sn + 1);
nx = we;
ny = sn;
var_units = 'W (m/s)';
maxlev = bt + 1;
maxval = 0.03;
end
if field_num > 3
start_var = start_var + (bt + 1)*we*sn;
var_units = 'GZ (m^2/s^2)';
maxval = 700.0;
end
if field_num > 4
start_var = start_var + (bt + 1)*we*sn;
var_units = 'T (K)';
maxlev = bt;
maxval = 4.0;
end
if field_num > 5
start_var = start_var + bt*we*sn;
var_units = 'MU (Pa)';
maxlev = 1;
maxval = 700.0;
end
if field_num > 6
start_var = start_var + we*sn;
var_units = 'QV (kg/kg)';
maxlev = bt;
maxval = 0.001;
end
if field_num > 7
start_var = start_var + bt*we*sn*(field_num-7);
var_units = 'QC (kg/kg)';
maxval = 0.00003;
end
if field_num > 8
var_units = 'QR (kg/kg)';
maxval = 0.00007;
end

x = [1:nx];
rmse = zeros(nx,num_true_times);
field1d = zeros(maxlev,1);

f_size = nx*ny*maxlev;

end_var = start_var + f_size - 1;

for itime = num_true_times-1:num_true_times

% Extract field

fname = 'True_State';
field_vec_truth = getnc(fname, 'state',[itime -1 start_var],[itime -1 end_var],[1 1 1]);

fname = 'Prior_Diag';
field_vec_prior = getnc(fname, 'state',[itime mean_ind start_var],[itime mean_ind end_var],[1 1 1]);

fname = 'Posterior_Diag';
field_vec_posterior = getnc(fname, 'state',[itime mean_ind start_var],[itime mean_ind end_var],[1 1 1]);

field_vec = field_vec_posterior - field_vec_truth;

field3d = reshape(field_vec, [nx, ny, maxlev]);

for ix = 1:nx
field1d(:) = field3d(ix,sn/2,:);
rmse(ix,itime) = sqrt((field1d'*field1d)/maxlev);
x(ix) = (ix-1)*dx;
end

end

subplot(3,3,field_num);

plot(x,rmse)

     plot_title = [var_units];

     title(plot_title)

%     room = (x(2*num_true_times)-x(1))/10;

%     axis ([(x(1)-room) (x(2*num_true_times)+room) 0.0 maxval])

end

xlabel('distance (km)')

% Loop for another try
%map_wrf;

