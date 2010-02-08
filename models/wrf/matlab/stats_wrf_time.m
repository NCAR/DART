%% stats_wrf_time

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

fname = 'Posterior_Diag';

nc = netcdf( [fname,'.nc'] , 'nowrite' ) ;

num_domains = size(nc('domain'),1);

if (num_domains > 1)

   disp(['Number of domains: ',int2str(num_domains)])
   id = input('Input domain id: ');

else

   id = 1;

end

sls = size(nc(['soil_layers_stag_d0',int2str(id)]),1);
we = size(nc(['west_east_d0',int2str(id)]), 1);
sn = size(nc(['south_north_d0',int2str(id)]), 1);
bt = size(nc(['bottom_top_d0',int2str(id)]), 1);
ens_size = size(nc('copy'), 1) - 2;
close(nc);

true_times = getnc(fname, 'time');
num_true_times = size(true_times, 1);
true_times = true_times - true_times(1);
time_unit = 'days';
if (true_times(num_true_times) < 1.0)
     true_times = true_times*24;
     time_unit = 'hours';
end
if (true_times(num_true_times) < 1.0)
     true_times = true_times*60;
     time_unit = 'minutes';
end
if (true_times(num_true_times) < 1.0)
     true_times = true_times*60;
     time_unit = 'seconds';
end

mean_ind = ens_size + 1;
std_ind = mean_ind + 1;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 0.7*scrsz(3) 0.7*scrsz(4)])

axes('FontSize',12)

for field_num = [1:9]

nx = we + 1;
ny = sn;
var_units = 'U (m/s)';
var_name = ['U_d0',int2str(id)];
maxlev = bt;
%maxval = 8.0;
maxval = 2.0;
if field_num > 1
nx = we;
ny = sn;
var_units = 'W (m/s)';
var_name = ['W_d0',int2str(id)];
maxlev = bt + 1;
%maxval = 0.03;
maxval = 0.3;
end
if field_num > 2
var_units = 'PH (m^2/s^2)';
var_name = ['PH_d0',int2str(id)];
%maxval = 700.0;
maxval = 100.0;
end
if field_num > 3
var_units = 'T (K)';
var_name = ['T_d0',int2str(id)];
maxlev = bt;
%maxval = 4.0;
maxval = 1.0;
end
if field_num > 4
var_units = 'MU (Pa)';
var_name = ['MU_d0',int2str(id)];
maxlev = 1;
%maxval = 700.0;
maxval = 100.0;
end
if field_num > 5
var_units = 'TSLB (K)';
var_name = ['TSLB_d0',int2str(id)];
maxlev = sls;
%maxval = 4.0;
maxval = 0.1;
end
if field_num > 6
var_units = 'QV (kg/kg)';
var_name = ['QVAPOR_d0',int2str(id)];
maxlev = bt;
%maxval = 0.001;
maxval = 0.0005;
end
if field_num > 7
var_units = 'QC (kg/kg)';
var_name = ['QCLOUD_d0',int2str(id)];
maxval = 0.00007;
end
if field_num > 8
var_units = 'QR (kg/kg)';
var_name = ['QRAIN_d0',int2str(id)];
%maxval = 0.00007;
maxval = 0.0002;
end

x = [1:2*num_true_times];
rmse = x;
spread = x;
E2 = x-x;

rms_mem = zeros(2*num_true_times,ens_size);

f_size = nx*ny*maxlev;

for itime = 1:num_true_times

if maxlev > 1
corner_m = [itime mean_ind -1 -1 -1];
end_point_m = [itime mean_ind -1 -1 -1];
stride = [1 1 1 1 1];
corner_s = [itime std_ind -1 -1 -1];
end_point_s = [itime std_ind -1 -1 -1];
corner_t = [itime -1 -1 -1 -1];
end_point_t = [itime -1 -1 -1 -1];
else
corner_m = [itime mean_ind -1 -1];
end_point_m = [itime mean_ind -1 -1];
stride = [1 1 1 1];
corner_s = [itime std_ind -1 -1];
end_point_s = [itime std_ind -1 -1];
corner_t = [itime -1 -1 -1];
end_point_t = [itime -1 -1 -1];
end

% Extract field

fname = 'True_State';
field_vec_truth = reshape(getnc(fname, var_name,corner_t,end_point_t,stride),f_size,1);

fname = 'Prior_Diag';

field_vec_prior = reshape(getnc(fname, var_name,corner_m,end_point_m,stride),f_size,1);

field_vec = reshape(getnc(fname, var_name,corner_s,end_point_s,stride),f_size,1);

spread(2*itime-1) = sqrt((field_vec'*field_vec)/(f_size));

field_vec = field_vec_prior - field_vec_truth;
rmse(2*itime-1) = sqrt((field_vec'*field_vec)/(f_size));

x(2*itime-1) = true_times(itime);

fname = 'Posterior_Diag';
field_vec_posterior = reshape(getnc(fname, var_name,corner_m,end_point_m,stride),f_size,1);

field_vec = reshape(getnc(fname, var_name,corner_s,end_point_s,stride),f_size,1);

spread(2*itime) = sqrt((field_vec'*field_vec)/(f_size));

field_vec = field_vec_posterior - field_vec_truth;
rmse(2*itime) = sqrt((field_vec'*field_vec)/(f_size));

x(2*itime) = true_times(itime);

  for imem = 1:ens_size

if maxlev > 1
corner = [itime imem -1 -1 -1];
end_point = [itime imem -1 -1 -1];
else
corner = [itime imem -1 -1];
end_point = [itime imem -1 -1];
end

fname = 'Prior_Diag';
field_vec_prior = reshape(getnc(fname, var_name,corner,end_point,stride),f_size,1);
field_vec = field_vec_prior - field_vec_truth;
rms_mem(2*itime-1,imem) = sqrt((field_vec'*field_vec)/(f_size));
    
fname = 'Posterior_Diag';
field_vec_posterior = reshape(getnc(fname, var_name,corner,end_point,stride),f_size,1);
field_vec = field_vec_posterior - field_vec_truth;
rms_mem(2*itime,imem) = sqrt((field_vec'*field_vec)/(f_size));

E2(2*itime-1) = E2(2*itime-1) + rms_mem(2*itime-1,imem);
E2(2*itime) = E2(2*itime) + rms_mem(2*itime,imem);

  end

end

subplot(3,3,field_num);

%if (ens_size > 0.0)
%     E2 = E2/ens_size;
%     plot(x,rmse,x,spread,x,E2,'--m')
%else
     plot(x,rmse,x,spread,'LineWidth',2)
%     plot(x,rmse,'--k',x,rmse_3dvar(:,field_num),':k',x,rmse_no_assim(:,field_num),'-.k','LineWidth',2)
%     plot(x,rmse,x,rmse_OSSE31(:,field_num),x,rmse_OSSE30(:,field_num),x,rmse_OSSE26(:,field_num),x,rmse_OSSE32(:,field_num),x,rmse_OSSE33(:,field_num),x,rmse_OSSE34(:,field_num),'LineWidth',2)
%end

%rmse_OSSE34(:,field_num) = rmse;

     plot_title = [var_units];

     title(plot_title,'Fontsize',12)

%     xlabel(time_unit)

     room = (x(2*num_true_times)-x(1)+1)/10;

%     axis ([(x(1)-room) (x(2*num_true_times)+room) min(min(rmse,spread)) max(max(rmse,spread))])
     axis ([(x(1)-room) (x(2*num_true_times)+room) 0.0 maxval])

%     if (field_num == 6)
%        legend('EnKF','3D-Var','no assim')
%        legend('EnKF','3D-Var')
%     end

end

%legend('Clim BC','20% 3D-Var BC','50% 3D-Var BC','100% 3D-Var BC','120% 3D-Var BC','150% 3D-Var BC','200% 3D-Var BC')

%if (ens_size > 0.0)
%     legend('RMS error','Spread','E2')
%else
     legend('RMS error','Spread')
%end
xlabel(time_unit,'Fontsize',12)

% Loop for another try
%map_wrf;
