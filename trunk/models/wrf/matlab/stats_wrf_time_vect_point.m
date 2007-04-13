% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
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

state_vec_prior = getnc(fname, 'state');

fname = 'Posterior_Diag';

state_vec_posterior = getnc(fname, 'state');

fname = 'True_State';

state_vec_truth = getnc(fname, 'state');

dim = size (state_vec_truth);

len = dim(2);

true_times = getnc(fname, 'time');
num_true_times = size(true_times, 1);

state_vec_truth = reshape(state_vec_truth, [num_true_times, 1, len]);

mean_ind = ens_size + 1;
std_ind = mean_ind + 1;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)])

start_var = 1;
start_var = start_var + bt*(we + 1)*sn ;
start_var = start_var + bt*we*(sn + 1);
%start_var = start_var + (bt + 1)*we*sn;
%start_var = start_var + (bt + 1)*we*sn;
%start_var = start_var + bt*we*sn;
%start_var = start_var + we*sn;
%start_var = start_var + bt*we*sn;
nx = we;
ny = sn;
maxlev = bt+1;
[maxqr,point] = max(state_vec_truth(num_true_times,1, start_var : start_var + nx*ny*maxlev - 1))

for field_num = [1:9]

start_var = 1;
nx = we + 1;
ny = sn;
var_units = 'U (m/s)';
maxlev = bt;
maxval = 7.0;
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
maxval = 0.02;
end
if field_num > 3
start_var = start_var + (bt + 1)*we*sn;
var_units = 'GZ (m^2/s^2)';
maxval = 500.0;
end
if field_num > 4
start_var = start_var + (bt + 1)*we*sn;
var_units = 'T (K)';
maxlev = bt;
maxval = 3.0;
end
if field_num > 5
start_var = start_var + bt*we*sn;
var_units = 'MU (Pa)';
maxlev = 1;
maxval = 600.0;
end
if field_num > 6
start_var = start_var + we*sn;
var_units = 'QV (kg/kg)';
maxlev = bt;
maxval = 0.0007;
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

x = [1:2*num_true_times];
rmse = x;
rmse_2 = x;
spread = x;
spread_2 = x;

start_var = start_var + point-1;

rms_mem = zeros(2*num_true_times,ens_size);

for itime = 1:num_true_times

% Extract field

field_vec = state_vec_prior(itime,mean_ind, start_var) - state_vec_truth(itime,1, start_var);
rmse(2*itime-1) = sqrt((field_vec*field_vec'));

spread(2*itime-1) = state_vec_prior(itime,std_ind, start_var);

x(2*itime-1) = true_times(itime,1)*24; % time in hours

field_vec = state_vec_posterior(itime,mean_ind, start_var) - state_vec_truth(itime,1, start_var);
rmse(2*itime) = sqrt((field_vec*field_vec'));

spread(2*itime) = state_vec_posterior(itime,std_ind, start_var);

x(2*itime) = true_times(itime,1)*24; % time in hours

mean_2 = state_vec_prior(itime,1, start_var);
mean_2 = mean_2 - mean_2;
spread_f = mean_2;

  for imem = 1:ens_size

mean_2 = mean_2 + state_vec_prior(itime,imem, start_var);

    field_vec = state_vec_prior(itime,imem, start_var) - state_vec_truth(itime,1, start_var);
rms_mem(2*itime-1,imem) = sqrt((field_vec*field_vec'));
    
    field_vec = state_vec_posterior(itime,imem, start_var) - state_vec_truth(itime,1, start_var);
rms_mem(2*itime,imem) = sqrt((field_vec*field_vec'));

  end

mean_2 = mean_2/ens_size;

field_vec = mean_2 - state_vec_truth(itime,1, start_var);
rmse_2(2*itime-1) = sqrt((field_vec*field_vec'));
rmse_2(2*itime) = sqrt((field_vec*field_vec'));

  for imem = 1:ens_size

spread_f = spread_f + (state_vec_prior(itime,imem, start_var) - mean_2).^2;

  end

spread_f = sqrt(spread_f/(ens_size-1));
spread_2(2*itime-1) = mean(spread_f);
spread_2(2*itime) = mean(spread_f);

end

subplot(3,3,field_num);
plot(x,rmse,x,spread,x,rms_mem,'-r',x,rmse_2,'--y',x,spread_2,'--m')

     plot_title = [var_units];

     title(plot_title)

%     xlabel('hours')

     room = (x(2*num_true_times)-x(1))/10;

%     axis ([(x(1)-room) (x(2*num_true_times)+room) min(min(rmse,spread)) max(max(rmse,spread))])
     axis ([(x(1)-room) (x(2*num_true_times)+room) 0.0 maxval])

end

legend('RMS error','Spread')
xlabel('hours')

% Loop for another try
%map_wrf;

