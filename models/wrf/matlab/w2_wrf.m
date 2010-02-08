%% w2_wrf

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
tlon = getnc(fname, 'west_east');
we = size(tlon, 2);
tlat = getnc(fname, 'south_north');
sn = size(tlat, 1);
level = getnc(fname, 'bottom_top');
bt = size(level, 1);
copy = getnc(fname, 'copy');
ens_size = size(copy, 1) - 2;

true_times = getnc(fname, 'time');
num_true_times = size(true_times, 1);

mean_ind = ens_size + 1;
std_ind = mean_ind + 1;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)])

     field_num = 3;

start_var = 1;
start_var = start_var + bt*(we + 1)*sn ;
start_var = start_var + bt*we*(sn + 1);
nx = we;
ny = sn;
var_units = 'W (m/s)';
maxlev = bt + 1;
maxval = 0.03;

x = [1:2*num_true_times];
estim = x;
truth = x;

rms_mem = zeros(2*num_true_times,ens_size);

f_size = nx*ny*maxlev;

end_var = start_var + f_size - 1;

for itime = 1:num_true_times

% Extract field

fname = 'True_State';
field_vec_truth = getnc(fname, 'state',[itime -1 start_var],[itime -1 end_var],[1 1 1]);

truth(2*itime-1) = sqrt((field_vec_truth'*field_vec_truth)/(f_size));
truth(2*itime) = truth(2*itime-1);

fname = 'Prior_Diag';
field_vec_prior = getnc(fname, 'state',[itime mean_ind start_var],[itime mean_ind end_var],[1 1 1]);

estim(2*itime-1) = sqrt((field_vec_prior'*field_vec_prior)/(f_size));

x(2*itime-1) = true_times(itime,1)*24; % time in hours

fname = 'Posterior_Diag';
field_vec_posterior = getnc(fname, 'state',[itime mean_ind start_var],[itime mean_ind end_var],[1 1 1]);

estim(2*itime) = sqrt((field_vec_posterior'*field_vec_posterior)/(f_size));

x(2*itime) = true_times(itime,1)*24; % time in hours

  for imem = 1:ens_size

fname = 'Prior_Diag';
field_vec_prior = getnc(fname, 'state',[itime imem start_var],[itime imem end_var],[1 1 1]);
field_vec = field_vec_prior - field_vec_truth;
rms_mem(2*itime-1,imem) = sqrt((field_vec'*field_vec)/(f_size));
    
fname = 'Posterior_Diag';
field_vec_posterior = getnc(fname, 'state',[itime imem start_var],[itime imem end_var],[1 1 1]);
field_vec = field_vec_posterior - field_vec_truth;
rms_mem(2*itime,imem) = sqrt((field_vec'*field_vec)/(f_size));

E2(2*itime-1) = E2(2*itime-1) + rms_mem(2*itime-1,imem);
E2(2*itime) = E2(2*itime) + rms_mem(2*itime,imem);

  end

end


%subplot(3,3,field_num);

if (ens_size > 0.0)
     E2 = E2/ens_size;
     plot(x,truth,x,estim,x,E2,'--m')
else
     plot(x,truth,x,estim)
end

     plot_title = [var_units];

     title(plot_title)

%     xlabel('hours')

     room = (x(2*num_true_times)-x(1))/10;

%     axis ([(x(1)-room) (x(2*num_true_times)+room) min(min(truth,estim)) max(max(truth,estim))])
     axis ([(x(1)-room) (x(2*num_true_times)+room) 0.0 maxval])

legend('Truth','State Estimate')
xlabel('hours')

% Loop for another try
%map_wrf;

