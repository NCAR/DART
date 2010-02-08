%% stats_wrf_prof_vect

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

mean_ind = ens_size + 1;
std_ind = mean_ind + 1;

itime = input('Time: ');

figure_title = ['Time: ',num2str(itime)];

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 0.9*scrsz(4)],'Name',figure_title)

for field_num = [1:5 7:9]

start_var = 1;
nx = we + 1;
ny = sn;
var_units = 'U (m/s)';
maxlev = bt;
maxval = 8.0;
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
maxval = 0.00004;
end
if field_num > 8
var_units = 'QR (kg/kg)';
maxval = 0.00007;
end

f_size = nx*ny;

y = [1:maxlev];
rmse_prior = y;
rmse_poste = y;
spread_prior = y;
spread_poste = y;

for field_level = 1:maxlev

% Extract field

end_var = start_var + nx*ny - 1;

fname = 'True_State';
state_vec_truth = getnc(fname, 'state',[itime -1 start_var],[itime -1 end_var],[1 1 1]);

fname = 'Prior_Diag';
field_vec_prior = getnc(fname, 'state',[itime mean_ind start_var],[itime mean_ind end_var],[1 1 1]);

field_vec = getnc(fname, 'state',[itime std_ind start_var],[itime std_ind end_var],[1 1 1]);

spread_prior(field_level) = sqrt((field_vec'*field_vec)/(f_size));

fname = 'Posterior_Diag';
field_vec_posterior = getnc(fname, 'state',[itime mean_ind start_var],[itime mean_ind end_var],[1 1 1]);

field_vec = getnc(fname, 'state',[itime std_ind start_var],[itime std_ind end_var],[1 1 1]);

spread_poste(field_level) = sqrt((field_vec'*field_vec)/(f_size));

field_vec = field_vec_prior - state_vec_truth;
rmse_prior(field_level) = sqrt((field_vec'*field_vec)/(f_size));

field_vec = field_vec_posterior - state_vec_truth;
rmse_poste(field_level) = sqrt((field_vec'*field_vec)/(f_size));

start_var = start_var + nx*ny;

end

subplot(3,3,field_num);
plot(rmse_prior,y,'b',rmse_poste,y,'g',spread_prior,y,'--b',spread_poste,y,'--g');

plot_title = [var_units];

title(plot_title);

     axis ([0.0 maxval 0.0 30.0])

end

legend('Prior','Posterior');

% Loop for another try
%map_wrf;

