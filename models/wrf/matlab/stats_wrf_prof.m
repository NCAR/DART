fname = 'Prior_Diag';
tlon = getnc(fname, 'west_east');
we = size(tlon, 1);
tlat = getnc(fname, 'south_north');
sn = size(tlat, 1);
level = getnc(fname, 'bottom_top');
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

itime = input('Time: ');

poste_tru = state_vec_posterior(itime,mean_ind,:) - state_vec_truth(itime,1,:);
prior_tru = state_vec_prior(itime,mean_ind, :) - state_vec_truth(itime,1,:);

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)])

for field_num = [1:5 7:9]

start_var = 1
nx = we + 1
ny = sn
var_units = 'U (m/s)'
maxlev = bt
if field_num > 1
   start_var = start_var + bt*(we + 1)*sn
   nx = we
   ny = sn + 1
var_units = 'V (m/s)'
end
if field_num > 2
   start_var = start_var + bt*we*(sn + 1)
   nx = we
   ny = sn
var_units = 'W (m/s)'
maxlev = bt + 1
end
if field_num > 3
   start_var = start_var + (bt + 1)*we*sn
var_units = 'GZ (m^2/s^2)'
end
if field_num > 4
   start_var = start_var + (bt + 1)*we*sn
var_units = 'T (K)'
maxlev = bt
end
if field_num > 5
   start_var = start_var + bt*we*sn
var_units = 'MU (Pa)'
maxlev = 1
end
if field_num > 6
   start_var = start_var + we*sn
var_units = 'QV (kg/kg)'
maxlev = bt
end
if field_num > 7
   start_var = start_var + bt*we*sn*(field_num-7)
var_units = 'QC (kg/kg)'
end
if field_num > 8
var_units = 'QR (kg/kg)'
end

y = [1:maxlev];
rmse_prior = y;
rmse_poste = y;
spread_prior = y;
spread_poste = y;

for field_level = 1:maxlev

% Extract field

field_vec = prior_tru(start_var : start_var + nx*ny - 1);
rmse_prior(field_level) = sqrt((field_vec*field_vec')/(nx*ny));

spread_prior(field_level) = mean(state_vec_prior(itime,std_ind, start_var : start_var + nx*ny - 1));

field_vec = poste_tru(start_var : start_var + nx*ny - 1);
rmse_poste(field_level) = sqrt((field_vec*field_vec')/(nx*ny));

spread_poste(field_level) = mean(state_vec_posterior(itime,std_ind, start_var : start_var + nx*ny - 1));

start_var = start_var + nx*ny;

end

subplot(3,3,field_num);
plot(rmse_prior,y,'b',rmse_poste,y,'g',spread_prior,y,'--b',spread_poste,y,'--g');

plot_title = [var_units];

title(plot_title);

end

legend('Prior','Posterior');

% Loop for another try
%map_wrf;

