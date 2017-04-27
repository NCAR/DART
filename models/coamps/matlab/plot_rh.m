%% plot_rh

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

clear all;
close all;

xdim = 100;
ydim = 100;

% Open up the various NetCDF files
%truth = netcdf('true_state.nc');
prior = netcdf('preassim.nc');
postr = netcdf('analysis.nc');

% Only consider a single state variable - note that the first two
% items are the ensemble mean and spread
avail_times = 1;
ijarea = xdim*ydim;
total_members = 100;

% Calculate the vector position for this particular state variable
% field
% THETA
varnum = 160;
spots = ((varnum-1)*ijarea+1):(varnum*ijarea);

% Read in the true value, the prior and posterior mean/spread
theta_prior_mean = squeeze(prior{'state'}(avail_times,1,spots));
theta_postr_mean = squeeze(prior{'state'}(avail_times,1,spots));

% Have to iterate to get all the ensemble members - offset by 2
% since the mean and spread take up the first slpts
for ii=[1:total_members]
  theta_prior_state(:,:,ii) = squeeze(prior{'state'}(avail_times,ii+2,spots));
  theta_postr_state(:,:,ii) = squeeze(postr{'state'}(avail_times,ii+2,spots));
end

% PERT. EXNER
varnum = 250;
spots = ((varnum-1)*ijarea+1):(varnum*ijarea);

% Read in the true value, the prior and posterior mean/spread
ppert_prior_mean = squeeze(prior{'state'}(avail_times,1,spots));
ppert_postr_mean = squeeze(prior{'state'}(avail_times,1,spots));

% Have to iterate to get all the ensemble members - offset by 2
% since the mean and spread take up the first slpts
for ii=[1:total_members]
  ppert_prior_state(:,:,ii) = squeeze(prior{'state'}(avail_times,ii+2,spots));
  ppert_postr_state(:,:,ii) = squeeze(postr{'state'}(avail_times,ii+2,spots));
end

% MEAN EXNER
varnum = 295;
spots = ((varnum-1)*ijarea+1):(varnum*ijarea);

% Read in the true value, the prior and posterior mean/spread
pmean_prior_mean = squeeze(prior{'state'}(avail_times,1,spots));
pmean_postr_mean = squeeze(prior{'state'}(avail_times,1,spots));

% Have to iterate to get all the ensemble members - offset by 2
% since the mean and spread take up the first slpts
for ii=[1:total_members]
  pmean_prior_state(:,:,ii) = squeeze(prior{'state'}(avail_times,ii+2,spots));
  pmean_postr_state(:,:,ii) = squeeze(postr{'state'}(avail_times,ii+2,spots));
end

% VAPOR
varnum = 386;
spots = ((varnum-1)*ijarea+1):(varnum*ijarea);

% Read in the true value, the prior and posterior mean/spread
qvapr_prior_mean = squeeze(prior{'state'}(avail_times,1,spots));
qvapr_postr_mean = squeeze(prior{'state'}(avail_times,1,spots));

% Have to iterate to get all the ensemble members - offset by 2
% since the mean and spread take up the first slpts
for ii=[1:total_members]
  qvapr_prior_state(:,:,ii) = squeeze(prior{'state'}(avail_times,ii+2,spots));
  qvapr_postr_state(:,:,ii) = squeeze(postr{'state'}(avail_times,ii+2,spots));
end

prior_mean = calculate_rh(theta_prior_mean, ppert_prior_mean, ...
                          pmean_prior_mean, qvapr_prior_mean);
postr_mean = calculate_rh(theta_postr_mean, ppert_postr_mean, ...
                          pmean_postr_mean, qvapr_postr_mean);
prior_state = calculate_rh(theta_prior_state, ppert_prior_state, ...
                          pmean_prior_state, qvapr_prior_state);
postr_state = calculate_rh(theta_postr_state, ppert_postr_state, ...
                          pmean_postr_state, qvapr_postr_state);

prior_spread = std(squeeze(prior_state),0,2);
postr_spread = std(squeeze(postr_state),0,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lat/Lon information comes from the truth run
lat = squeeze(prior{'lat'}(spots));
lon = squeeze(prior{'lon'}(spots));

% Convert the vector to something we can easily plot
% NEED TO DO THIS AS [ydim xdim] when using lat/lon!!!
%truth_state_field = reshape(truth_state,[ydim xdim]);
for ii=1:total_members
  prior_state_field(:,:,ii) = reshape(prior_state(:,:,ii),[ydim xdim]);
  postr_state_field(:,:,ii) = reshape(postr_state(:,:,ii),[ydim xdim]);
end
prior_spread_field = reshape(prior_spread,[ydim xdim]);
postr_spread_field = reshape(postr_spread,[ydim xdim]);
prior_mean_field = reshape(prior_mean,[ydim xdim]);
postr_mean_field = reshape(postr_mean,[ydim xdim]);

lat_field = reshape(lat,[ydim xdim]);
lon_field = reshape(lon,[ydim xdim]);

% Observation information
obs_spacing = 4;
lat_obs = lat(1:obs_spacing:end);
lon_obs = lon(1:obs_spacing:end);

% Standard contour levels
n_contours = 10;
%stdlevels =
%linspace(min(truth_state),max(truth_state),n_contours);
stdlevels = linspace(min(prior_mean),max(prior_mean),n_contours);

% Plot the ensemble spread, mean, and truth
disp('Creating prior plot...')
priorfig = figure('Visible','on');
contourf(lon_field,lat_field,prior_spread_field);
cur_ax = gca;
caxis(cur_ax,[min(min(prior_spread_field)) max(max(prior_spread_field))]);
shading flat;
alpha(0.3);
colorbar('h');
hold on;
% for ii=1:total_members
%   [c_mem,h_mem] = contour(lon_field, lat_field, prior_state_field(: ...
% 						  ,:,ii),stdlevels); 
%   set(h_mem, 'LineWidth',1.0,'LineStyle','-','LineColor',[1 0 0])
% end
[c_truth,h_truth] = contour(lon_field, lat_field, ...
			    prior_mean_field,stdlevels, 'LineColor', ...
			    [0 0 0]); 
clabel(c_truth,h_truth);
axis image;
set(h_truth,'LineWidth',2.0,'LineStyle','-','LineColor',[0 0 0]);
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees N)');
title (['Prior Estimate of Ensemble Mean and Spread: 31 August 1999' ...
	' : \tau = ' sprintf('%dhrs',avail_times + 11)]);
% leg = legend([h_mem(1) h_truth(1)],'Ensemble Members','Ensemble Mean', ...
%              'Location','SouthOutside');
%saveas(priorfig,'prior','png');

% Plot the ensemble spread, mean, and truth
disp('Creating posterior plot...')
postrfig = figure('Visible','on');
contourf(lon_field,lat_field,postr_spread_field);
cur_ax = gca;
caxis(cur_ax,[min(min(postr_spread_field)) max(max(postr_spread_field))]);
shading flat;
alpha(0.3);
colorbar('h');
hold on;
% for ii=1:total_members
%   [c_mem,h_mem] = contour(lon_field, lat_field, postr_state_field(: ...
% 						  ,:,ii),stdlevels); 
%   set(h_mem, 'LineWidth',0.0,'LineStyle','-','LineColor',[1 0 0])
% end
[c_truth,h_truth] = contour(lon_field, lat_field, ...
			    postr_mean_field,stdlevels, 'LineColor', ...
			    [0 0 0]);
clabel(c_truth,h_truth);
axis image;
%set(h_truth,'LineWidth',2.0,'LineStyle','-','LineColor',[0 1 0]);
set(h_truth,'LineWidth',2.0,'LineStyle','-','LineColor',[0 0 0]);
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees N)');
title (['Posterior Estimate of Ensemble Mean and Spread: 31 August 1999 : \tau = ' sprintf('%dhrs',avail_times + 11)]);
% leg = legend([h_mem(1) h_truth(1)],'Ensemble Members','Truth', ...
%              'Location','SouthOutside');
%saveas(postrfig,'posterior','png');

% Plot the difference between prior and posterior
disp('Creating difference plot...')
difffig = figure('visible','on');
disp('Calculating difference field')
diff_field = sqrt(sum((postr_state_field-prior_state_field).^2,3))/100;
disp('Plotting difference field')
pcolor(lon_field,lat_field,diff_field);
disp('Trying to shade...')
shading flat;
colorbar('h')
hold on;
axis image;
plot(lon_obs,lat_obs,'xk');
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees N)');
title (['Difference between Prior and Posterior Estimates: ' ...
        '31 August 1999 : \tau = ' sprintf('%dhrs',avail_times ...
                                           + 11)]);
%saveas(difffig,'difference','png');

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
